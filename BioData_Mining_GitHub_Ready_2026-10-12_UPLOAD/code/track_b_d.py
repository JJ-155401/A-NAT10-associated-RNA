"""Track D: Observational mediation Braak -> NAT10 -> ac4C-target.
Also Track B: single-gene directionality for top AD-risk ac4C targets.
All on real AD brains (GSE125583, n=289). No wet lab, public data.
"""
import pandas as pd, numpy as np, gzip, os
from scipy import stats
import statsmodels.api as sm
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('AD_NAT10_DATA_DIR', os.path.abspath(os.path.join(HERE, '..', 'data', 'external')))
OUT = os.path.join(HERE,'..','results'); FIG = os.path.join(HERE,'..','figures')
os.makedirs(OUT, exist_ok=True)

meta = pd.read_csv(os.path.join(ROOT,'processed','GSE125583_metadata.csv')).reset_index(drop=True)
counts = pd.read_parquet(os.path.join(ROOT,'processed','GSE125583_counts.parquet'))
counts.index=[str(x).split('_sample')[0] for x in counts.index]
counts = counts.loc[meta['GSM'].values].reset_index(drop=True)

ac2sym={}; ac2entrez={}
with gzip.open(os.path.join(ROOT, 'reference', 'uniprot_human_idmapping.dat.gz'),'rt') as f:
    for line in f:
        p=line.rstrip('\n').split('\t')
        if len(p)<3: continue
        if p[1]=='Gene_Name': ac2sym[p[0]]=p[2]
        elif p[1]=='GeneID': ac2entrez[p[0]]=p[2]
entrez2sym={ac2entrez[a]:ac2sym[a] for a in ac2entrez if a in ac2sym}
sym2entrez={v:k for k,v in entrez2sym.items()}
ac4c=set(pd.read_csv(os.path.join(HERE, '..', 'data', 'human_ac4c_targets_curated.csv'))['gene_symbol'])

# NAT10 expression
meta['NAT10']=np.log2(counts['55226']+1).values
# Braak ordinal
meta['braak_ord']=meta['braak'].map({'I':1,'II':2,'III':3,'IV':4,'V':5,'VI':6})
meta['age']=pd.to_numeric(meta['age'], errors='coerce')
meta['sex_m']=(meta['sex']=='M').astype(int)
cov=['age','sex_m']

# ---- Track D: mediation Braak -> NAT10 -> target score ----
ac4c_cols=[str(c) for c in counts.columns if str(c) in [sym2entrez[s] for s in ac4c if s in sym2entrez]]
meta['ac4c_score']=np.log2(counts[ac4c_cols]+1).mean(axis=1).values

def mediate(df, X, M, Y, covs):
    d=df[[X,M,Y]+covs].dropna()
    # total effect X->Y
    Xt=sm.add_constant(d[[X]+covs]); yt=d[Y]
    t_mod=sm.OLS(yt, Xt).fit(); c=d[X].corr(d[Y])
    # X->M
    Xm=sm.add_constant(d[[X]+covs]); mt=d[M]
    a_mod=sm.OLS(mt, Xm).fit(); a=a_mod.params[X]
    # M->Y | X
    Ym=sm.add_constant(d[[X,M]+covs]); ym=d[Y]
    b_mod=sm.OLS(ym, Ym).fit(); b=b_mod.params[M]
    # indirect
    indirect=a*b
    # direct X->Y | M
    direct=b_mod.params[X]
    return {'path_X_M(a)':a,'path_M_Y(b)':b,'indirect(axb)':indirect,'direct_cprime':direct,'total_c':t_mod.params[X],'n':len(d)}

med=mediate(meta,'braak_ord','NAT10','ac4c_score',cov)
print('=== Mediation: Braak -> NAT10 -> ac4C-target score ===')
for k,v in med.items(): print(f'  {k}: {v:.4f}' if isinstance(v,float) else f'  {k}: {v}')
pd.DataFrame([med]).to_csv(os.path.join(OUT,'mediation_braak_nat10_ac4c.csv'), index=False)

# ---- Track B: top-target directionality ----
top=['APP','BIN1','CLU','SORL1','TREM2','MAPT','PSEN1','APOE','MEF2C','CD33','ABCA7','PICALM','NAT10']
rows=[]
for s in top:
    if s not in sym2entrez: continue
    col=str(sym2entrez[s])
    if col not in counts.columns: continue
    expr=np.log2(counts[col]+1)
    mask_ad=meta['diagnosis'].values=='AD'; mask_ct=meta['diagnosis'].values=='Control'
    ad=expr[mask_ad]; ct=expr[mask_ct]
    t,p=stats.ttest_ind(ad,ct,equal_var=False)
    # ac4C status
    is_ac4c = s in ac4c
    # Naive ac4C role inference: if target down in AD (where ac4C/NAT10 also down),
    # ac4C likely ACTIVATES it (loss -> down). If target up in AD, ac4C likely REPRESSES it.
    direction='DOWN' if (ad.mean()-ct.mean())<0 else 'UP'
    ac4c_role = 'activator (loss->target down)' if direction=='DOWN' else 'repressor (loss->target up)'
    rows.append({'symbol':s,'entrez':sym2entrez[s],'ac4C_target':is_ac4c,'AD_log2FC':round(float(ad.mean()-ct.mean()),3),
                 'p':round(float(p),4),'AD_direction':direction,'inferred_ac4C_role':ac4c_role})
tb=pd.DataFrame(rows)
tb.to_csv(os.path.join(OUT,'top_targets_directionality.csv'), index=False)
print('\n=== Top targets: AD direction & inferred ac4C role ===')
print(tb.to_string())

# Figure: mediation path diagram (text-based) + top target directions
fig,ax=plt.subplots(1,2,figsize=(11,4))
# left: path coefficients
ax[0].axis('off')
txt=(f"Mediation (n={int(med['n'])})\n"
      f"Braak → NAT10 : a={med['path_X_M(a)']:.3f}\n"
      f"NAT10 → ac4C-target : b={med['path_M_Y(b)']:.3f}\n"
      f"Indirect (a×b) : {med['indirect(axb)']:.4f}\n"
      f"Direct (c') : {med['direct_cprime']:.4f}\n"
      f"Total (c) : {med['total_c']:.4f}")
ax[0].text(0.05,0.5,txt,fontsize=11,va='center',family='monospace')
# right: top target bar
cols={'DOWN':'#d95f02','UP':'#2c7fb8'}
ax[1].barh(tb['symbol'], tb['AD_log2FC'], color=[cols[d] for d in tb['AD_direction']])
ax[1].axvline(0,color='k',lw=0.5); ax[1].set_xlabel('AD log2FC'); ax[1].set_title('Top AD-risk ac4C targets in AD brains')
plt.tight_layout(); plt.savefig(os.path.join(FIG,'mediation_directionality.png'),dpi=150,bbox_inches='tight')
plt.savefig(os.path.join(FIG,'mediation_directionality.pdf'),bbox_inches='tight')
print('\nSaved mediation_directionality.png')
print('DONE')
