"""Independent replication of NAT10 downregulation in AD — GSE95587 (2nd cohort).
Free alternative to Synapse protein data: replicates the NAT10-mRNA finding in an
independent AD brain RNA-seq cohort, strengthening the causal claim.
"""
import pandas as pd, numpy as np, os
from scipy import stats
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

HERE=os.path.dirname(os.path.abspath(__file__))
ROOT=os.environ.get('AD_NAT10_DATA_DIR', os.path.abspath(os.path.join(HERE, '..', 'data', 'external')))
OUT=os.path.join(HERE,'..','results'); FIG=os.path.join(HERE,'..','figures')
os.makedirs(OUT, exist_ok=True)

meta=pd.read_csv(os.path.join(ROOT,'processed','GSE95587_metadata.csv'))
counts=pd.read_parquet(os.path.join(ROOT,'processed','GSE95587_raw_counts.parquet'))
counts.index=[str(x).split('_sample')[0] if '_sample' in str(x) else str(x) for x in counts.index]
common=[g for g in meta['GSM'] if g in set(counts.index)]
c=counts.loc[common]; m=meta[meta['GSM'].isin(common)].reset_index(drop=True)
c=c.loc[m['GSM'].values].reset_index(drop=True)
print('aligned', c.shape, m.shape)

# find group column
grp_col=None
for col in m.columns:
    vals=m[col].astype(str).str.lower()
    if vals.str.contains('control|alzheimer|disease|ad|case').any() and m[col].nunique()<=4:
        grp_col=col; break
print('group col:', grp_col, '->', m[grp_col].unique()[:6] if grp_col else 'NONE')

if grp_col and '55226' in [str(x) for x in c.columns]:
    expr=np.log2(c['55226']+1)
    m['nat10']=expr.values
    # normalize group labels
    def lbl(v):
        v=str(v).lower()
        if 'control' in v or 'healthy' in v or v=='0': return 'Control'
        return 'AD'
    m['grp']=m[grp_col].map(lbl)
    ad=m[m.grp=='AD']['nat10']; ct=m[m.grp=='Control']['nat10']
    t,p=stats.ttest_ind(ad,ct,equal_var=False)
    u,pu=stats.mannwhitneyu(ad,ct)
    lfc=float(ad.mean()-ct.mean())
    print(f'GSE95587 NAT10 AD vs Ctrl: AD={ad.mean():.3f} Ctrl={ct.mean():.3f} log2FC={lfc:.3f} t-p={p:.2e} MWU-p={pu:.2e} (nAD={len(ad)},nC={len(ct)})')
    pd.DataFrame([{'cohort':'GSE95587','AD_mean':ad.mean(),'Ctrl_mean':ct.mean(),'log2FC':lfc,
                   't_p':p,'MWU_p':pu,'n_AD':len(ad),'n_Ctrl':len(ct)}]).to_csv(
                   os.path.join(OUT,'nat10_gse95587_replication.csv'), index=False)
    # combine with GSE125583 for meta-effect
    g1=pd.read_csv(os.path.join(OUT,'nat10_gse125583_ad.csv'))
    comb=pd.DataFrame([
        {'cohort':'GSE125583','log2FC':g1['log2FC'].iloc[0],'p':g1['t_test_p'].iloc[0],'n_AD':g1['n_AD'].iloc[0],'n_Ctrl':g1['n_Ctrl'].iloc[0]},
        {'cohort':'GSE95587','log2FC':lfc,'p':p,'n_AD':len(ad),'n_Ctrl':len(ct)}])
    # Stouffer meta p
    from scipy.stats import norm
    z=np.sign(comb['log2FC'])*norm.ppf(1-comb['p']/2)
    # weight by n
    w=np.sqrt(comb['n_AD']+comb['n_Ctrl'])
    z_meta=np.sum(z*w)/np.sqrt(np.sum(w**2))
    p_meta=2*(1-norm.cdf(abs(z_meta)))
    print(f'META (Stouffer, n-weighted): log2FC~{comb.log2FC.mean():.3f}, z={z_meta:.2f}, p_meta={p_meta:.2e}')
    comb['meta_z']=z_meta; comb['meta_p']=p_meta
    comb.to_csv(os.path.join(OUT,'nat10_two_cohort_meta.csv'), index=False)
    # figure
    fig,ax=plt.subplots(1,2,figsize=(7,3.4))
    ax[0].boxplot([ct,ad],labels=['Control','AD']); ax[0].set_ylabel('log2 NAT10'); ax[0].set_title(f'GSE95587 NAT10\nlog2FC={lfc:.2f}, p={pu:.1e}')
    ax[1].bar(comb['cohort'],comb['log2FC'],color=['#2c7fb8','#d95f02'],yerr=[[0,0],[0,0]])
    ax[1].axhline(0,color='k',lw=0.5); ax[1].set_ylabel('log2FC (AD/CTRL)'); ax[1].set_title(f'Two-cohort meta\np={p_meta:.1e}')
    plt.tight_layout(); plt.savefig(os.path.join(FIG,'nat10_replication.png'),dpi=150,bbox_inches='tight')
    plt.savefig(os.path.join(FIG,'nat10_replication.pdf'),bbox_inches='tight')
    print('Saved nat10_replication.png')
print('DONE')
