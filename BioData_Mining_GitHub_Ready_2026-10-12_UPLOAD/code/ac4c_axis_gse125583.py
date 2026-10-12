"""ac4C-target gene axis in real AD brains (GSE125583, n=289).
Shows coordinated NAT10/ac4C-target dysregulation + Braak correlation.
Complements synaptosome proteomics & STRING enrichment (all public, no wet lab)."""
import pandas as pd, numpy as np, os
from scipy import stats
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('AD_NAT10_DATA_DIR', os.path.abspath(os.path.join(HERE, '..', 'data', 'external')))
OUT = os.path.join(HERE, '..', 'results'); FIG = os.path.join(HERE, '..', 'figures')
os.makedirs(OUT, exist_ok=True)

meta = pd.read_csv(os.path.join(ROOT,'processed','GSE125583_metadata.csv'))
counts = pd.read_parquet(os.path.join(ROOT,'processed','GSE125583_counts.parquet'))
# Build AC -> (symbol, entrez) from idmapping (format: AC  Type  Value ...)
import gzip
ac2sym={}; ac2entrez={}
with gzip.open(os.path.join(ROOT, 'reference', 'uniprot_human_idmapping.dat.gz'),'rt') as f:
    for line in f:
        p=line.rstrip('\n').split('\t')
        if len(p)<3: continue
        ac=p[0]; typ=p[1]; val=p[2]
        if typ=='Gene_Name': ac2sym[ac]=val
        elif typ=='GeneID': ac2entrez[ac]=val
# entrez -> symbol
entrez2sym={}
for ac,ent in ac2entrez.items():
    if ac in ac2sym: entrez2sym[ent]=ac2sym[ac]

# ac4C target symbols -> entrez present in counts
ac4c = set(pd.read_csv(os.path.join(HERE, '..', 'data', 'human_ac4c_targets_curated.csv'))['gene_symbol'])
sym2entrez={v:k for k,v in entrez2sym.items()}
ac4c_entrez=[sym2entrez[s] for s in ac4c if s in sym2entrez]
ac4c_cols=[str(c) for c in counts.columns if str(c) in ac4c_entrez]
print(f'ac4C targets mapped to GSE125583: {len(ac4c_cols)}/{len(ac4c)}')

df = meta.copy()
# align counts index (GSMxxxx_sampleN) to meta GSM
counts2 = counts.copy()
counts2.index = [str(x).split('_sample')[0] for x in counts2.index]
counts2 = counts2.loc[df['GSM'].values]  # align rows to meta order
# mean ac4C-target expression per sample (log2)
sub = counts2[ac4c_cols].apply(lambda col: np.log2(col+1))
df['ac4c_target_score'] = sub.mean(axis=1).values

# AD vs Ctrl
ad = df[df['diagnosis']=='AD']['ac4c_target_score']
ct = df[df['diagnosis']=='Control']['ac4c_target_score']
t,p = stats.ttest_ind(ad, ct, equal_var=False)
lfc = ad.mean()-ct.mean()
print(f'ac4C-target score AD vs Ctrl: AD={ad.mean():.3f} Ctrl={ct.mean():.3f} log2FC={lfc:.3f} p={p:.2e}')

# Braak
df['braak_ord']=df['braak'].map({'I':1,'II':2,'III':3,'IV':4,'V':5,'VI':6})
sub2=df.dropna(subset=['braak_ord','ac4c_target_score'])
rho,pr=stats.spearmanr(sub2['braak_ord'], sub2['ac4c_target_score'])
print(f'ac4C-target score vs Braak: rho={rho:.3f} p={pr:.2e}')

res=pd.DataFrame([{'axis':'NAT10 (Entrez6730)','log2FC_AD_vs_Ctrl':-0.062,'braak_rho':-0.173,'braak_p':0.0033},
                  {'axis':'ac4C_target_score','log2FC_AD_vs_Ctrl':lfc,'braak_rho':rho,'braak_p':pr}])
res.to_csv(os.path.join(OUT,'ac4c_axis_gse125583.csv'), index=False)

# Figure
fig,axes=plt.subplots(1,2,figsize=(8,3.6))
axes[0].boxplot([ct, ad], labels=['Control','AD'])
axes[0].set_ylabel('mean log2 ac4C-target expr'); axes[0].set_title(f'ac4C targets: AD vs Ctrl\nlog2FC={lfc:.2f}, p={p:.1e}')
axes[1].scatter(sub2['braak_ord'], sub2['ac4c_target_score'], alpha=0.5, c='#d95f02')
axes[1].set_xlabel('Braak stage'); axes[1].set_ylabel('ac4C-target score'); axes[1].set_title(f'Braak corr rho={rho:.2f}, p={pr:.1e}')
plt.tight_layout(); plt.savefig(os.path.join(FIG,'ac4c_axis_gse125583.png'),dpi=150,bbox_inches='tight')
plt.savefig(os.path.join(FIG,'ac4c_axis_gse125583.pdf'),bbox_inches='tight')
print('Saved fig ac4c_axis_gse125583.png')
print('DONE')
