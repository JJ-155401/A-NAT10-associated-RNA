"""NAT10 mRNA in real AD brains — GSE125583 RNA-seq (289 samples, AD/Control + Braak).
Fills the mRNA pillar (real AD brain) to complement protein/synaptosome evidence.
Also quantifies ac4C-target genes' coordinated downregulation in AD brains.
"""
import pandas as pd, numpy as np, os
from scipy import stats
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('AD_NAT10_DATA_DIR', os.path.abspath(os.path.join(HERE, '..', 'data', 'external')))
OUT = os.path.join(HERE, '..', 'results'); FIG = os.path.join(HERE, '..', 'figures')
os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)

meta = pd.read_csv(os.path.join(ROOT,'processed','GSE125583_metadata.csv'))
counts = pd.read_parquet(os.path.join(ROOT,'processed','GSE125583_counts.parquet'))
print('counts', counts.shape, '(samples x genes)')

# genes are columns (Entrez IDs); NAT10 Entrez = 55226 (NOT 6730=SRP68)
cols = set(str(x) for x in counts.columns)
nat10_col = '55226' if '55226' in cols else None
print('NAT10 column (Entrez 55226):', nat10_col)
counts = counts[str(nat10_col)] if nat10_col else None
if counts is not None:
    df = meta.copy()
    df['NAT10_expr'] = counts.values
    # AD vs Control
    ad = df[df['diagnosis']=='AD']['NAT10_expr']
    ct = df[df['diagnosis']=='Control']['NAT10_expr']
    # use log2(raw+1)
    ad_l = np.log2(ad+1); ct_l = np.log2(ct+1)
    t,p = stats.ttest_ind(ad_l, ct_l, equal_var=False)
    u,pu = stats.mannwhitneyu(ad_l, ct_l)
    lfc = ad_l.mean()-ct_l.mean()
    print(f'NAT10 AD vs Ctrl: AD={ad_l.mean():.3f} Ctrl={ct_l.mean():.3f} log2FC={lfc:.3f} t-p={p:.2e} MWU-p={pu:.2e}')
    # Braak correlation (ordered)
    df['braak_ord'] = df['braak'].map({'I':1,'II':2,'III':3,'IV':4,'V':5,'VI':6})
    sub = df.dropna(subset=['braak_ord','NAT10_expr'])
    rho,pr = stats.spearmanr(sub['braak_ord'], np.log2(sub['NAT10_expr']+1))
    print(f'NAT10 vs Braak: rho={rho:.3f} p={pr:.2e} (n={len(sub)})')
    # Save
    res = {'gene':'NAT10','AD_mean_log2':ad_l.mean(),'Ctrl_mean_log2':ct_l.mean(),
           'log2FC':lfc,'t_test_p':p,'MWU_p':pu,'braak_rho':rho,'braak_p':pr,'n_AD':len(ad),'n_Ctrl':len(ct)}
    pd.DataFrame([res]).to_csv(os.path.join(OUT,'nat10_gse125583_ad.csv'), index=False)

    # Figure: NAT10 AD vs Ctrl + Braak
    fig,axes=plt.subplots(1,2,figsize=(8,3.6))
    axes[0].boxplot([ct_l, ad_l], labels=['Control','AD'])
    axes[0].set_ylabel('log2(NAT10 expr+1)'); axes[0].set_title(f'NAT10: AD vs Ctrl\nlog2FC={lfc:.2f}, p={pu:.1e}')
    axes[1].scatter(sub['braak_ord'], np.log2(sub['NAT10_expr']+1), alpha=0.5, c='#2c7fb8')
    axes[1].set_xlabel('Braak stage'); axes[1].set_ylabel('log2 NAT10'); axes[1].set_title(f'Braak corr rho={rho:.2f}')
    plt.tight_layout(); plt.savefig(os.path.join(FIG,'nat10_gse125583.png'),dpi=150,bbox_inches='tight')
    plt.savefig(os.path.join(FIG,'nat10_gse125583.pdf'),bbox_inches='tight')
    print('Saved fig nat10_gse125583.png')
print('DONE')
