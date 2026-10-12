"""
Cleaner NVS test: bulk-based vulnerability (AD differential expression)
- Are ac4C target genes more vulnerable (larger |log2FC|) than non-targets in GSE125583?
- This is independent of gene-set overlap (no circularity)
"""
import pandas as pd
import numpy as np
from scipy import stats
import os

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get('AD_NAT10_DATA_DIR', os.path.join(HERE, '..', 'data', 'external'))
results_dir = os.path.join(HERE, '..', 'results')

# Load ac4C targets
targets = pd.read_csv(os.path.join(HERE, '..', 'data', 'human_ac4c_targets_curated.csv'))
idmap = pd.read_csv(os.path.join(BASE, 'processed', 'entrez_to_symbol.tsv'), sep='\t')
id_col = [c for c in idmap.columns if 'entrez' in c.lower() or 'geneid' in c.lower()][0]
sym_col = [c for c in idmap.columns if 'symbol' in c.lower() or 'gene' in c.lower()][0]
mapping = dict(zip(idmap[id_col].astype(str), idmap[sym_col].astype(str)))
ac4c_genes_raw = targets['gene_symbol'].dropna().astype(str).tolist()
ac4c_genes = set([mapping.get(g, g) for g in ac4c_genes_raw])

# Load GSE125583 DE results
de = pd.read_csv('results/ac4c_axis_gse125583.csv')
print("DE columns:", de.columns.tolist())
print(de.head())

# Check what columns we have
if 'gene_symbol' in de.columns:
    de_genes = set(de['gene_symbol'].dropna().astype(str))
elif 'gene' in de.columns:
    de_genes = set(de['gene'].dropna().astype(str))
else:
    # Assume first column is gene
    de_genes = set(de.iloc[:, 0].dropna().astype(str))

# Find log2FC column
fc_col = None
for c in de.columns:
    if 'log2' in c.lower() or 'fc' in c.lower() or 'lfc' in c.lower():
        fc_col = c
        break
if fc_col is None:
    # Try to find a column with numeric values that could be FC
    for c in de.columns:
        if de[c].dtype in ['float64', 'float32'] and 'p' not in c.lower():
            fc_col = c
            break

print(f"\nUsing FC column: {fc_col}")

# Test: ac4C targets vs non-targets
de_targets = de[de['gene_symbol'].isin(ac4c_genes)] if 'gene_symbol' in de.columns else de[de.iloc[:,0].isin(ac4c_genes)]
de_nontargets = de[~de['gene_symbol'].isin(ac4c_genes)] if 'gene_symbol' in de.columns else de[~de.iloc[:,0].isin(ac4c_genes)]

ac4c_fc = de_targets[fc_col].dropna().abs()
nontarget_fc = de_nontargets[fc_col].dropna().abs()

# Mann-Whitney U test
u_stat, u_p = stats.mannwhitneyu(ac4c_fc, nontarget_fc, alternative='greater')
# t-test
t_stat, t_p = stats.ttest_ind(ac4c_fc, nontarget_fc, equal_var=False)

print(f"\n=== Bulk vulnerability test (GSE125583) ===")
print(f"ac4C targets: n={len(ac4c_fc)}, mean|log2FC|={ac4c_fc.mean():.4f}")
print(f"Non-targets:  n={len(nontarget_fc)}, mean|log2FC|={nontarget_fc.mean():.4f}")
print(f"Mann-Whitney U: p={u_p:.2e}")
print(f"Welch t-test: p={t_p:.2e}")

# Also test signed (downregulation = vulnerability)
ac4c_fc_signed = de_targets[fc_col].dropna()
nontarget_fc_signed = de_nontargets[fc_col].dropna()
# More negative = more vulnerable
u_stat2, u_p2 = stats.mannwhitneyu(-ac4c_fc_signed, -nontarget_fc_signed, alternative='greater')
print(f"\nSigned (down=more vulnerable): Mann-Whitney p={u_p2:.2e}")
print(f"  ac4C targets mean log2FC: {ac4c_fc_signed.mean():.4f}")
print(f"  Non-targets mean log2FC: {nontarget_fc_signed.mean():.4f}")

# Save
result = {
    'test': 'bulk_DE_vulnerability',
    'ac4c_n': len(ac4c_fc),
    'ac4c_mean_abs_fc': ac4c_fc.mean(),
    'nontarget_n': len(nontarget_fc),
    'nontarget_mean_abs_fc': nontarget_fc.mean(),
    'mwu_p': u_p,
    'ttest_p': t_p,
    'signed_mwu_p': u_p2,
    'ac4c_mean_signed_fc': ac4c_fc_signed.mean(),
    'nontarget_mean_signed_fc': nontarget_fc_signed.mean(),
}
pd.DataFrame([result]).to_csv(os.path.join(results_dir, 'nvs_bulk_de_test.csv'), index=False)
print("\nSaved: nvs_bulk_de_test.csv")
