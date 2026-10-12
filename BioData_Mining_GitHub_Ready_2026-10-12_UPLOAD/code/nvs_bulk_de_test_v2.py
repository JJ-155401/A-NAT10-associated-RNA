"""
Cleaner NVS test: bulk-based vulnerability (AD differential expression from GSE125583)
- Are ac4C target genes more vulnerable (larger |log2FC|) than non-targets?
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

# Load DE results (per-gene)
de = pd.read_csv(os.path.join(BASE, 'processed', 'DE_results_with_symbol.csv'))
print("DE columns:", de.columns.tolist()[:10])
print("DE shape:", de.shape)

# Find gene symbol and log2FC columns
gene_col = None
for c in de.columns:
    if 'symbol' in c.lower() or 'gene' == c.lower():
        gene_col = c
        break
if gene_col is None:
    gene_col = de.columns[0]

fc_col = None
for c in de.columns:
    if 'log2' in c.lower():
        fc_col = c
        break
if fc_col is None:
    # find numeric column
    for c in de.columns:
        if de[c].dtype in ['float64', 'float32'] and 'p' not in c.lower():
            fc_col = c
            break

print(f"Gene col: {gene_col}, FC col: {fc_col}")

# Map genes
de['_in_ac4c'] = de[gene_col].astype(str).isin(ac4c_genes)
de_targets = de[de['_in_ac4c']]
de_nontargets = de[~de['_in_ac4c']]

ac4c_fc = de_targets[fc_col].dropna().abs()
nontarget_fc = de_nontargets[fc_col].dropna().abs()

print(f"\n=== Bulk vulnerability test (GSE125583) ===")
print(f"ac4C targets: n={len(ac4c_fc)}, mean|log2FC|={ac4c_fc.mean():.4f}")
print(f"Non-targets:  n={len(nontarget_fc)}, mean|log2FC|={nontarget_fc.mean():.4f}")

# Mann-Whitney U test
u_stat, u_p = stats.mannwhitneyu(ac4c_fc, nontarget_fc, alternative='greater')
t_stat, t_p = stats.ttest_ind(ac4c_fc, nontarget_fc, equal_var=False)

print(f"Mann-Whitney U (targets > non-targets): p={u_p:.2e}")
print(f"Welch t-test: p={t_p:.2e}")

# Signed test (down = vulnerable)
ac4c_signed = de_targets[fc_col].dropna()
nontarget_signed = de_nontargets[fc_col].dropna()
u_stat2, u_p2 = stats.mannwhitneyu(-ac4c_signed, -nontarget_signed, alternative='greater')
print(f"\nSigned (AD-down = vulnerable): Mann-Whitney p={u_p2:.2e}")
print(f"  ac4C targets mean log2FC: {ac4c_signed.mean():.4f}")
print(f"  Non-targets mean log2FC: {nontarget_signed.mean():.4f}")

# Also: fraction of significantly downregulated genes
sig_threshold = 0.05
if 'p' in de.columns[-1].lower() or 'padj' in de.columns[-1].lower() or 'FDR' in de.columns[-1]:
    p_col = de.columns[-1]
    ac4c_sig_down = (de_targets[fc_col] < -0.1) & (de_targets[p_col] < sig_threshold)
    nontarget_sig_down = (de_nontargets[fc_col] < -0.1) & (de_nontargets[p_col] < sig_threshold)
    print(f"\nFraction downregulated (|log2FC|>0.1, p<{sig_threshold}):")
    print(f"  ac4C targets: {ac4c_sig_down.mean():.3f} ({ac4c_sig_down.sum()}/{len(de_targets)})")
    print(f"  Non-targets:  {nontarget_sig_down.mean():.3f} ({nontarget_sig_down.sum()}/{len(de_nontargets)})")

# Save
result = {
    'test': 'bulk_DE_vulnerability_GSE125583',
    'ac4c_n': len(ac4c_fc),
    'ac4c_mean_abs_fc': ac4c_fc.mean(),
    'nontarget_n': len(nontarget_fc),
    'nontarget_mean_abs_fc': nontarget_fc.mean(),
    'mwu_p': u_p,
    'ttest_p': t_p,
    'signed_mwu_p': u_p2,
    'ac4c_mean_signed_fc': ac4c_signed.mean(),
    'nontarget_mean_signed_fc': nontarget_signed.mean(),
}
pd.DataFrame([result]).to_csv(os.path.join(results_dir, 'nvs_bulk_de_test.csv'), index=False)
print("\nSaved: nvs_bulk_de_test.csv")
