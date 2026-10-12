"""
Reviewer-proof GWAS enrichment permutation with multi-feature matching.
Confounders controlled: transcript abundance (expression), transcript length
(and GC content when available). RNA-modification target sets are known to be
biased by these technical features, so we match the permuted background on them.

Methods reported in manuscript:
  A. simple_random          - 303 genes drawn from all expressed genes
  B. expression_matched     - background restricted to genes within 0.5 log2 of
                              any ac4C target's expression
  C. expr_length_matched    - background restricted to genes matching BOTH
                              expression (0.5 log2) AND length (0.5 log10 bp)
  D. expr_length_gc_matched - same as C plus GC content (if GTEx/annotation GC
                              provided; otherwise falls back to C)
"""
import pandas as pd
import numpy as np
import os

np.random.seed(42)

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get('AD_NAT10_DATA_DIR', os.path.join(HERE, '..', 'data', 'external'))
results_dir = os.path.join(HERE, '..', 'results')

# ---------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------
targets = pd.read_csv(os.path.join(HERE, '..', 'data', 'human_ac4c_targets_curated.csv'))
idmap = pd.read_csv(os.path.join(BASE, 'processed', 'entrez_to_symbol.tsv'), sep='\t')
id_col = [c for c in idmap.columns if 'entrez' in c.lower() or 'geneid' in c.lower()][0]
sym_col = [c for c in idmap.columns if 'symbol' in c.lower() or 'gene' in c.lower()][0]
mapping = dict(zip(idmap[id_col].astype(str), idmap[sym_col].astype(str)))
ac4c_genes_raw = targets['gene_symbol'].dropna().astype(str).tolist()
ac4c_genes = list(set([mapping.get(g, g) for g in ac4c_genes_raw]))

de = pd.read_csv(os.path.join(BASE, 'processed', 'DE_results_with_symbol.csv'))
de['_symbol'] = de['gene_symbol'].astype(str)
expr_map = dict(zip(de['_symbol'], de['baseMean']))

gwas = pd.read_csv(os.path.join(results_dir, 'ad_gwas_enrichment.csv'))
gwas_genes = set(g.strip() for g in str(gwas.iloc[0]['overlap_genes']).split(';') if g.strip())

all_genes = de['_symbol'].tolist()
all_expr = de['baseMean'].values
gene_expr_dict = dict(zip(all_genes, all_expr))

# Gene length (from GENCODE v44, longest protein-coding transcript)
length_file = os.path.join(BASE, 'processed', 'gene_lengths.tsv')
if os.path.exists(length_file):
    lengths = pd.read_csv(length_file, sep='\t')
    length_map = dict(zip(lengths['gene_symbol'], lengths['length']))
    print(f"Gene length file: {len(length_map)} genes")
else:
    length_map = {}
    print("No gene length file - length matching skipped")

# Optional GC content (if provided separately); not available from GTF -> skip
gc_map = {}

print(f"ac4C targets: {len(ac4c_genes)}")
print(f"Genes with expression: {len(all_genes)}")
print(f"AD GWAS genes: {len(gwas_genes)}")

# ---------------------------------------------------------------
# 2. Feature-matched sampling pools
# ---------------------------------------------------------------
ac4c_with_expr = [(g, expr_map[g]) for g in ac4c_genes if g in expr_map]
ac4c_expr = np.array([e for _, e in ac4c_with_expr])
n_targets = len(ac4c_genes)
n_overlap_observed = len(set(ac4c_genes) & gwas_genes)
print(f"ac4C targets with expression: {len(ac4c_with_expr)}")
print(f"Observed overlap: {n_overlap_observed}")

# Expression-matched pool
expr_pool = set()
for e in ac4c_expr:
    for g, ge in gene_expr_dict.items():
        if abs(np.log2(ge + 1) - np.log2(e + 1)) < 0.5:
            expr_pool.add(g)
print(f"Expression-matched pool: {len(expr_pool)}")

# Expression + length matched pool
len_pool = set()
ac4c_len = np.array([length_map[g] for g in ac4c_genes if g in length_map and length_map[g] > 0])
if len(ac4c_len) > 0:
    for g in expr_pool:
        if g not in length_map or length_map[g] <= 0:
            continue
        for tl in ac4c_len:
            if abs(np.log10(length_map[g]) - np.log10(tl)) < 0.5:
                len_pool.add(g)
                break
print(f"Expr+Length matched pool: {len(len_pool)}")

# Expr + length + GC matched pool (GC optional)
gc_pool = set(len_pool)
if gc_map:
    ac4c_gc = np.array([gc_map[g] for g in ac4c_genes if g in gc_map])
    for g in len_pool:
        if g not in gc_map:
            continue
        for tgc in ac4c_gc:
            if abs(gc_map[g] - tgc) < 0.05:
                gc_pool.add(g)
                break
    print(f"Expr+Length+GC matched pool: {len(gc_pool)}")
else:
    print("GC map not available - using expr+length pool for matched method")

# ---------------------------------------------------------------
# 3. Permutations
# ---------------------------------------------------------------
N_ITER = 10000

def permute(pool, label):
    res = []
    pool_list = list(pool)
    if len(pool_list) < n_targets:
        print(f"  WARNING {label}: pool {len(pool_list)} < {n_targets}, falling back to all genes")
        pool_list = all_genes
    for _ in range(N_ITER):
        sample = set(np.random.choice(pool_list, size=n_targets, replace=False))
        res.append(len(sample & gwas_genes))
    arr = np.array(res)
    z = (n_overlap_observed - arr.mean()) / arr.std()
    p = np.mean(arr >= n_overlap_observed)
    return arr.mean(), arr.std(), z, p

print("\n" + "=" * 60)
print("MULTI-FEATURE MATCHED GWAS PERMUTATION")
print("=" * 60)

simple_mean, simple_std, simple_z, simple_p = permute(set(all_genes), "simple")
expr_mean, expr_std, expr_z, expr_p = permute(expr_pool, "expr")
len_mean, len_std, len_z, len_p = permute(len_pool, "expr+len")

print(f"\nSimple random      : exp={simple_mean:.2f} z={simple_z:.1f} p={simple_p:.2e}")
print(f"Expr-matched       : exp={expr_mean:.2f} z={expr_z:.1f} p={expr_p:.2e}")
print(f"Expr+Length-matched: exp={len_mean:.2f} z={len_z:.1f} p={len_p:.2e}")

rows = [
    ('simple_random', n_overlap_observed, simple_mean, simple_std, simple_z, simple_p),
    ('expression_matched', n_overlap_observed, expr_mean, expr_std, expr_z, expr_p),
    ('expr_length_matched', n_overlap_observed, len_mean, len_std, len_z, len_p),
]
if gc_map:
    gc_mean, gc_std, gc_z, gc_p = permute(gc_pool, "expr+len+gc")
    print(f"Expr+Length+GC     : exp={gc_mean:.2f} z={gc_z:.1f} p={gc_p:.2e}")
    rows.append(('expr_length_gc_matched', n_overlap_observed, gc_mean, gc_std, gc_z, gc_p))

result_df = pd.DataFrame(rows, columns=['method', 'observed', 'expected', 'std', 'z_score', 'p_value'])
result_df.to_csv(os.path.join(results_dir, 'gwas_enrichment_expression_matched.csv'), index=False)
print("\nSaved: gwas_enrichment_expression_matched.csv")
