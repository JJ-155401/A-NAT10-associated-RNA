"""
GWAS leave-one-out sensitivity (Reviewer must-fix #3):
Recompute Fisher exact OR/p for the 303 ac4C targets vs the 88 AD GWAS genes
(Bellenguez 2022) in a background of 20,162 genes, after removing each
overlapping gene (and key-locus groups) one at a time.
gwas_n is held at the true 88 (only the removed overlap gene leaves the GWAS set).
"""
import pandas as pd
import numpy as np
from scipy import stats
import os

HERE = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(HERE, '..', 'results')
gwas = pd.read_csv(os.path.join(results_dir, 'ad_gwas_enrichment.csv')).iloc[0]

ac4c_n = int(gwas['ac4c_targets_n'])     # 303
gwas_n = int(gwas['ad_gwas_genes_n'])    # 88 (true full set)
bg = int(gwas['background_n'])           # 20162
overlap_genes = [g.strip() for g in str(gwas['overlap_genes']).split(';') if g.strip()]
print(f"ac4C targets={ac4c_n}, AD GWAS genes={gwas_n}, background={bg}, overlap={len(overlap_genes)}")

def fisher_or(k, a, b, N):
    table = [[k, a - k], [b - k, N - a - b + k]]
    orr, p = stats.fisher_exact(table)
    return orr, p

base_or, base_p = fisher_or(len(overlap_genes), ac4c_n, gwas_n, bg)
print(f"Baseline: overlap={len(overlap_genes)} OR={base_or:.2f} p={base_p:.2e}")

rows = [('baseline (none removed)', len(overlap_genes), base_or, base_p)]

for g in overlap_genes:
    k = len(overlap_genes) - 1
    b = gwas_n - 1
    orr, p = fisher_or(k, ac4c_n, b, bg)
    rows.append((f'remove {g}', k, orr, p))

groups = {
    'remove APOE': ['APOE'],
    'remove APOE+APP': ['APOE', 'APP'],
    'remove top4 (APOE,APP,MAPT,PSEN1)': ['APOE', 'APP', 'MAPT', 'PSEN1'],
    'remove APOE+APP+MAPT+PSEN1+BIN1': ['APOE', 'APP', 'MAPT', 'PSEN1', 'BIN1'],
}
for label, grp in groups.items():
    grp = [g for g in grp if g in overlap_genes]
    if not grp:
        continue
    k = len(overlap_genes) - len(grp)
    b = gwas_n - len(grp)
    orr, p = fisher_or(k, ac4c_n, b, bg)
    rows.append((label, k, orr, p))

out = pd.DataFrame(rows, columns=['scenario', 'overlap_n', 'odds_ratio', 'p_value'])
out.to_csv(os.path.join(results_dir, 'gwas_leave_one_out.csv'), index=False)
print("\nLeave-one-out sensitivity:")
print(out.to_string(index=False))
print("\nSaved: gwas_leave_one_out.csv")
