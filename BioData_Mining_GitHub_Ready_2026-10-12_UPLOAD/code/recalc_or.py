import pandas as pd
import numpy as np
from scipy import stats
import os

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get('AD_NAT10_DATA_DIR', os.path.join(HERE, '..', 'data', 'external'))
results_dir = os.path.join(HERE, '..', 'results')

targets = pd.read_csv(os.path.join(HERE, '..', 'data', 'human_ac4c_targets_curated.csv'))
target_genes = set(targets['gene_symbol'].dropna().astype(str))
print(f'ac4C targets: {len(target_genes)}')

mod = pd.read_csv(os.path.join(results_dir, 'vulnerability_module_enrichment.csv'))
print('\n=== ORIGINAL CSV VALUES ===')
print(mod[['module','module_gene_count','ac4c_target_count','overlap_count','odds_ratio','p_value','fold_enrichment']].to_string())

N_bg = 20162
print('\n=== RECALCULATED OR (N_bg=20162) ===')
recalc = []
for _, row in mod.iterrows():
    a = int(row['overlap_count'])
    b = int(row['ac4c_target_count']) - a
    c = int(row['module_gene_count']) - a
    d = N_bg - a - b - c
    odds, pval = stats.fisher_exact([[a,b],[c,d]], alternative='greater')
    fold = (a / row['ac4c_target_count']) / (row['module_gene_count'] / N_bg)
    recalc.append({
        'module': row['module'],
        'overlap': a,
        'module_size': int(row['module_gene_count']),
        'OR_recalculated': odds,
        'p_recalculated': pval,
        'fold_enrichment': fold,
    })
    print(f"  {row['module']}: {a}/{int(row['ac4c_target_count'])} overlap, OR={odds:.1f}, fold={fold:.1f}, p={pval:.1e}")

# Compare with original
print('\n=== COMPARISON: original vs recalculated ===')
for _, row in mod.iterrows():
    for r in recalc:
        if r['module'] == row['module']:
            diff = abs(row['odds_ratio'] - r['OR_recalculated'])
            print(f"  {row['module']}: orig={row['odds_ratio']:.1f}, recalculated={r['OR_recalculated']:.1f}, diff={diff:.2f}")

# Save recalculated
pd.DataFrame(recalc).to_csv(os.path.join(results_dir, 'module_enrichment_recalculated_v2.csv'), index=False)
print('\nSaved: module_enrichment_recalculated_v2.csv')
