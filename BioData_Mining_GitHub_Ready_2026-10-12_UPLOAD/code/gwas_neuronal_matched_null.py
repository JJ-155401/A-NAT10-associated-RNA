"""
Supplementary Table S9: Neuronal-expression matched null model for GWAS enrichment.
Reviewer concern: ac4C targets are enriched in neuronal genes; AD GWAS enrichment
may simply reflect this neuronal bias. We test: do random 88-gene sets drawn from
neuronally-expressed genes (high in neurons) also show high ac4C overlap?
If AD GWAS overlap remains far above this neuronal null, the enrichment is specific
to AD genetic architecture, not just neuronal gene bias.
"""
import pandas as pd
import numpy as np
from scipy import stats
import os

np.random.seed(42)
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get('AD_NAT10_DATA_DIR', os.path.join(HERE, '..', 'data', 'external'))
results_dir = os.path.join(HERE, '..', 'results')

# Load ac4C targets
targets = pd.read_csv(os.path.join(HERE, '..', 'data', 'human_ac4c_targets_curated.csv'))
idmap = pd.read_csv(os.path.join(BASE, 'processed', 'entrez_to_symbol.tsv'), sep='\t')
id_col = [c for c in idmap.columns if 'entrez' in c.lower() or 'geneid' in c.lower()][0]
sym_col = [c for c in idmap.columns if 'symbol' in c.lower() or 'gene' in c.lower()][0]
mapping = dict(zip(idmap[id_col].astype(str), idmap[sym_col].astype(str)))
ac4c_set = set(mapping.get(g, g) for g in targets['gene_symbol'].dropna().astype(str))

# Load DE results to get the gene universe
de = pd.read_csv(os.path.join(BASE, 'processed', 'DE_results_with_symbol.csv'))
de['baseMean'] = de['baseMean'].astype(float)
# Use baseMean as proxy for brain expression level
bg_genes = de[['gene_symbol', 'baseMean']].dropna()
bg_genes['gene_symbol'] = bg_genes['gene_symbol'].astype(str)

# Define "neuronally expressed" as top 50% by baseMean in brain tissue
median_expr = bg_genes['baseMean'].median()
neuronal_genes = set(bg_genes[bg_genes['baseMean'] >= median_expr]['gene_symbol'])
print(f"Neuronal-expression universe: {len(neuronal_genes)} genes (top 50% by brain expression)")

# Also define "highly neuronal" as top 30%
top30_threshold = bg_genes['baseMean'].quantile(0.70)
highly_neuronal = set(bg_genes[bg_genes['baseMean'] >= top30_threshold]['gene_symbol'])
print(f"Highly neuronal universe (top 30%): {len(highly_neuronal)} genes")

# AD GWAS genes
gwas = pd.read_csv(os.path.join(results_dir, 'ad_gwas_enrichment.csv')).iloc[0]
gwas_n = int(gwas['ad_gwas_genes_n'])  # 88
overlap_genes = set(g.strip() for g in str(gwas['overlap_genes']).split(';') if g.strip())

# Fishers exact
def fisher_or(k, a, b, N):
    table = [[k, a - k], [b - k, N - a - b + k]]
    return stats.fisher_exact(table)[0]

# For each universe, compute observed + random
results = []
for label, uni in [('neuronal_top50', neuronal_genes), ('neuronal_top30', highly_neuronal)]:
    uni_list = list(uni)
    N_uni = len(uni)
    ac4c_in_uni = ac4c_set & uni
    a_obs = len(ac4c_in_uni)
    
    # Observed: overlap of AD GWAS with ac4C in this universe
    gwas_in_uni = overlap_genes & uni
    k_obs = len(ac4c_in_uni & gwas_in_uni)
    obs_or = fisher_or(k_obs, a_obs, len(gwas_in_uni), N_uni)
    obs_p = stats.fisher_exact([[k_obs, a_obs - k_obs], [len(gwas_in_uni) - k_obs, N_uni - a_obs - len(gwas_in_uni) + k_obs]])[1]
    
    print(f"\n{label}: universe={N_uni}, ac4C_in_uni={a_obs}, gwas_in_uni={len(gwas_in_uni)}, overlap={k_obs}, OR={obs_or:.2f}, p={obs_p:.2e}")
    
    # Random: draw 88-gene sets from this neuronal universe
    N_ITER = 10000
    rand_ors = []
    for _ in range(N_ITER):
        samp = set(np.random.choice(uni_list, gwas_n, replace=False))
        k = len(ac4c_in_uni & samp)
        rand_ors.append(fisher_or(k, a_obs, gwas_n, N_uni))
    rand_ors = np.array(rand_ors)
    p_rand = np.mean(rand_ors >= obs_or)
    z = (obs_or - rand_ors.mean()) / rand_ors.std() if rand_ors.std() > 0 else np.inf
    print(f"  Random neuronal sets: mean OR={rand_ors.mean():.2f}, max={rand_ors.max():.2f}, z={z:.1f}, p={p_rand:.2e}")
    
    results.append({
        'universe': label,
        'universe_size': N_uni,
        'ac4c_targets_in_universe': a_obs,
        'observed_overlap': k_obs,
        'observed_OR': obs_or,
        'observed_p': obs_p,
        'random_mean_OR': rand_ors.mean(),
        'random_max_OR': rand_ors.max(),
        'z_vs_random': z,
        'empirical_p': p_rand,
    })

out = pd.DataFrame(results)
out.to_csv(os.path.join(results_dir, 'gwas_neuronal_matched_null.csv'), index=False)
print("\nSaved: gwas_neuronal_matched_null.csv")
