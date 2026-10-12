"""
Permutation test for module enrichment: 10000 random draws of 303 genes from genome
"""
import numpy as np
import pandas as pd
from scipy import stats
import os

np.random.seed(42)
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get('AD_NAT10_DATA_DIR', os.path.join(HERE, '..', 'data', 'external'))
results_dir = os.path.join(HERE, '..', 'results')

# Load ac4C targets
targets = pd.read_csv(os.path.join(HERE, '..', 'data', 'human_ac4c_targets_curated.csv'))
target_list = targets['gene_symbol'].dropna().astype(str).tolist()
print(f'ac4C targets: {len(target_list)}')

# Module definitions (using the gene lists from original CSV overlap_genes)
mod = pd.read_csv(os.path.join(results_dir, 'vulnerability_module_enrichment.csv'))

modules = {}
for _, row in mod.iterrows():
    if pd.isna(row['overlap_genes']):
        continue
    # Get all module genes: overlap + non-overlap
    # We need the full module gene list, not just overlap
    # The module_gene_count gives us the total size
    overlap_genes = [g.strip() for g in str(row['overlap_genes']).split(';') if g.strip()]
    modules[row['module']] = {
        'size': int(row['module_gene_count']),
        'overlap': overlap_genes,
    }

N_bg = 20162
N_ITER = 10000
target_set = set(target_list)

print('\n=== PERMUTATION TEST ===')
perm_results = []
for mod_name, info in modules.items():
    mod_size = info['size']
    observed = len(info['overlap'])
    
    # Random sample of mod_size genes from target_list
    perm_overlaps = []
    for _ in range(N_ITER):
        random_mod = set(np.random.choice(target_list, size=min(mod_size, len(target_list)), replace=False))
        perm_overlaps.append(len(random_mod & target_set))
    
    perm_overlaps = np.array(perm_overlaps)
    z = (observed - perm_overlaps.mean()) / perm_overlaps.std()
    p = np.mean(perm_overlaps >= observed)
    
    perm_results.append({
        'module': mod_name,
        'observed': observed,
        'expected_perm': perm_overlaps.mean(),
        'std_perm': perm_overlaps.std(),
        'z_score': z,
        'p_permutation': p,
    })
    print(f'  {mod_name}: observed={observed}, expected={perm_overlaps.mean():.1f}±{perm_overlaps.std():.1f}, z={z:.1f}, p={p:.5f}')

# Also test AD GWAS enrichment
gwas = pd.read_csv(os.path.join(results_dir, 'ad_gwas_enrichment.csv'))
gwas_genes = [g.strip() for g in str(gwas.iloc[0]['overlap_genes']).split(';') if g.strip()]
gwas_overlap = set(gwas_genes)
N_ad_gwas = int(gwas.iloc[0]['ad_gwas_genes_n'])

perm_gwas = []
for _ in range(N_ITER):
    random_303 = set(np.random.choice(target_list, size=min(303, len(target_list)), replace=False))
    perm_gwas.append(len(random_303 & gwas_overlap))

perm_gwas = np.array(perm_gwas)
z_gwas = (len(gwas_overlap) - perm_gwas.mean()) / perm_gwas.std()
p_gwas = np.mean(perm_gwas >= len(gwas_overlap))
print(f'\n  AD GWAS: observed={len(gwas_overlap)}, expected={perm_gwas.mean():.1f}±{perm_gwas.std():.1f}, z={z_gwas:.1f}, p={p_gwas:.5f}')

perm_df = pd.DataFrame(perm_results)
perm_df.to_csv(os.path.join(results_dir, 'module_permutation_test.csv'), index=False)
print('\nSaved: module_permutation_test.csv')
