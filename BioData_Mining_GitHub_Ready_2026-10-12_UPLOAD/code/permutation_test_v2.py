"""
Permutation test v2: for each iteration, randomly sample 303 genes from the genome
and count overlap with each module.
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
target_set = set(targets['gene_symbol'].dropna().astype(str))
print(f'ac4C targets: {len(target_set)}')

# Module gene sets from CSV
mod = pd.read_csv(os.path.join(results_dir, 'vulnerability_module_enrichment.csv'))
modules = {}
for _, row in mod.iterrows():
    if pd.isna(row['overlap_genes']):
        continue
    overlap_genes = set(g.strip() for g in str(row['overlap_genes']).split(';') if g.strip())
    modules[row['module']] = {
        'size': int(row['module_gene_count']),
        'overlap': overlap_genes,
        'all_genes': overlap_genes,  # We only have overlap genes, not full module
    }

# For the permutation, we need the full genome gene list
# Use the overlap genes + non-overlap as approximation
# Better: load a comprehensive gene list
# Let's use the STRING background
string_file = os.path.join(results_dir, 'string_enrichment_nat10_ac4c_ad.csv')
string_data = pd.read_csv(string_file)

# Extract all unique genes from STRING data as background
all_string_genes = set()
for col in string_data.columns:
    if string_data[col].dtype == object:
        for val in string_data[col].dropna():
            genes = str(val).split(';')
            for g in genes:
                g = g.strip()
                if g and g != 'nan':
                    all_string_genes.add(g)

print(f'Background genes from STRING: {len(all_string_genes)}')

# Use a genome-wide gene list instead
# All protein-coding genes: ~20162
# We'll simulate by sampling from target_set size (303) and counting overlap with each module
# The key insight: we're testing whether 303 randomly chosen genes overlap with modules more than chance

# For a proper permutation:
# H0: ac4C targets are randomly distributed across the genome
# Test: sample 303 random genes, count how many fall in each module

# We don't have the full genome gene list, so we'll use a different approach:
# For each module, test whether its overlap with ac4C targets exceeds what's expected
# by sampling from the universe of all measured genes

# Alternative: use the module gene lists directly
# The Fisher exact test already gives us the correct p-values
# The permutation test should validate these

# Let me use a simpler approach: permute the gene labels
N_bg = 20162
N_iter = 10000
n_targets = len(target_set)

print('\n=== PERMUTATION TEST (sampling from genome) ===')
perm_results = []
for mod_name, info in modules.items():
    mod_size = info['size']
    observed_overlap = len(info['overlap'])
    
    # Randomly sample mod_size genes from a simulated genome of N_bg genes
    # and count how many overlap with a random sample of n_targets genes
    perm_overlaps = []
    for _ in range(N_iter):
        # Sample n_targets random gene indices from genome
        random_targets = set(np.random.choice(N_bg, size=n_targets, replace=False))
        # Sample mod_size random gene indices from genome
        random_mod = set(np.random.choice(N_bg, size=mod_size, replace=False))
        perm_overlaps.append(len(random_targets & random_mod))
    
    perm_overlaps = np.array(perm_overlaps)
    if perm_overlaps.std() > 0:
        z = (observed_overlap - perm_overlaps.mean()) / perm_overlaps.std()
        p = np.mean(perm_overlaps >= observed_overlap)
    else:
        z = np.inf
        p = 0.0
    
    perm_results.append({
        'module': mod_name,
        'observed': observed_overlap,
        'expected_perm': perm_overlaps.mean(),
        'std_perm': perm_overlaps.std(),
        'z_score': z,
        'p_permutation': p,
    })
    print(f'  {mod_name}: observed={observed_overlap}, expected={perm_overlaps.mean():.1f}±{perm_overlaps.std():.1f}, z={z:.1f}, p={p:.2e}')

# AD GWAS
gwas = pd.read_csv(os.path.join(results_dir, 'ad_gwas_enrichment.csv'))
gwas_overlap_count = int(gwas.iloc[0]['overlap_n'])
n_gwas = int(gwas.iloc[0]['ad_gwas_genes_n'])

perm_gwas = []
for _ in range(N_iter):
    random_targets = set(np.random.choice(N_bg, size=n_targets, replace=False))
    random_gwas = set(np.random.choice(N_bg, size=n_gwas, replace=False))
    perm_gwas.append(len(random_targets & random_gwas))

perm_gwas = np.array(perm_gwas)
z_gwas = (gwas_overlap_count - perm_gwas.mean()) / perm_gwas.std()
p_gwas = np.mean(perm_gwas >= gwas_overlap_count)
print(f'\n  AD GWAS: observed={gwas_overlap_count}, expected={perm_gwas.mean():.1f}±{perm_gwas.std():.1f}, z={z_gwas:.1f}, p={p_gwas:.2e}')

perm_df = pd.DataFrame(perm_results)
perm_df.to_csv(os.path.join(results_dir, 'module_permutation_test_v2.csv'), index=False)
print('\nSaved: module_permutation_test_v2.csv')
