"""
Issue 3+4: ac4C target sensitivity analysis + OR recalculation + permutation test
"""
import numpy as np
import pandas as pd
from scipy import stats
import os, warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get('AD_NAT10_DATA_DIR', os.path.join(HERE, '..', 'data', 'external'))
results_dir = os.path.join(HERE, '..', 'results')

# ============================================================
# Load data
# ============================================================
targets = pd.read_csv(os.path.join(HERE, '..', 'data', 'human_ac4c_targets_curated.csv'))
print(f"Loaded {len(targets)} ac4C targets")
print(f"Columns: {targets.columns.tolist()}")

# Assume 'GeneID' or 'entrez' column
gene_col = None
for c in targets.columns:
    if 'gene' in c.lower() or 'entrez' in c.lower() or 'symbol' in c.lower():
        gene_col = c
        break
if gene_col is None:
    gene_col = targets.columns[0]
print(f"Using gene column: {gene_col}")
target_genes = set(targets[gene_col].dropna().astype(str))

# Load gene universe (background)
# Use the string_enrichment file to get the background
string_file = os.path.join(results_dir, 'string_enrichment_nat10_ac4c_ad.csv')
if os.path.exists(string_file):
    string_data = pd.read_csv(string_file)
    print(f"STRING data loaded: {len(string_data)} rows")

# ============================================================
# 1. Define ac4C target tiers
# ============================================================
# Tier 1: Core targets (high-confidence, from primary literature)
# Tier 2: Extended targets (moderate confidence)
# Tier 3: All targets

# For sensitivity, we'll use different subsets
# Since we don't have confidence scores, we'll do random subsampling

# ============================================================
# 2. Module enrichment OR recalculation
# ============================================================
print("\n" + "="*60)
print("MODULE ENRICHMENT OR RECALCULATION")
print("="*60)

# Define modules manually based on known gene sets
modules = {
    'Synaptic_transmission': [
        'GRIN1','GRIN2A','GRIN2B','GRIN2C','GRIN2D','GRIN3A','GRIN3B',
        'GRIA1','GRIA2','GRIA3','GRIA4','GRIA5',
        'GABRA1','GABRA2','GABRA3','GABRA4','GABRA5','GABRA6',
        'GABRB1','GABRB2','GABRB3','GABRG1','GABRG2','GABRG3',
        'CHRNA7','CHRNB2','CHRNB3','CHRNB4','CHRNA4',
        'DRD1','DRD2','DRD3','DRD4','DRD5',
        'HTR1A','HTR1B','HTR2A','HTR2C','HTR3A','HTR4','HTR6','HTR7',
        'SYN1','SYN2','SYN3','SYP','SNAP25','STX1A','STX1B','STXBP1',
        'VAMP1','VAMP2','VAMP3','VAMP5',
        'SLC17A6','SLC17A7','SLC17A8','SLC32A1','SLC6A3','SLC6A4',
        'TH','DBH','TPH2',
    ],
    'Mitochondrial_function': [
        'ATP5F1A','ATP5F1B','ATP5F1C','ATP5F1D','ATP5F1E','ATP5F1G',
        'ATP5MC1','ATP5MC2','ATP5MC3','ATP5PD','ATP5PF',
        'COX5A','COX5B','COX6A','COX6B','COX7A','COX7B','COX8A',
        'NDUFA1','NDUFA2','NDUFA3','NDUFA4','NDUFA5','NDUFA6','NDUFA7','NDUFA8','NDUFA9','NDUFA10','NDUFA11','NDUFA12','NDUFA13',
        'NDUFB1','NDUFB2','NDUFB3','NDUFB4','NDUFB5','NDUFB6','NDUFB7','NDUFB8','NDUFB9','NDUFB10','NDUFB11',
        'NDUFS1','NDUFS2','NDUFS3','NDUFS4','NDUFS5','NDUFS6','NDUFS7','NDUFS8',
        'NDUFV1','NDUFV2','NDUFV3',
        'SDHA','SDHB','SDHC','SDHD',
        'UQCR10','UQCR11','UQCRB','UQCRC1','UQCRC2','UQCRFS1','UQCRH','UQCRQ',
        'CS','IDH1','IDH2','MDH1','MDH2','FH',
    ],
    'Neuronal_resilience': [
        'ATF4','ATF5','ATF6','CREB1','CREB3','CREM',
        'HSPA8','HSPA5','HSP90AA1','HSPB1',
        'NFKB1','NFKB2','RELA','RELB',
        'BCL2','BCL2L1','BAX','CASP3','CASP9',
        'BDNF','NGF','NT3','NT4',
    ],
    'Vesicle_cycling_transport': [
        'SNAP25','STX1A','STX1B','STXBP1',
        'VAMP1','VAMP2','VAMP3','VAMP5',
        'SYT1','SYT2','CPLX1','CPLX2',
        'NSF','SCAMP1','SCAMP3',
    ],
    'Axonal_transport': [
        'KIF1A','KIF1B','KIF1C','KIF21A','KIF21B','KIF2A','KIF3A','KIF3B','KIF5A','KIF5B','KIF5C',
        'DYNC1H1','DYNC1I1','DYNC1I2','DYNE1',
        'MAPT','MAP1A','MAP1B','MAP2',
        'BACE1','BACE2',
        'TUBA1A','TUBA1B','TUBB','TUBB3','TUBB4B',
    ],
}

# Background gene count (total protein-coding genes in human genome)
N_background = 20162

results = []
for mod_name, mod_genes in modules.items():
    mod_set = set(mod_genes)
    overlap = target_genes & mod_set
    a = len(overlap)
    b = len(target_genes - mod_set)
    c = len(mod_set - target_genes)
    d = N_background - a - b - c
    
    # Fisher exact test
    odds, pval = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
    
    # Recalculate manually to verify
    OR_manual = (a * d) / (b * c) if (b * c) > 0 else float('inf')
    fold = (a / len(target_genes)) / (len(mod_set) / N_background)
    
    results.append({
        'module': mod_name,
        'module_size': len(mod_set),
        'overlap': a,
        'OR_fisher': odds,
        'OR_manual': OR_manual,
        'fold_enrichment': fold,
        'p_value': pval,
        'overlap_genes': ';'.join(sorted(overlap)),
    })

df_results = pd.DataFrame(results)
print(df_results[['module','module_size','overlap','OR_fisher','OR_manual','fold_enrichment','p_value']].to_string())

# ============================================================
# 3. Permutation test (10,000 random draws)
# ============================================================
print("\n" + "="*60)
print("PERMUTATION TEST (10,000 iterations)")
print("="*60)

target_list = list(target_genes)
perm_results = []

for mod_name, mod_genes in modules.items():
    mod_set = set(mod_genes)
    observed_overlap = len(target_genes & mod_set)
    
    perm_overlaps = []
    for _ in range(10000):
        random_sample = set(np.random.choice(target_list, size=len(mod_set), replace=False))
        perm_overlaps.append(len(random_sample & mod_set))
    
    perm_overlaps = np.array(perm_overlaps)
    z_score = (observed_overlap - perm_overlaps.mean()) / perm_overlaps.std()
    p_perm = np.mean(perm_overlaps >= observed_overlap)
    
    perm_results.append({
        'module': mod_name,
        'observed_overlap': observed_overlap,
        'expected_perm': perm_overlaps.mean(),
        'std_perm': perm_overlaps.std(),
        'z_score': z_score,
        'p_permutation': p_perm,
        'OR_empirical': observed_overlap / max(perm_overlaps.mean(), 0.1),
    })
    print(f"  {mod_name}: observed={observed_overlap}, expected={perm_overlaps.mean():.1f} +/- {perm_overlaps.std():.1f}, z={z_score:.2f}, p={p_perm:.4f}")

perm_df = pd.DataFrame(perm_results)

# ============================================================
# 4. Sensitivity: random subsampling of ac4C targets
# ============================================================
print("\n" + "="*60)
print("SENSITIVITY: ac4C target subsampling (1000 iterations)")
print("="*60)

synapse_mod = set(modules['Synaptic_transmission'])
mito_mod = set(modules['Mitochondrial_function'])

sensitivity_syn = []
sensitivity_mito = []
target_sizes = [50, 100, 150, 200, 250, 303]

for target_n in target_sizes:
    if target_n > len(target_list):
        continue
    syn_ors = []
    mito_ors = []
    for _ in range(1000):
        sample = set(np.random.choice(target_list, size=target_n, replace=False))
        
        # Synaptic
        a = len(sample & synapse_mod)
        b = target_n - a
        c = len(synapse_mod - sample)
        d = N_background - a - b - c
        if b > 0 and c > 0:
            o, _ = stats.fisher_exact([[a,b],[c,d]], alternative='greater')
            syn_ors.append(o)
        
        # Mito
        a2 = len(sample & mito_mod)
        b2 = target_n - a2
        c2 = len(mito_mod - sample)
        d2 = N_background - a2 - b2 - c2
        if b2 > 0 and c2 > 0:
            o2, _ = stats.fisher_exact([[a2,b2],[c2,d2]], alternative='greater')
            mito_ors.append(o2)
    
    sensitivity_syn.append({
        'target_n': target_n,
        'median_OR': np.median(syn_ors),
        'ci_low': np.percentile(syn_ors, 2.5),
        'ci_high': np.percentile(syn_ors, 97.5),
    })
    sensitivity_mito.append({
        'target_n': target_n,
        'median_OR': np.median(mito_ors),
        'ci_low': np.percentile(mito_ors, 2.5),
        'ci_high': np.percentile(mito_ors, 97.5),
    })
    print(f"  n={target_n}: Synapse OR={np.median(syn_ors):.1f} [{np.percentile(syn_ors,2.5):.1f}-{np.percentile(syn_ors,97.5):.1f}], Mito OR={np.median(mito_ors):.1f} [{np.percentile(mito_ors,2.5):.1f}-{np.percentile(mito_ors,97.5):.1f}]")

# Save
df_results.to_csv(os.path.join(results_dir, 'module_enrichment_recalculated.csv'), index=False)
perm_df.to_csv(os.path.join(results_dir, 'module_enrichment_permutation.csv'), index=False)
pd.DataFrame(sensitivity_syn).to_csv(os.path.join(results_dir, 'sensitivity_synapse_or.csv'), index=False)
pd.DataFrame(sensitivity_mito).to_csv(os.path.join(results_dir, 'sensitivity_mito_or.csv'), index=False)

print("\n=== SAVED ===")
print("module_enrichment_recalculated.csv")
print("module_enrichment_permutation.csv")
print("sensitivity_synapse_or.csv")
print("sensitivity_mito_or.csv")
