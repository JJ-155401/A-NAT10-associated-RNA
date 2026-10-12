"""
Sensitivity analysis: Re-run GWAS enrichment with alternative ac4C target lists.

Tests robustness of the ac4C-AD GWAS enrichment finding to choice of target gene list.
"""
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
from statsmodels.stats.multitest import multipletests
import os

# ============================================================
# CONFIG
# ============================================================
HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("AD_NAT10_DATA_DIR", HERE.parent / "data" / "external"))
REPO_ROOT = HERE.parent
RESULTS_DIR = REPO_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Load full AD GWAS gene list (89 genes from Bellenguez 2022).
# The derived input is kept under the external-data root so this repository
# does not depend on the original local project layout.
ad_gwas_file = ROOT / "derived" / "ad_gwas_risk_genes_real.csv"
ad_gwas_df = pd.read_csv(ad_gwas_file)
ad_gwas_full = set(ad_gwas_df['gene'].dropna().unique())
print(f"AD GWAS genes (Bellenguez 2022): {len(ad_gwas_full)}")

# Load background gene list
# The background is all genes expressed in the brain dataset
# From the stored result: background_n = 20162
# Let me load the DE results to get the full gene list
de_file = ROOT / "processed" / "DE_results_with_symbol.csv"
de_df = pd.read_csv(de_file)
all_genes = set(de_df['gene_symbol'].dropna().unique())
print(f"Total genes in dataset: {len(all_genes)}")

# ============================================================
# DEFINE ALTERNATIVE ac4C TARGET LISTS
# ============================================================
print("\n" + "=" * 70)
print("DEFINING ALTERNATIVE ac4C TARGET LISTS")
print("=" * 70)

# List 1: Current curated list (303 genes)
current_file = REPO_ROOT / "data" / "human_ac4c_targets_curated.csv"
current_df = pd.read_csv(current_file)
current_targets = set(current_df['gene_symbol'].dropna().unique())
print(f"\nList 1 - Current curated: {len(current_targets)} genes")

# List 2: Arango 2018 acRIP-seq peaks (2,139 genes)
arango_file = ROOT / "downloads" / "ac4c_targets" / "arango2018_ac4c_targets.csv"
arango_df = pd.read_csv(arango_file)
arango_targets = set(arango_df['gene_symbol'].dropna().unique())
print(f"List 2 - Arango 2018 acRIP-seq: {len(arango_targets)} genes")

# List 3: Core overlap (44 genes) - most conservative
core_targets = current_targets & arango_targets
print(f"List 3 - Core overlap: {len(core_targets)} genes")

# List 4: Arango strict (top peaks by score)
peaks_file = ROOT / "downloads" / "ac4c_targets" / "GSE102113_acRIP_peaks.csv"
peaks_df = pd.read_csv(peaks_file)
# Select peaks with score >= 50 (top peaks)
strict_peaks = peaks_df[peaks_df['score'] >= 50]
arango_strict = set(strict_peaks['geneSymbol'].dropna().unique())
print(f"List 4 - Arango strict (score>=50): {len(arango_strict)} genes")

# ============================================================
# AD GWAS GENE LIST (already loaded above)
# ============================================================

# Intersect with background genes
ad_gwas_in_background = ad_gwas_full & all_genes
print(f"AD GWAS genes in background: {len(ad_gwas_in_background)}")

# ============================================================
# RUN ENRICHMENT FOR EACH TARGET LIST
# ============================================================
print("\n" + "=" * 70)
print("GWAS ENRICHMENT ANALYSIS FOR EACH TARGET LIST")
print("=" * 70)

target_lists = {
    'Current curated (303)': current_targets,
    'Arango 2018 acRIP-seq (2139)': arango_targets,
    'Core overlap (44)': core_targets,
    'Arango strict score>=50': arango_strict,
}

results_list = []

for name, targets in target_lists.items():
    print(f"\n--- {name} ---")
    
    # Intersect targets with background
    targets_in_bg = targets & all_genes
    print(f"  Targets in background: {len(targets_in_bg)}")
    
    # Fisher's exact test
    a = len(targets_in_bg & ad_gwas_in_background)  # overlap
    b = len(targets_in_bg - ad_gwas_in_background)  # targets not in GWAS
    c = len(ad_gwas_in_background - targets_in_bg)  # GWAS not in targets
    d = len(all_genes - targets_in_bg - ad_gwas_in_background)  # neither
    
    print(f"  2x2 table: a={a}, b={b}, c={c}, d={d}")
    
    if a > 0:
        odds_ratio, p_value = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
        fold_enrichment = (a / (a + b)) / ((a + c) / (a + b + c + d)) if (a + b) > 0 and (a + c) > 0 else 0
        
        # Get overlap genes
        overlap_genes = targets_in_bg & ad_gwas_in_background
        overlap_str = '; '.join(sorted(overlap_genes))
        
        print(f"  Overlap: {a} genes")
        print(f"  Odds ratio: {odds_ratio:.2f}")
        print(f"  P-value: {p_value:.2e}")
        print(f"  Fold enrichment: {fold_enrichment:.1f}")
        print(f"  Overlap genes: {overlap_str}")
        
        results_list.append({
            'target_list': name,
            'n_targets': len(targets_in_bg),
            'n_ad_gwas': len(ad_gwas_in_background),
            'overlap_n': a,
            'odds_ratio': odds_ratio,
            'p_value': p_value,
            'fold_enrichment': fold_enrichment,
            'overlap_genes': overlap_str,
        })
    else:
        print(f"  No overlap found")
        results_list.append({
            'target_list': name,
            'n_targets': len(targets_in_bg),
            'n_ad_gwas': len(ad_gwas_in_background),
            'overlap_n': 0,
            'odds_ratio': np.nan,
            'p_value': 1.0,
            'fold_enrichment': np.nan,
            'overlap_genes': '',
        })

# ============================================================
# SAVE RESULTS
# ============================================================
results_df = pd.DataFrame(results_list)
results_df.to_csv(RESULTS_DIR / 'gwas_enrichment_sensitivity.csv', index=False)
print(f"\nResults saved to {RESULTS_DIR / 'gwas_enrichment_sensitivity.csv'}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SENSITIVITY ANALYSIS SUMMARY")
print("=" * 70)
print(f"\n{'Target List':<35} {'N':>5} {'Overlap':>8} {'OR':>8} {'P-value':>12} {'Fold':>6}")
print("-" * 80)
for _, row in results_df.iterrows():
    print(f"{row['target_list']:<35} {row['n_targets']:>5} {row['overlap_n']:>8} {row['odds_ratio']:>8.1f} {row['p_value']:>12.2e} {row['fold_enrichment']:>6.1f}")

# Check if results are robust
n_sig = (results_df['p_value'] < 0.05).sum()
print(f"\nSignificant (p<0.05): {n_sig}/{len(results_df)} target lists")
if n_sig == len(results_df):
    print("RESULT: Enrichment is ROBUST across all target list choices")
elif n_sig > 0:
    print("RESULT: Enrichment is PARTIALLY robust (some lists significant)")
else:
    print("RESULT: Enrichment is NOT robust to target list choice")
