"""
Cross-validation of NAT10/ac4C single-cell findings:
Zhou et al. 2020 (GSE147528) vs Mathys et al. 2019 (GSE138852)

Runs IDENTICAL pipeline on both datasets and reports concordance.
"""
import warnings
warnings.filterwarnings('ignore')

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path
from statsmodels.stats.multitest import multipletests
import os

# ============================================================
# CONFIG
# ============================================================
REPO_ROOT = Path(__file__).resolve().parents[1]
EXTERNAL_DATA = Path(os.environ.get("AD_NAT10_DATA_DIR", REPO_ROOT / "data" / "external"))
TARGETS_PATH = REPO_ROOT / "data" / "human_ac4c_targets_curated.csv"
MATHYS_COUNTS = EXTERNAL_DATA / "scrna" / "GSE138852_counts.csv.gz"
MATHYS_COVARS = EXTERNAL_DATA / "scrna" / "GSE138852_covariates.csv.gz"
ZHOU_RESULTS = REPO_ROOT / "results" / "single_cell_analysis.csv"
OUT_DIR = REPO_ROOT / "results"
FIG_DIR = REPO_ROOT / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 1
sc.settings.figdir = str(FIG_DIR)

print("=" * 70)
print("SINGLE-CELL CROSS-VALIDATION: Zhou 2020 vs Mathys 2019")
print("=" * 70)

# ============================================================
# STEP 1: Load Mathys 2019 data
# ============================================================
print("\n[1] Loading Mathys 2019 data...")

counts = pd.read_csv(MATHYS_COUNTS, index_col=0, compression='gzip')
covars = pd.read_csv(MATHYS_COVARS, index_col=0, compression='gzip')
print(f"  Counts: {counts.shape} (genes x cells)")
print(f"  Covariates: {covars.shape}")

# Align cells
common_cells = counts.columns.intersection(covars.index)
counts = counts[common_cells]
covars = covars.loc[common_cells]
print(f"  Aligned cells: {len(common_cells)}")

# Create AnnData (transpose to cells x genes)
adata = ad.AnnData(X=counts.T.values.astype(np.float32),
                   obs=covars.copy(),
                   var=pd.DataFrame(index=counts.index))
adata.var_names_make_unique()
print(f"  AnnData: {adata.shape}")

# ============================================================
# STEP 2: QC and normalization (identical to Zhou)
# ============================================================
print("\n[2] QC and normalization...")

# Rename column for consistency
adata.obs['cell_type'] = adata.obs['oupSample.cellType']
adata.obs['group'] = adata.obs['oupSample.batchCond'].map({'AD': 'AD', 'ct': 'Control'})
print(f"  Groups: {adata.obs['group'].value_counts().to_dict()}")
print(f"  Cell types: {adata.obs['cell_type'].value_counts().to_dict()}")

# QC metrics
sc.pp.calculate_qc_metrics(adata, inplace=True, log1p=False)
print(f"  Before filtering: {adata.shape}")
print(f"  Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")

# Apply same filters as Zhou
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
mask = adata.obs['total_counts'] >= 500
adata = adata[mask].copy()
print(f"  After filtering: {adata.shape}")

# Normalize (identical to Zhou: target_sum=1e4 + log1p)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['counts'] = adata.X.copy()
print("  Normalized (target_sum=1e4, log1p)")

# ============================================================
# STEP 3: NAT10 expression by cell type
# ============================================================
print("\n[3] NAT10 expression analysis...")

if 'NAT10' in adata.var_names:
    nat10_expr = adata[:, 'NAT10'].X
    if hasattr(nat10_expr, 'toarray'):
        nat10_expr = nat10_expr.toarray().flatten()
    elif hasattr(nat10_expr, 'A1'):
        nat10_expr = nat10_expr.A1
    else:
        nat10_expr = np.array(nat10_expr).flatten()
    adata.obs['NAT10_expr'] = nat10_expr
    print(f"  NAT10 found. Expression range: {nat10_expr.min():.2f} - {nat10_expr.max():.2f}")
    print(f"  Fraction expressing: {(nat10_expr > 0).mean():.2%}")
else:
    print("  WARNING: NAT10 not found!")
    adata.obs['NAT10_expr'] = 0.0

# Mean NAT10 by cell type
nat10_by_ct = adata.obs.groupby('cell_type')['NAT10_expr'].agg(['mean', 'median', 'std', 'count'])
print("\n  NAT10 by cell type:")
print(nat10_by_ct.to_string())

# ============================================================
# STEP 4: ac4C target module score (identical to Zhou)
# ============================================================
print("\n[4] ac4C target module score...")

targets_df = pd.read_csv(TARGETS_PATH)
ac4c_targets = targets_df['gene_symbol'].unique().tolist()
targets_in_data = [g for g in ac4c_targets if g in adata.var_names]
print(f"  ac4C targets: {len(ac4c_targets)} total, {len(targets_in_data)} found ({len(targets_in_data)/len(ac4c_targets)*100:.1f}%)")

# Compute module score (identical to Zhou: ctrl_size=200)
sc.tl.score_genes(adata, targets_in_data, score_name='ac4c_module', ctrl_size=200)
print(f"  ac4C module range: {adata.obs['ac4c_module'].min():.3f} - {adata.obs['ac4c_module'].max():.3f}")

module_by_ct = adata.obs.groupby('cell_type')['ac4c_module'].agg(['mean', 'median', 'std', 'count'])
print("\n  ac4C module by cell type:")
print(module_by_ct.to_string())

# ============================================================
# STEP 5: AD vs Control comparison per cell type (identical to Zhou)
# ============================================================
print("\n[5] AD vs Control comparison per cell type...")

results_list = []
cell_types = [ct for ct in adata.obs['cell_type'].unique() if ct not in ('Unknown', 'unID')]

for ct in sorted(cell_types):
    mask_ct = adata.obs['cell_type'] == ct
    n_ctrl = (mask_ct & (adata.obs['group'] == 'Control')).sum()
    n_ad = (mask_ct & (adata.obs['group'] == 'AD')).sum()

    if n_ctrl < 10 or n_ad < 10:
        print(f"\n  {ct}: Too few cells (Ctrl={n_ctrl}, AD={n_ad}), skipping")
        continue

    ctrl_data = adata.obs.loc[mask_ct & (adata.obs['group'] == 'Control')]
    ad_data = adata.obs.loc[mask_ct & (adata.obs['group'] == 'AD')]

    # NAT10 expression comparison (Mann-Whitney U, identical to Zhou)
    nat10_ctrl = ctrl_data['NAT10_expr']
    nat10_ad = ad_data['NAT10_expr']
    stat_nat10, pval_nat10 = stats.mannwhitneyu(nat10_ad, nat10_ctrl, alternative='two-sided')
    fold_nat10 = nat10_ad.mean() / (nat10_ctrl.mean() + 1e-10)

    # ac4C module comparison
    mod_ctrl = ctrl_data['ac4c_module']
    mod_ad = ad_data['ac4c_module']
    stat_mod, pval_mod = stats.mannwhitneyu(mod_ad, mod_ctrl, alternative='two-sided')
    fold_mod = mod_ad.mean() / (mod_ctrl.mean() + 1e-10)

    # Vulnerability score
    adata.obs.loc[mask_ct, 'nat10_z'] = (adata.obs.loc[mask_ct, 'NAT10_expr'] - adata.obs['NAT10_expr'].mean()) / (adata.obs['NAT10_expr'].std() + 1e-10)
    adata.obs.loc[mask_ct, 'ac4c_z'] = (adata.obs.loc[mask_ct, 'ac4c_module'] - adata.obs['ac4c_module'].mean()) / (adata.obs['ac4c_module'].std() + 1e-10)
    adata.obs.loc[mask_ct, 'vulnerability'] = (adata.obs.loc[mask_ct, 'nat10_z'] + adata.obs.loc[mask_ct, 'ac4c_z']) / 2

    ctrl_vuln = ctrl_data['vulnerability'] if 'vulnerability' in ctrl_data.columns else adata.obs.loc[mask_ct & (adata.obs['group'] == 'Control'), 'vulnerability']
    ad_vuln = ad_data['vulnerability'] if 'vulnerability' in ad_data.columns else adata.obs.loc[mask_ct & (adata.obs['group'] == 'AD'), 'vulnerability']
    stat_vuln, pval_vuln = stats.mannwhitneyu(ad_vuln, ctrl_vuln, alternative='two-sided')

    row = {
        'cell_type': ct,
        'n_ctrl': n_ctrl,
        'n_ad': n_ad,
        'NAT10_mean_ctrl': nat10_ctrl.mean(),
        'NAT10_mean_ad': nat10_ad.mean(),
        'NAT10_fold_change': fold_nat10,
        'NAT10_pvalue': pval_nat10,
        'ac4C_module_mean_ctrl': mod_ctrl.mean(),
        'ac4C_module_mean_ad': mod_ad.mean(),
        'ac4C_module_fold_change': fold_mod,
        'ac4C_module_pvalue': pval_mod,
        'vuln_mean_ctrl': ctrl_vuln.mean(),
        'vuln_mean_ad': ad_vuln.mean(),
        'vuln_pvalue': pval_vuln,
    }
    results_list.append(row)

    print(f"\n  {ct} (Ctrl={n_ctrl}, AD={n_ad}):")
    print(f"    NAT10: Ctrl={nat10_ctrl.mean():.4f}, AD={nat10_ad.mean():.4f}, FC={fold_nat10:.3f}, p={pval_nat10:.2e}")
    print(f"    ac4C module: Ctrl={mod_ctrl.mean():.4f}, AD={mod_ad.mean():.4f}, FC={fold_mod:.3f}, p={pval_mod:.2e}")
    print(f"    Vulnerability: Ctrl={ctrl_vuln.mean():.4f}, AD={ad_vuln.mean():.4f}, p={pval_vuln:.2e}")

results_df = pd.DataFrame(results_list)

# FDR correction
for pcol in ['NAT10_pvalue', 'ac4C_module_pvalue', 'vuln_pvalue']:
    if pcol in results_df.columns:
        valid = results_df[pcol].notna()
        if valid.sum() > 0:
            _, fdr, _, _ = multipletests(results_df.loc[valid, pcol].values, method='fdr_bh')
            results_df.loc[valid, pcol.replace('_pvalue', '_fdr')] = fdr

# Save Mathys results
results_df.to_csv(OUT_DIR / 'single_cell_mathys2019.csv', index=False)
print(f"\n  Mathys results saved to {OUT_DIR / 'single_cell_mathys2019.csv'}")

# ============================================================
# STEP 6: Load Zhou 2020 results and compare
# ============================================================
print("\n" + "=" * 70)
print("[6] CONCORDANCE ANALYSIS: Zhou 2020 vs Mathys 2019")
print("=" * 70)

zhou_df = pd.read_csv(ZHOU_RESULTS)
print(f"\n  Zhou 2020 results: {len(zhou_df)} cell types")
print(f"  Mathys 2019 results: {len(results_df)} cell types")

# Map cell type names for comparison
zhou_ct_map = {
    'Excitatory neuron': 'neuron',
    'Inhibitory neuron': 'neuron',
    'Astrocyte': 'astro',
    'Microglia': 'mg',
    'Oligodendrocyte': 'oligo',
    'OPC': 'OPC',
    'Endothelial': 'endo',
}

# Aggregate Zhou excitatory + inhibitory into 'neuron' for comparison
zhou_neuron = zhou_df[zhou_df['cell_type'].isin(['Excitatory neuron', 'Inhibitory neuron'])]
if len(zhou_neuron) > 0:
    zhou_neuron_agg = pd.DataFrame([{
        'cell_type': 'neuron',
        'n_ctrl': zhou_neuron['n_ctrl'].sum(),
        'n_ad': zhou_neuron['n_ad'].sum(),
        'NAT10_mean_ctrl': np.average(zhou_neuron['NAT10_mean_ctrl'], weights=zhou_neuron['n_ctrl']),
        'NAT10_mean_ad': np.average(zhou_neuron['NAT10_mean_ad'], weights=zhou_neuron['n_ad']),
        'ac4C_module_mean_ctrl': np.average(zhou_neuron['ac4C_module_mean_ctrl'], weights=zhou_neuron['n_ctrl']),
        'ac4C_module_mean_ad': np.average(zhou_neuron['ac4C_module_mean_ad'], weights=zhou_neuron['n_ad']),
    }])
    zhou_comp = pd.concat([zhou_df[~zhou_df['cell_type'].isin(['Excitatory neuron', 'Inhibitory neuron'])], zhou_neuron_agg], ignore_index=True)
else:
    zhou_comp = zhou_df.copy()

# Standardize cell type names in Zhou
zhou_comp['cell_type_std'] = zhou_comp['cell_type'].map(zhou_ct_map).fillna(zhou_comp['cell_type'])

# Compare NAT10 expression direction (AD vs Ctrl)
print("\n  --- NAT10 Expression Direction (AD vs Ctrl) ---")
print(f"  {'Cell Type':<15} {'Zhou FC':>10} {'Mathys FC':>10} {'Same Dir?':>10}")
print(f"  {'-'*45}")

concordant_nat10 = 0
total_ct = 0
for ct in results_df['cell_type'].values:
    zhou_row = zhou_comp[zhou_comp['cell_type_std'] == ct]
    mathys_row = results_df[results_df['cell_type'] == ct]
    if len(zhou_row) > 0 and len(mathys_row) > 0:
        zhou_fc = zhou_row['NAT10_fold_change'].values[0]
        mathys_fc = mathys_row['NAT10_fold_change'].values[0]
        same_dir = (zhou_fc > 1 and mathys_fc > 1) or (zhou_fc < 1 and mathys_fc < 1) or (abs(zhou_fc - 1) < 0.05 and abs(mathys_fc - 1) < 0.05)
        if same_dir:
            concordant_nat10 += 1
        total_ct += 1
        print(f"  {ct:<15} {zhou_fc:>10.3f} {mathys_fc:>10.3f} {'YES' if same_dir else 'NO':>10}")

print(f"\n  NAT10 direction concordance: {concordant_nat10}/{total_ct} cell types")

# Compare ac4C module direction
print("\n  --- ac4C Module Direction (AD vs Ctrl) ---")
print(f"  {'Cell Type':<15} {'Zhou FC':>10} {'Mathys FC':>10} {'Same Dir?':>10}")
print(f"  {'-'*45}")

concordant_ac4c = 0
total_ac4c = 0
for ct in results_df['cell_type'].values:
    zhou_row = zhou_comp[zhou_comp['cell_type_std'] == ct]
    mathys_row = results_df[results_df['cell_type'] == ct]
    if len(zhou_row) > 0 and len(mathys_row) > 0:
        zhou_fc = zhou_row['ac4C_module_fold_change'].values[0]
        mathys_fc = mathys_row['ac4C_module_fold_change'].values[0]
        same_dir = (zhou_fc > 1 and mathys_fc > 1) or (zhou_fc < 1 and mathys_fc < 1) or (abs(zhou_fc - 1) < 0.05 and abs(mathys_fc - 1) < 0.05)
        if same_dir:
            concordant_ac4c += 1
        total_ac4c += 1
        print(f"  {ct:<15} {zhou_fc:>10.3f} {mathys_fc:>10.3f} {'YES' if same_dir else 'NO':>10}")

print(f"\n  ac4C module direction concordance: {concordant_ac4c}/{total_ac4c} cell types")

# Correlation of effect sizes
print("\n  --- Effect Size Correlation ---")
common_cts = set(results_df['cell_type'].values) & set(zhou_comp['cell_type_std'].values)
if len(common_cts) >= 3:
    zhou_nat10_fcs = []
    mathys_nat10_fcs = []
    zhou_ac4c_fcs = []
    mathys_ac4c_fcs = []
    for ct in common_cts:
        zr = zhou_comp[zhou_comp['cell_type_std'] == ct].iloc[0]
        mr = results_df[results_df['cell_type'] == ct].iloc[0]
        zhou_nat10_fcs.append(zr['NAT10_fold_change'])
        mathys_nat10_fcs.append(mr['NAT10_fold_change'])
        zhou_ac4c_fcs.append(zr['ac4C_module_fold_change'])
        mathys_ac4c_fcs.append(mr['ac4C_module_fold_change'])

    rho_nat10, p_nat10 = stats.spearmanr(zhou_nat10_fcs, mathys_nat10_fcs)
    rho_ac4c, p_ac4c = stats.spearmanr(zhou_ac4c_fcs, mathys_ac4c_fcs)
    print(f"  NAT10 FC correlation: Spearman rho={rho_nat10:.3f}, p={p_nat10:.3e} (n={len(common_cts)})")
    print(f"  ac4C FC correlation:  Spearman rho={rho_ac4c:.3f}, p={p_ac4c:.3e} (n={len(common_cts)})")
else:
    print(f"  Only {len(common_cts)} common cell types, cannot compute correlation")

# ============================================================
# STEP 7: Summary concordance metrics
# ============================================================
print("\n" + "=" * 70)
print("[7] CONCORDANCE SUMMARY")
print("=" * 70)

# Compute overall concordance score
overall_score = (concordant_nat10 + concordant_ac4c) / (total_ct + total_ac4c) * 100 if (total_ct + total_ac4c) > 0 else 0
print(f"\n  Overall direction concordance: {overall_score:.1f}%")
print(f"  NAT10 concordance: {concordant_nat10}/{total_ct}")
print(f"  ac4C concordance: {concordant_ac4c}/{total_ac4c}")

# Effect of cell-type composition on bulk signal
print("\n  --- Cell-Type Composition Effect ---")
# From Mathys: neuron fraction in AD vs Control
neuron_mask = adata.obs['cell_type'] == 'neuron'
neuron_frac_ad = (neuron_mask & (adata.obs['group'] == 'AD')).sum() / (adata.obs['group'] == 'AD').sum()
neuron_frac_ctrl = (neuron_mask & (adata.obs['group'] == 'Control')).sum() / (adata.obs['group'] == 'Control').sum()
print(f"  Neuron fraction: AD={neuron_frac_ad:.3f}, Control={neuron_frac_ctrl:.3f}")
print(f"  Neuron depletion in AD: delta={neuron_frac_ad - neuron_frac_ctrl:.3f}")

# ============================================================
# STEP 8: Generate cross-validation figure
# ============================================================
print("\n[8] Generating cross-validation figure...")

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# Panel A: NAT10 FC comparison (Zhou vs Mathys)
ax = axes[0, 0]
if len(common_cts) >= 3:
    ax.scatter(zhou_nat10_fcs, mathys_nat10_fcs, s=100, c='#2c7fb8', edgecolors='black', linewidth=0.5, zorder=5)
    for i, ct in enumerate(list(common_cts)):
        ax.annotate(ct, (zhou_nat10_fcs[i], mathys_nat10_fcs[i]), fontsize=8, ha='center', va='bottom')
    ax.axhline(1, color='gray', linestyle='--', linewidth=0.5)
    ax.axvline(1, color='gray', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Zhou 2020 NAT10 FC (AD/Ctrl)')
    ax.set_ylabel('Mathys 2019 NAT10 FC (AD/Ctrl)')
    ax.set_title(f'NAT10 FC Concordance\nSpearman rho={rho_nat10:.2f}, p={p_nat10:.1e}')
    lim = max(abs(min(zhou_nat10_fcs + mathys_nat10_fcs) - 1), abs(max(zhou_nat10_fcs + mathys_nat10_fcs) - 1)) * 1.2
    ax.set_xlim(1 - lim, 1 + lim)
    ax.set_ylim(1 - lim, 1 + lim)

# Panel B: ac4C FC comparison
ax = axes[0, 1]
if len(common_cts) >= 3:
    ax.scatter(zhou_ac4c_fcs, mathys_ac4c_fcs, s=100, c='#d95f02', edgecolors='black', linewidth=0.5, zorder=5)
    for i, ct in enumerate(list(common_cts)):
        ax.annotate(ct, (zhou_ac4c_fcs[i], mathys_ac4c_fcs[i]), fontsize=8, ha='center', va='bottom')
    ax.axhline(1, color='gray', linestyle='--', linewidth=0.5)
    ax.axvline(1, color='gray', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Zhou 2020 ac4C Module FC')
    ax.set_ylabel('Mathys 2019 ac4C Module FC')
    ax.set_title(f'ac4C Module FC Concordance\nSpearman rho={rho_ac4c:.2f}, p={p_ac4c:.1e}')

# Panel C: Concordance summary bar
ax = axes[0, 2]
categories = ['NAT10\nDirection', 'ac4C\nDirection', 'Overall']
scores = [concordant_nat10/total_ct*100 if total_ct > 0 else 0,
          concordant_ac4c/total_ac4c*100 if total_ac4c > 0 else 0,
          overall_score]
colors = ['#2c7fb8', '#d95f02', '#666666']
bars = ax.bar(categories, scores, color=colors, edgecolor='black', linewidth=0.5)
ax.set_ylabel('Concordance (%)')
ax.set_title('Cross-Dataset Concordance')
ax.set_ylim(0, 110)
for bar, score in zip(bars, scores):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2, f'{score:.0f}%', ha='center', fontsize=11, fontweight='bold')
ax.axhline(50, color='red', linestyle='--', linewidth=0.5, alpha=0.5, label='Chance level')

# Panel D: NAT10 expression by cell type (Mathys)
ax = axes[1, 0]
ct_order = sorted([ct for ct in adata.obs['cell_type'].unique() if ct not in ('Unknown', 'unID', 'doublet')])
if 'neuron' not in ct_order and 'neuron' in adata.obs['cell_type'].values:
    ct_order.append('neuron')
ct_order = sorted(ct_order)
nat10_data = [adata.obs.loc[adata.obs['cell_type'] == ct, 'NAT10_expr'].values for ct in ct_order]
vp = ax.violinplot(nat10_data, showmeans=True, showmedians=True)
ax.set_xticks(range(1, len(ct_order) + 1))
ax.set_xticklabels(ct_order, rotation=45, ha='right', fontsize=9)
ax.set_ylabel('NAT10 expression')
ax.set_title('Mathys 2019: NAT10 by Cell Type')

# Panel E: ac4C module by cell type (Mathys)
ax = axes[1, 1]
mod_data = [adata.obs.loc[adata.obs['cell_type'] == ct, 'ac4c_module'].values for ct in ct_order]
vp = ax.violinplot(mod_data, showmeans=True, showmedians=True)
ax.set_xticks(range(1, len(ct_order) + 1))
ax.set_xticklabels(ct_order, rotation=45, ha='right', fontsize=9)
ax.set_ylabel('ac4C module score')
ax.set_title('Mathys 2019: ac4C Module by Cell Type')

# Panel F: AD vs Ctrl per cell type (Mathys)
ax = axes[1, 2]
if len(results_df) > 0:
    plot_data = results_df[['cell_type', 'NAT10_fold_change', 'ac4C_module_fold_change']].set_index('cell_type')
    plot_data.plot(kind='bar', ax=ax, color=['#E91E63', '#9C27B0'], edgecolor='black', linewidth=0.5)
    ax.axhline(1, color='gray', linestyle='--', linewidth=0.5)
    ax.set_ylabel('Fold Change (AD/Ctrl)')
    ax.set_title('Mathys 2019: AD vs Ctrl FC')
    ax.tick_params(axis='x', rotation=45)
    ax.legend(fontsize=8)

plt.tight_layout()
fig.savefig(FIG_DIR / 'cross_validation_sc.pdf', dpi=300, bbox_inches='tight')
fig.savefig(FIG_DIR / 'cross_validation_sc.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"  Figure saved to {FIG_DIR / 'cross_validation_sc.pdf'}")

# ============================================================
# STEP 9: Save concordance summary
# ============================================================
summary = {
    'analysis': 'single_cell_cross_validation',
    'zhou_dataset': 'GSE147528 (Zhou 2020)',
    'mathys_dataset': 'GSE138852 (Mathys 2019)',
    'zhou_n_cells': len(zhou_df),
    'mathys_n_cells': adata.n_obs,
    'n_common_cell_types': len(common_cts),
    'nat10_direction_concordance': f'{concordant_nat10}/{total_ct}',
    'ac4c_direction_concordance': f'{concordant_ac4c}/{total_ac4c}',
    'overall_concordance_pct': round(overall_score, 1),
    'nat10_fc_spearman_rho': round(rho_nat10, 3) if len(common_cts) >= 3 else 'N/A',
    'ac4c_fc_spearman_rho': round(rho_ac4c, 3) if len(common_cts) >= 3 else 'N/A',
}
pd.DataFrame([summary]).to_csv(OUT_DIR / 'cross_validation_summary.csv', index=False)
print(f"\n  Summary saved to {OUT_DIR / 'cross_validation_summary.csv'}")

print("\n" + "=" * 70)
print("CROSS-VALIDATION COMPLETE")
print("=" * 70)
