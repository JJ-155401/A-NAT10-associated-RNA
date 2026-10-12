"""RNA-protein concordance analysis: GSE125583 (mRNA) vs PXD055784 (synaptosome protein).
Validates that NAT10 and ac4C target genes show correlated mRNA-protein changes in AD.
"""
import pandas as pd, numpy as np, os, gzip
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('AD_NAT10_DATA_DIR', os.path.abspath(os.path.join(HERE, '..', 'data', 'external')))
REF  = os.path.join(ROOT, 'reference')
OUT  = os.path.join(HERE, '..', 'results')
os.makedirs(OUT, exist_ok=True)

# ========== 1. Load data ==========

# ac4C targets
ac4c = set(pd.read_csv(os.path.join(HERE, '..', 'data', 'human_ac4c_targets_curated.csv'))['gene_symbol'])
print(f'ac4C targets: {len(ac4c)}')

# UniProt -> gene symbol mapping
print('Loading UniProt mapping...')
ac2gene = {}
with gzip.open(os.path.join(REF, 'uniprot_human_idmapping.dat.gz'), 'rt', encoding='utf-8') as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 3 and p[1] == 'Gene_Name':
            ac2gene[p[0]] = p[2]
print(f'  {len(ac2gene)} UniProt->gene mappings')

# Synaptosome protein DE (PXD055784, 5 AD vs 5 HC)
prot_de = pd.read_csv(os.path.join(OUT, 'synaptosome_protein_DE.csv'))
print(f'Synaptosome proteins: {len(prot_de)}')
# map UniProt -> gene
prot_de['gene'] = prot_de['Protein'].map(ac2gene)
prot_mapped = prot_de.dropna(subset=['gene']).copy()
print(f'  mapped to gene symbols: {len(prot_mapped)}')

# GSE125583 RNA-seq - compute per-gene mRNA DE from raw counts (AD vs Ctrl)
meta = pd.read_csv(os.path.join(ROOT, 'processed', 'GSE125583_metadata.csv'))
counts = pd.read_parquet(os.path.join(ROOT, 'processed', 'GSE125583_counts.parquet'))
counts.index = [str(x).split('_sample')[0] for x in counts.index]
counts = counts.loc[meta['GSM'].values].reset_index(drop=True)

# Compute mRNA DE per gene (AD vs Ctrl, log2(mean_AD+1) - log2(mean_CTRL+1))
ad_mask = meta['diagnosis'] == 'AD'
ctrl_mask = meta['diagnosis'] == 'Control'
mrna_de_list = []
for col in counts.columns:
    ad_vals = np.log2(counts.loc[ad_mask, col].values + 1)
    ctrl_vals = np.log2(counts.loc[ctrl_mask, col].values + 1)
    ad_mean = ad_vals.mean()
    ctrl_mean = ctrl_vals.mean()
    lfc = ad_mean - ctrl_mean
    t, p = stats.ttest_ind(ad_vals, ctrl_vals, equal_var=False)
    mrna_de_list.append({'entrez_id': col, 'mrna_log2FC': lfc, 'mrna_pvalue': p})

mrna_de = pd.DataFrame(mrna_de_list)

# Map entrez to symbol
entrez_sym = {}
for _, r in pd.read_csv(os.path.join(ROOT, 'processed', 'entrez_to_symbol.tsv'), sep='\t').iterrows():
    entrez_sym[str(r['entrez_id'])] = str(r['symbol'])
mrna_de['gene'] = mrna_de['entrez_id'].astype(str).map(entrez_sym)
mrna_de = mrna_de.dropna(subset=['gene'])
print(f'GSE125583 mRNA DE: {len(mrna_de)} genes mapped')

# NAT10 specific
nat10_mrna = mrna_de[mrna_de['gene'] == 'NAT10']
print(f'\nNAT10 mRNA DE: log2FC={nat10_mrna["mrna_log2FC"].values[0]:.4f}, p={nat10_mrna["mrna_pvalue"].values[0]:.2e}')

nat10_prot = prot_mapped[prot_mapped['gene'] == 'NAT10']
if len(nat10_prot):
    print(f'NAT10 protein DE: log2FC={nat10_prot["log2FC"].values[0]:.4f}, q={nat10_prot["q"].values[0]:.2e}')
else:
    print('NAT10 protein: not detected in synaptosome dataset (expected for low-abundance enzyme)')

# ========== 2. Merge mRNA and protein data ==========

# Keep genes present in both datasets
mrna_sub = mrna_de[['gene', 'mrna_log2FC', 'mrna_pvalue']]
prot_sub = prot_mapped[['gene', 'log2FC', 'p', 'q']].rename(
    columns={'log2FC': 'prot_log2FC', 'p': 'prot_p', 'q': 'prot_q'})
merged = mrna_sub.merge(prot_sub, on='gene', how='inner')
print(f'\nGenes with both mRNA and protein data: {len(merged)}')

# ========== 3. NAT10 mRNA-protein correlation ==========

# NAT10 expression per individual (meta and counts already loaded above)
meta['NAT10_mrna'] = np.log2(counts['55226'].values + 1)
meta['braak_ord'] = meta['braak'].map({'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6})

# AD vs Ctrl for NAT10 mRNA
ad_nat10 = meta[meta.diagnosis == 'AD']['NAT10_mrna']
ctrl_nat10 = meta[meta.diagnosis == 'Control']['NAT10_mrna']
nat10_t, nat10_p = stats.ttest_ind(ad_nat10, ctrl_nat10, equal_var=False)
nat10_fc = ad_nat10.mean() - ctrl_nat10.mean()

# Braak correlation
sb = meta.dropna(subset=['braak_ord', 'NAT10_mrna'])
nat10_rho, nat10_rho_p = stats.spearmanr(sb['braak_ord'], sb['NAT10_mrna'])

print(f'\n=== NAT10 mRNA in GSE125583 (n=289) ===')
print(f'  AD vs Ctrl: log2FC={nat10_fc:.4f}, t-test p={nat10_p:.2e}')
print(f'  Braak: rho={nat10_rho:.4f}, p={nat10_rho_p:.2e}')

# ========== 4. Overall mRNA-protein concordance ==========

merged['mrna_dir'] = np.sign(merged['mrna_log2FC'])
merged['prot_dir'] = np.sign(merged['prot_log2FC'])
merged['concordant'] = (merged['mrna_dir'] == merged['prot_dir']).astype(int)
# Exclude zero-zero (no change in either)
both_zero = (merged['mrna_log2FC'] == 0) & (merged['prot_log2FC'] == 0)
merged_eval = merged[~both_zero].copy()
n_concordant = merged_eval['concordant'].sum()
n_total = len(merged_eval)
concordance_pct = n_concordant / n_total * 100

print(f'\n=== Overall mRNA-protein concordance (AD direction) ===')
print(f'  Evaluated genes: {n_total}')
print(f'  Concordant: {n_concordant}/{n_total} ({concordance_pct:.1f}%)')

# Test against chance (50%)
binom_p = stats.binomtest(n_concordant, n_total, 0.5, alternative='greater').pvalue
print(f'  vs 50% chance: binomial p={binom_p:.2e}')

# Pearson and Spearman correlation of log2FC
r_pearson, r_p_pearson = stats.pearsonr(merged_eval['mrna_log2FC'], merged_eval['prot_log2FC'])
r_spearman, r_p_spearman = stats.spearmanr(merged_eval['mrna_log2FC'], merged_eval['prot_log2FC'])
print(f'  Pearson r={r_pearson:.4f}, p={r_p_pearson:.2e}')
print(f'  Spearman rho={r_spearman:.4f}, p={r_p_spearman:.2e}')

# ========== 5. ac4C targets vs non-targets ==========

merged_eval['is_ac4c_target'] = merged_eval['gene'].isin(ac4c)

ac4c_sub = merged_eval[merged_eval['is_ac4c_target']]
non_sub = merged_eval[~merged_eval['is_ac4c_target']]

ac4c_n = len(ac4c_sub)
ac4c_cons = ac4c_sub['concordant'].sum()
ac4c_pct = ac4c_cons / ac4c_n * 100 if ac4c_n > 0 else np.nan

non_n = len(non_sub)
non_cons = non_sub['concordant'].sum()
non_pct = non_cons / non_n * 100 if non_n > 0 else np.nan

print(f'\n=== ac4C targets vs non-targets concordance ===')
print(f'  ac4C targets: {ac4c_cons}/{ac4c_n} ({ac4c_pct:.1f}%)')
print(f'  non-targets:  {non_cons}/{non_n} ({non_pct:.1f}%)')

# Fisher's exact test: ac4C targets more concordant?
a = int(ac4c_cons)  # ac4C concordant
b = ac4c_n - a       # ac4C discordant
c = int(non_cons)    # non concordant
d = non_n - c        # non discordant
odds_ratio, fisher_p = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
print(f'  Fisher exact (greater): OR={odds_ratio:.3f}, p={fisher_p:.3e}')

# Also test if ac4C targets have stronger mRNA-protein correlation
if ac4c_n >= 5:
    r_ac, p_ac = stats.spearmanr(ac4c_sub['mrna_log2FC'], ac4c_sub['prot_log2FC'])
else:
    r_ac, p_ac = np.nan, np.nan
r_non, p_non = stats.spearmanr(non_sub['mrna_log2FC'], non_sub['prot_log2FC'])
print(f'  ac4C targets Spearman rho={r_ac:.4f} (p={p_ac:.2e})')
print(f'  non-targets Spearman rho={r_non:.4f} (p={p_non:.2e})')

# Mean absolute log2FC (effect size)
ac4c_mean_abs = ac4c_sub[['mrna_log2FC', 'prot_log2FC']].abs().mean().mean()
non_mean_abs = non_sub[['mrna_log2FC', 'prot_log2FC']].abs().mean().mean()
print(f'  ac4C targets mean |log2FC|={ac4c_mean_abs:.4f}')
print(f'  non-targets mean |log2FC|={non_mean_abs:.4f}')

# ========== 5b. Concordance among significantly changed genes ==========
# Focus on genes significant in at least one dataset (p<0.05)
sig_mrna = merged_eval['mrna_pvalue'] < 0.05
sig_prot = merged_eval['prot_p'] < 0.05
sig_either = sig_mrna | sig_prot
sig_both = sig_mrna & sig_prot

print(f'\n=== Concordance among significantly changed genes ===')
print(f'  Significant in mRNA only: {sig_mrna.sum() - sig_both.sum()}')
print(f'  Significant in protein only: {sig_prot.sum() - sig_both.sum()}')
print(f'  Significant in both: {sig_both.sum()}')

# Among genes significant in mRNA (direction from mRNA)
if sig_mrna.sum() > 0:
    sig_mrna_sub = merged_eval[sig_mrna]
    conc_sig_mrna = sig_mrna_sub['concordant'].sum()
    print(f'  Among sig-mRNA genes: {conc_sig_mrna}/{len(sig_mrna_sub)} concordant ({conc_sig_mrna/len(sig_mrna_sub)*100:.1f}%)')
    # ac4C vs non-targets within sig-mRNA
    ac_sig = sig_mrna_sub[sig_mrna_sub['is_ac4c_target']]
    non_sig = sig_mrna_sub[~sig_mrna_sub['is_ac4c_target']]
    if len(ac_sig) > 0 and len(non_sig) > 0:
        ac_cons = ac_sig['concordant'].sum()
        non_cons = non_sig['concordant'].sum()
        print(f'    ac4C targets: {ac_cons}/{len(ac_sig)} ({ac_cons/len(ac_sig)*100:.1f}%)')
        print(f'    non-targets: {non_cons}/{len(non_sig)} ({non_cons/len(non_sig)*100:.1f}%)')

# Among genes significant in both
if sig_both.sum() > 0:
    sig_both_sub = merged_eval[sig_both]
    conc_both = sig_both_sub['concordant'].sum()
    print(f'  Among sig-both genes: {conc_both}/{len(sig_both_sub)} concordant ({conc_both/len(sig_both_sub)*100:.1f}%)')
    ac_both = sig_both_sub[sig_both_sub['is_ac4c_target']]
    non_both = sig_both_sub[~sig_both_sub['is_ac4c_target']]
    if len(ac_both) > 0:
        print(f'    ac4C targets: {ac_both["concordant"].sum()}/{len(ac_both)} ({ac_both["concordant"].mean()*100:.1f}%)')
    if len(non_both) > 0:
        print(f'    non-targets: {non_both["concordant"].sum()}/{len(non_both)} ({non_both["concordant"].mean()*100:.1f}%)')

# ========== 5c. Concordance in upregulated proteins (synaptosome) ==========
# Many synaptosome proteins are up in AD; check if mRNA also goes up
prot_up = merged_eval['prot_log2FC'] > 0
prot_down = merged_eval['prot_log2FC'] < 0
print(f'\n=== Direction breakdown ===')
print(f'  Protein UP in AD: {prot_up.sum()}, concordant (mRNA also UP): {merged_eval[prot_up]["concordant"].sum()}')
print(f'  Protein DOWN in AD: {prot_down.sum()}, concordant (mRNA also DOWN): {merged_eval[prot_down]["concordant"].sum()}')

# ========== 5d. Direction consistency specifically for DOWN-ac4C-activator model ==========
# Hypothesis: NAT10 down -> ac4C down -> activator targets down in AD
# So ac4C activator targets should be DOWN in both mRNA and protein in AD
# Activator targets = those whose expression DECREASES when NAT10/ac4C is lost
# We can approximate by looking at genes DOWN in mRNA (assuming mRNA reflects ac4C effect)
mrna_down = merged_eval['mrna_log2FC'] < 0
print(f'\n=== Activator model: mRNA-down genes ===')
print(f'  Total mRNA-down: {mrna_down.sum()}')
ac_mrna_down = merged_eval[mrna_down & merged_eval['is_ac4c_target']]
non_mrna_down = merged_eval[mrna_down & ~merged_eval['is_ac4c_target']]
if len(ac_mrna_down) > 0:
    ac_prot_down = (ac_mrna_down['prot_log2FC'] < 0).sum()
    print(f'  ac4C targets mRNA-down: {len(ac_mrna_down)}, also protein-down: {ac_prot_down} ({ac_prot_down/len(ac_mrna_down)*100:.1f}%)')
if len(non_mrna_down) > 0:
    non_prot_down = (non_mrna_down['prot_log2FC'] < 0).sum()
    print(f'  non-targets mRNA-down: {len(non_mrna_down)}, also protein-down: {non_prot_down} ({non_prot_down/len(non_mrna_down)*100:.1f}%)')

# ========== 6. Specific ac4C target AD risk genes ==========

ad_risk_genes = {'APP', 'BIN1', 'CLU', 'SORL1', 'TREM2', 'MAPT', 'PSEN1', 'APOE', 'CD33', 'ABCA7', 'PICALM'}
both_ac4c_ad = ac4c & ad_risk_genes
ad_overlap_sub = merged_eval[merged_eval['gene'].isin(both_ac4c_ad)]
print(f'\n=== AD-risk ac4C target overlap (n={len(ad_overlap_sub)}) ===')
for _, row in ad_overlap_sub.iterrows():
    print(f'  {row["gene"]}: mRNA log2FC={row["mrna_log2FC"]:.4f}, protein log2FC={row["prot_log2FC"]:.4f}, '
          f'concordant={bool(row["concordant"])}')

# ========== 7. Save results ==========

# Main results table
results = []

# NAT10 row
results.append({
    'gene': 'NAT10',
    'category': 'writer_enzyme',
    'is_ac4c_target': False,
    'mRNA_log2FC_AD_vs_Ctrl': round(nat10_fc, 4),
        'mRNA_p': nat10_p,
        'mRNA_Braak_rho': round(nat10_rho, 4),
        'mRNA_Braak_p': nat10_rho_p,
    'protein_log2FC_AD_vs_HC': round(nat10_prot['log2FC'].values[0], 4) if len(nat10_prot) else np.nan,
    'protein_p': nat10_prot['p'].values[0] if len(nat10_prot) else np.nan,
    'protein_detected': len(nat10_prot) > 0,
    'concordant': np.nan,
    'note': 'writer enzyme; protein below detection in synaptosome MS'
})

# All evaluated genes
for _, row in merged_eval.iterrows():
    cat = 'ac4C_target' if row['is_ac4c_target'] else 'non_target'
    results.append({
        'gene': row['gene'],
        'category': cat,
        'is_ac4c_target': row['is_ac4c_target'],
        'mRNA_log2FC_AD_vs_Ctrl': round(row['mrna_log2FC'], 4),
        'mRNA_p': row['mrna_pvalue'],
        'mRNA_Braak_rho': np.nan,
        'mRNA_Braak_p': np.nan,
        'protein_log2FC_AD_vs_HC': round(row['prot_log2FC'], 4),
        'protein_p': row['prot_p'],
        'protein_detected': True,
        'concordant': bool(row['concordant']),
        'note': ''
    })

res_df = pd.DataFrame(results)
res_df.to_csv(os.path.join(OUT, 'rna_protein_concordance.csv'), index=False)
print(f'\nSaved: rna_protein_concordance.csv ({len(res_df)} rows)')

# Summary statistics table
summary = pd.DataFrame([
    {'metric': 'NAT10_mRNA_log2FC_AD_vs_Ctrl', 'value': nat10_fc},
    {'metric': 'NAT10_mRNA_p', 'value': nat10_p},
    {'metric': 'NAT10_mRNA_Braak_rho', 'value': nat10_rho},
    {'metric': 'NAT10_mRNA_Braak_p', 'value': nat10_rho_p},
    {'metric': 'NAT10_protein_detected_in_synaptosome', 'value': int(len(nat10_prot) > 0)},
    {'metric': 'NAT10_protein_log2FC', 'value': nat10_prot['log2FC'].values[0] if len(nat10_prot) else np.nan},
    {'metric': 'NAT10_protein_p', 'value': nat10_prot['p'].values[0] if len(nat10_prot) else np.nan},
    {'metric': 'n_genes_with_mRNA_and_protein', 'value': n_total},
    {'metric': 'overall_concordant', 'value': n_concordant},
    {'metric': 'overall_concordance_pct', 'value': round(concordance_pct, 1)},
    {'metric': 'overall_binomial_p_vs_50pct', 'value': binom_p},
    {'metric': 'overall_Pearson_r', 'value': round(r_pearson, 4)},
    {'metric': 'overall_Pearson_p', 'value': r_p_pearson},
    {'metric': 'overall_Spearman_rho', 'value': round(r_spearman, 4)},
    {'metric': 'overall_Spearman_p', 'value': r_p_spearman},
    {'metric': 'ac4C_target_n', 'value': ac4c_n},
    {'metric': 'ac4C_target_concordant', 'value': int(ac4c_cons)},
    {'metric': 'ac4C_target_concordance_pct', 'value': round(ac4c_pct, 1)},
    {'metric': 'ac4C_target_Spearman_rho', 'value': round(r_ac, 4) if not np.isnan(r_ac) else np.nan},
    {'metric': 'ac4C_target_Spearman_p', 'value': p_ac if not np.isnan(p_ac) else np.nan},
    {'metric': 'ac4C_target_mean_abs_log2FC', 'value': round(ac4c_mean_abs, 4)},
    {'metric': 'non_target_n', 'value': non_n},
    {'metric': 'non_target_concordant', 'value': int(non_cons)},
    {'metric': 'non_target_concordance_pct', 'value': round(non_pct, 1)},
    {'metric': 'non_target_Spearman_rho', 'value': round(r_non, 4)},
    {'metric': 'non_target_Spearman_p', 'value': p_non},
    {'metric': 'non_target_mean_abs_log2FC', 'value': round(non_mean_abs, 4)},
    {'metric': 'Fisher_OR_ac4C_vs_non', 'value': round(odds_ratio, 3)},
    {'metric': 'Fisher_p_ac4C_vs_non', 'value': fisher_p},
])
summary.to_csv(os.path.join(OUT, 'rna_protein_concordance_summary.csv'), index=False)
print(f'Saved: rna_protein_concordance_summary.csv')
print('\nDONE')
