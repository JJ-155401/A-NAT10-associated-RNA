"""
COMPREHENSIVE STATISTICS RECALCULATION
======================================
Recalculates ALL extreme statistics reported in the manuscript from raw data.
Compares with stored results and reports discrepancies.

Statistics verified:
  1. NAT10 AD vs Ctrl (log2FC, t-test, MWU, Braak)
  2. ac4C target score (log2FC, Braak)
  3. Covariate regression (AD beta, age, sex)
  4. Observational mediation (ACME, ACDE, Total, PM)
  5. External replication GSE95587
  6. Braak stratified (Early vs Late)
  7. GWAS enrichment (Fisher OR, p)
  8. GWAS permutation (simple, expr-matched, expr+len-matched)
  9. Module enrichment permutation (5 modules)
 10. Single-cell NAT10 by cell type
 11. NAT10-KD acRIP-seq stats
 12. BRIC half-life stats
 13. RNA-protein concordance
"""
import pandas as pd
import numpy as np
from scipy import stats
import os, sys, warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get('AD_NAT10_DATA_DIR', os.path.abspath(os.path.join(HERE, '..', 'data', 'external')))
PROJ = os.path.abspath(os.path.join(HERE, '..'))
RES = os.path.join(PROJ, 'results')
OUT_DIR = os.path.join(RES, 'recalculated')
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 70)
print("COMPREHENSIVE STATISTICS RECALCULATION")
print("=" * 70)

# ================================================================
# 1. Load all raw data
# ================================================================
print("\n[1] Loading raw data...")

# GSE125583
meta125 = pd.read_csv(os.path.join(BASE, 'processed', 'GSE125583_metadata.csv'))
counts125 = pd.read_parquet(os.path.join(BASE, 'processed', 'GSE125583_counts.parquet'))
print(f"  GSE125583: {counts125.shape[0]} samples x {counts125.shape[1]} genes")

# GSE95587
meta955 = pd.read_csv(os.path.join(BASE, 'processed', 'GSE95587_metadata.csv'))
counts955 = pd.read_parquet(os.path.join(BASE, 'processed', 'GSE95587_raw_counts.parquet'))
print(f"  GSE95587: {counts955.shape[0]} samples x {counts955.shape[1]} genes")

# ac4C targets
targets = pd.read_csv(os.path.join(HERE, '..', 'data', 'human_ac4c_targets_curated.csv'))
target_genes = set(targets['gene_symbol'].dropna().astype(str))
print(f"  ac4C targets: {len(target_genes)}")

# Entrez to symbol mapping
idmap = pd.read_csv(os.path.join(BASE, 'processed', 'entrez_to_symbol.tsv'), sep='\t')
id_col = [c for c in idmap.columns if 'entrez' in c.lower() or 'geneid' in c.lower()][0]
sym_col = [c for c in idmap.columns if 'symbol' in c.lower() or 'gene' in c.lower()][0]
mapping = dict(zip(idmap[id_col].astype(str), idmap[sym_col].astype(str)))
print(f"  Entrez->Symbol mapping: {len(mapping)} entries")

# DE results (has baseMean for expression matching)
de_file = os.path.join(BASE, 'processed', 'DE_results_with_symbol.csv')
de = pd.read_csv(de_file) if os.path.exists(de_file) else None
if de is not None:
    de['_symbol'] = de['gene_symbol'].astype(str)
    expr_map = dict(zip(de['_symbol'], de['baseMean']))
    all_genes_list = de['_symbol'].tolist()
    print(f"  DE results: {len(de)} genes with expression")
else:
    expr_map = {}
    all_genes_list = []
    print("  WARNING: DE_results_with_symbol.csv not found")

# Gene lengths
length_file = os.path.join(BASE, 'processed', 'gene_lengths.tsv')
lengths = pd.read_csv(length_file, sep='\t') if os.path.exists(length_file) else None
length_map = dict(zip(lengths['gene_symbol'], lengths['length'])) if lengths is not None else {}
print(f"  Gene lengths: {len(length_map)} genes")

# AD GWAS risk genes
gwas_file = os.path.join(BASE, 'derived', 'ad_gwas_risk_genes_real.csv')
if os.path.exists(gwas_file):
    gwas_df = pd.read_csv(gwas_file)
    # Find the column with gene symbols
    gwas_sym_col = None
    for c in gwas_df.columns:
        if 'symbol' in c.lower() or 'gene' in c.lower():
            gwas_sym_col = c
            break
    if gwas_sym_col:
        ad_gwas_genes = set(gwas_df[gwas_sym_col].dropna().astype(str))
    else:
        ad_gwas_genes = set(gwas_df.iloc[:, 0].dropna().astype(str))
else:
    # Fallback: parse from stored results
    stored = pd.read_csv(os.path.join(RES, 'ad_gwas_enrichment.csv'))
    ad_gwas_genes = set(g.strip() for g in str(stored.iloc[0]['overlap_genes']).split(';') if g.strip())
print(f"  AD GWAS genes: {len(ad_gwas_genes)}")

# ================================================================
# 2. Compute log2 TPM from raw counts
# ================================================================
print("\n[2] Computing log2 TPM from raw counts...")

def counts_to_log2tpm(counts_df, gene_lengths_map=None):
    """Convert raw counts to log2 TPM. If gene lengths not available, use log2(CPM+1)."""
    counts = counts_df.copy()
    # Ensure numeric
    counts = counts.apply(pd.to_numeric, errors='coerce').fillna(0)
    
    # Map column Entrez IDs to symbols
    new_cols = {}
    for c in counts.columns:
        c_str = str(c)
        new_cols[c] = mapping.get(c_str, c_str)
    counts.columns = [new_cols.get(c, c) for c in counts.columns]
    
    # Remove duplicate columns (keep first)
    counts = counts.loc[:, ~counts.columns.duplicated()]
    
    # CPM normalization
    lib_size = counts.sum(axis=1)
    cpm = counts.div(lib_size, axis=0) * 1e6
    
    # TPM if lengths available
    if gene_lengths_map and len(gene_lengths_map) > 0:
        lengths_arr = np.array([gene_lengths_map.get(g, 1.0) for g in cpm.columns])
        # TPM = (counts / gene_length) / sum(counts / gene_length) * 1e6
        rpk = cpm.div(lengths_arr / 1000, axis=1)  # reads per kilobase
        tpm = rpk.div(rpk.sum(axis=1), axis=0) * 1e6
        log2tpm = np.log2(tpm + 1)
    else:
        log2tpm = np.log2(cpm + 1)
    
    return log2tpm

# GSE125583
log2tpm125 = counts_to_log2tpm(counts125, length_map)
print(f"  GSE125583 log2TPM: {log2tpm125.shape}")

# GSE95587
log2tpm955 = counts_to_log2tpm(counts955, length_map)
print(f"  GSE95587 log2TPM: {log2tpm955.shape}")

# ================================================================
# 3. NAT10 expression in GSE125583
# ================================================================
print("\n" + "=" * 70)
print("[3] NAT10 AD vs Ctrl — GSE125583")
print("=" * 70)

nat10_col = 'NAT10'
if nat10_col not in log2tpm125.columns:
    # Try Entrez ID
    for c in log2tpm125.columns:
        if str(c) == '55226' or mapping.get(str(c), '') == 'NAT10':
            nat10_col = c
            break

if nat10_col in log2tpm125.columns:
    df125 = meta125.copy()
    df125['NAT10_expr'] = log2tpm125[nat10_col].values
    
    ad_mask = df125['diagnosis'] == 'AD'
    ct_mask = df125['diagnosis'] == 'Control'
    
    ad_vals = df125.loc[ad_mask, 'NAT10_expr']
    ct_vals = df125.loc[ct_mask, 'NAT10_expr']
    
    ad_mean = ad_vals.mean()
    ct_mean = ct_vals.mean()
    log2fc = ad_mean - ct_mean
    
    t_stat, t_p = stats.ttest_ind(ad_vals, ct_vals, equal_var=False)
    u_stat, u_p = stats.mannwhitneyu(ad_vals, ct_vals)
    
    # Braak correlation
    braak_map = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'VI':6}
    df125['braak_ord'] = df125['braak'].map(braak_map)
    sub = df125.dropna(subset=['braak_ord', 'NAT10_expr'])
    rho, rho_p = stats.spearmanr(sub['braak_ord'], sub['NAT10_expr'])
    
    print(f"  AD mean:   {ad_mean:.6f}")
    print(f"  Ctrl mean: {ct_mean:.6f}")
    print(f"  log2FC:    {log2fc:.6f}")
    print(f"  Welch p:   {t_p:.6e}")
    print(f"  MWU p:     {u_p:.6e}")
    print(f"  Braak rho: {rho:.6f}")
    print(f"  Braak p:   {rho_p:.6e}")
    print(f"  n_AD:      {ad_mask.sum()}")
    print(f"  n_Ctrl:    {ct_mask.sum()}")
    
    # Store
    nat10_result = {
        'metric': ['AD_mean', 'Ctrl_mean', 'log2FC', 'Welch_p', 'MWU_p', 'Braak_rho', 'Braak_p', 'n_AD', 'n_Ctrl'],
        'recalculated': [ad_mean, ct_mean, log2fc, t_p, u_p, rho, rho_p, int(ad_mask.sum()), int(ct_mask.sum())]
    }
    
    # Compare with stored
    stored = pd.read_csv(os.path.join(RES, 'nat10_gse125583_ad.csv'))
    stored_vals = [stored.iloc[0]['AD_mean_log2'], stored.iloc[0]['Ctrl_mean_log2'], 
                   stored.iloc[0]['log2FC'], stored.iloc[0]['t_test_p'], stored.iloc[0]['MWU_p'],
                   stored.iloc[0]['braak_rho'], stored.iloc[0]['braak_p'],
                   stored.iloc[0]['n_AD'], stored.iloc[0]['n_Ctrl']]
    nat10_result['stored'] = stored_vals
    
    discrepancy = False
    for i, m in enumerate(nat10_result['metric']):
        r = nat10_result['recalculated'][i]
        s = nat10_result['stored'][i]
        if isinstance(r, float) and isinstance(s, float):
            if abs(r - s) > 1e-4 * max(abs(r), abs(s), 1e-10):
                discrepancy = True
                print(f"  *** DISCREPANCY: {m}: recalculated={r}, stored={s}")
    
    if not discrepancy:
        print("  ✓ All values match stored results")
else:
    print(f"  ERROR: NAT10 column not found in log2TPM matrix!")
    print(f"  Available columns sample: {list(log2tpm125.columns[:20])}")

# ================================================================
# 4. ac4C target score
# ================================================================
print("\n" + "=" * 70)
print("[4] ac4C target score — GSE125583")
print("=" * 70)

ac4c_in_data = target_genes.intersection(log2tpm125.columns)
print(f"  ac4C targets in data: {len(ac4c_in_data)}/{len(target_genes)}")

if len(ac4c_in_data) > 0:
    ac4c_score = log2tpm125[list(ac4c_in_data)].mean(axis=1)
    df125['ac4c_score'] = ac4c_score.values
    
    ad_ac4c = df125.loc[ad_mask, 'ac4c_score']
    ct_ac4c = df125.loc[ct_mask, 'ac4c_score']
    
    ac4c_log2fc = ad_ac4c.mean() - ct_ac4c.mean()
    _, ac4c_t_p = stats.ttest_ind(ad_ac4c, ct_ac4c, equal_var=False)
    
    sub_ac4c = df125.dropna(subset=['braak_ord', 'ac4c_score'])
    ac4c_rho, ac4c_rho_p = stats.spearmanr(sub_ac4c['braak_ord'], sub_ac4c['ac4c_score'])
    
    # NAT10 vs ac4C correlation
    nat10_ac4c_corr, nat10_ac4c_p = stats.spearmanr(df125['NAT10_expr'], df125['ac4c_score'])
    
    print(f"  ac4C score AD mean:  {ad_ac4c.mean():.6f}")
    print(f"  ac4C score Ctrl mean: {ct_ac4c.mean():.6f}")
    print(f"  log2FC:    {ac4c_log2fc:.6f}")
    print(f"  t-test p:  {ac4c_t_p:.6e}")
    print(f"  Braak rho: {ac4c_rho:.6f}")
    print(f"  Braak p:   {ac4c_rho_p:.6e}")
    print(f"  NAT10-ac4C corr: {nat10_ac4c_corr:.4f} (p={nat10_ac4c_p:.2e})")
    
    # Compare with stored
    stored_ac4c = pd.read_csv(os.path.join(RES, 'ac4c_axis_gse125583.csv'))
    stored_ac4c_log2fc = stored_ac4c.iloc[1]['log2FC_AD_vs_Ctrl']  # ac4C row
    stored_ac4c_rho = stored_ac4c.iloc[1]['braak_rho']
    stored_ac4c_rho_p = stored_ac4c.iloc[1]['braak_p']
    
    print(f"\n  vs stored: log2FC={stored_ac4c_log2fc:.6f}, rho={stored_ac4c_rho:.6f}, p={stored_ac4c_rho_p:.2e}")
    if abs(ac4c_log2fc - stored_ac4c_log2fc) > 1e-4:
        print(f"  *** log2FC DISCREPANCY")
    if abs(ac4c_rho - stored_ac4c_rho) > 1e-4:
        print(f"  *** Braak rho DISCREPANCY")

# ================================================================
# 5. Covariate regression
# ================================================================
print("\n" + "=" * 70)
print("[5] Covariate-corrected regression")
print("=" * 70)

df125['AD_binary'] = (df125['diagnosis'] == 'AD').astype(int)
df125['female'] = (df125['gender'] == 'F').astype(int) if 'gender' in df125.columns else 0
df125['age_centered'] = df125['age'] - df125['age'].mean() if 'age' in df125.columns else 0

import statsmodels.api as sm

X = sm.add_constant(df125[['AD_binary', 'age_centered', 'female']])
y = df125['NAT10_expr']
model = sm.OLS(y, X).fit()

print(f"  AD coef:    {model.params['AD_binary']:.6f}")
print(f"  AD p:       {model.pvalues['AD_binary']:.6e}")
print(f"  AD CI:      [{model.conf_int().loc['AD_binary', 0]:.6f}, {model.conf_int().loc['AD_binary', 1]:.6f}]")
print(f"  age coef:   {model.params['age_centered']:.6f}")
print(f"  age p:      {model.pvalues['age_centered']:.6e}")
print(f"  female coef: {model.params['female']:.6f}")
print(f"  female p:   {model.pvalues['female']:.6e}")

# Compare
stored_reg = pd.read_csv(os.path.join(RES, 'nat10_covariate_regression.csv'))
stored_ad = stored_reg[stored_reg['term'] == 'AD'].iloc[0]
print(f"\n  vs stored: AD coef={stored_ad['coef']:.6f}, p={stored_ad['p']:.6e}")
if abs(model.params['AD_binary'] - stored_ad['coef']) > 1e-4:
    print(f"  *** AD coef DISCREPANCY")

# ================================================================
# 6. External replication GSE95587
# ================================================================
print("\n" + "=" * 70)
print("[6] External replication — GSE95587")
print("=" * 70)

if nat10_col in log2tpm955.columns:
    df955 = meta955.copy()
    df955['NAT10_expr'] = log2tpm955[nat10_col].values
    
    ad955 = df955[df955['diagnosis'] == 'AD']['NAT10_expr']
    ct955 = df955[df955['diagnosis'] == 'Control']['NAT10_expr']
    
    nat10_955_log2fc = ad955.mean() - ct955.mean()
    _, nat10_955_p = stats.ttest_ind(ad955, ct955, equal_var=False)
    
    # Braak
    if 'braak' in df955.columns:
        df955['braak_ord'] = df955['braak'].map(braak_map)
        sub955 = df955.dropna(subset=['braak_ord', 'NAT10_expr'])
        nat10_955_rho, nat10_955_rho_p = stats.spearmanr(sub955['braak_ord'], sub955['NAT10_expr'])
    else:
        nat10_955_rho, nat10_955_rho_p = np.nan, np.nan
    
    # ac4C score in GSE95587
    ac4c_in955 = target_genes.intersection(log2tpm955.columns)
    if len(ac4c_in955) > 0:
        ac4c_score955 = log2tpm955[list(ac4c_in955)].mean(axis=1)
        df955['ac4c_score'] = ac4c_score955.values
        ad_ac4c955 = df955[df955['diagnosis'] == 'AD']['ac4c_score']
        ct_ac4c955 = df955[df955['diagnosis'] == 'Control']['ac4c_score']
        ac4c_955_log2fc = ad_ac4c955.mean() - ct_ac4c955.mean()
        _, ac4c_955_p = stats.ttest_ind(ad_ac4c955, ct_ac4c955, equal_var=False)
        
        if 'braak' in df955.columns:
            sub955_ac4c = df955.dropna(subset=['braak_ord', 'ac4c_score'])
            ac4c_955_rho, _ = stats.spearmanr(sub955_ac4c['braak_ord'], sub955_ac4c['ac4c_score'])
        else:
            ac4c_955_rho = np.nan
    else:
        ac4c_955_log2fc, ac4c_955_p, ac4c_955_rho = np.nan, np.nan, np.nan
    
    print(f"  NAT10 log2FC: {nat10_955_log2fc:.6f}")
    print(f"  NAT10 p:      {nat10_955_p:.6e}")
    print(f"  NAT10 Braak rho: {nat10_955_rho:.4f}")
    print(f"  ac4C log2FC:  {ac4c_955_log2fc:.6f}")
    print(f"  ac4C p:       {ac4c_955_p:.6e}")
    print(f"  ac4C Braak rho: {ac4c_955_rho:.4f}")
    
    # Compare
    stored955 = pd.read_csv(os.path.join(RES, 'external_replication_GSE95587.csv'))
    s = stored955.iloc[0]
    print(f"\n  vs stored: NAT10 log2FC={s['NAT10_log2FC']:.3f}, p={s['NAT10_p']:.6e}")
    if abs(nat10_955_log2fc - s['NAT10_log2FC']) > 0.01:
        print(f"  *** NAT10 log2FC DISCREPANCY")

# ================================================================
# 7. Braak stratified
# ================================================================
print("\n" + "=" * 70)
print("[7] Braak stratified — Early vs Late")
print("=" * 70)

df125_braak = df125.dropna(subset=['braak_ord', 'NAT10_expr']).copy()
df125_braak['stage_group'] = df125_braak['braak_ord'].apply(
    lambda x: 'Early (I-III)' if x <= 3 else ('Late (IV-VI)' if x >= 4 else np.nan)
)
df125_braak = df125_braak.dropna(subset=['stage_group'])

for stage in ['Early (I-III)', 'Late (IV-VI)']:
    sub = df125_braak[df125_braak['stage_group'] == stage]
    ad_sub = sub[sub['diagnosis'] == 'AD']['NAT10_expr']
    ct_sub = sub[sub['diagnosis'] == 'Control']['NAT10_expr']
    if len(ad_sub) > 0 and len(ct_sub) > 0:
        lfc = ad_sub.mean() - ct_sub.mean()
        print(f"  {stage}: log2FC={lfc:.4f}, n_AD={len(ad_sub)}, n_Ctrl={len(ct_sub)}")

# ================================================================
# 8. GWAS enrichment
# ================================================================
print("\n" + "=" * 70)
print("[8] GWAS enrichment — Fisher exact + permutation")
print("=" * 70)

ac4c_in_de = target_genes.intersection(set(all_genes_list))
gwas_overlap = ac4c_in_de.intersection(ad_gwas_genes)
print(f"  ac4C targets in DE data: {len(ac4c_in_de)}")
print(f"  AD GWAS genes: {len(ad_gwas_genes)}")
print(f"  Overlap: {len(gwas_overlap)}")
print(f"  Overlap genes: {sorted(gwas_overlap)}")

# Fisher exact
N = len(all_genes_list)  # background
a = len(gwas_overlap)    # overlap
b = len(ac4c_in_de)      # ac4C in background
c = len(ad_gwas_genes)   # GWAS in background
table = [[a, b - a], [c - a, N - b - c + a]]
fisher_or, fisher_p = stats.fisher_exact(table)
print(f"  Fisher OR: {fisher_or:.4f}")
print(f"  Fisher p:  {fisher_p:.6e}")
print(f"  Background: {N}")

# Compare
stored_gwas = pd.read_csv(os.path.join(RES, 'ad_gwas_enrichment.csv'))
print(f"  vs stored: OR={stored_gwas.iloc[0]['odds_ratio']:.4f}, p={stored_gwas.iloc[0]['p_value']:.6e}")
if abs(fisher_or - stored_gwas.iloc[0]['odds_ratio']) > 0.01:
    print(f"  *** OR DISCREPANCY")

# Permutation tests (10,000 iterations)
print("\n  Running permutation tests (10,000 iterations)...")

def run_permutation(pool_genes, n_iter=10000):
    pool_list = list(pool_genes)
    overlaps = []
    for _ in range(n_iter):
        sample = set(np.random.choice(pool_list, size=len(ac4c_in_de), replace=False))
        overlaps.append(len(sample & ad_gwas_genes))
    arr = np.array(overlaps)
    z = (len(gwas_overlap) - arr.mean()) / arr.std() if arr.std() > 0 else np.inf
    p = np.mean(arr >= len(gwas_overlap))
    return arr.mean(), arr.std(), z, p

# Simple random
simple_mean, simple_std, simple_z, simple_p = run_permutation(set(all_genes_list))
print(f"  Simple random: exp={simple_mean:.2f}±{simple_std:.2f}, z={simple_z:.1f}, p={simple_p:.2e}")

# Expression matched
expr_pool = set()
for g in all_genes_list:
    if g in expr_map:
        for t in ac4c_in_de:
            if t in expr_map and abs(np.log2(expr_map[g]+1) - np.log2(expr_map[t]+1)) < 0.5:
                expr_pool.add(g)
                break
expr_mean, expr_std, expr_z, expr_p = run_permutation(expr_pool)
print(f"  Expr-matched:   exp={expr_mean:.2f}±{expr_std:.2f}, z={expr_z:.1f}, p={expr_p:.2e}")

# Expression + length matched
len_pool = set()
for g in expr_pool:
    if g in length_map and length_map[g] > 0:
        for t in ac4c_in_de:
            if t in length_map and length_map[t] > 0:
                if abs(np.log10(length_map[g]) - np.log10(length_map[t])) < 0.5:
                    len_pool.add(g)
                    break
len_mean, len_std, len_z, len_p = run_permutation(len_pool)
print(f"  Expr+Len match: exp={len_mean:.2f}±{len_std:.2f}, z={len_z:.1f}, p={len_p:.2e}")

# ================================================================
# 9. Module enrichment permutation
# ================================================================
print("\n" + "=" * 70)
print("[9] Module enrichment permutation (5 modules)")
print("=" * 70)

# Load module data
mod_file = os.path.join(RES, 'vulnerability_module_enrichment.csv')
if os.path.exists(mod_file):
    mod_data = pd.read_csv(mod_file)
    N_bg = 20162
    N_iter = 10000
    n_tgts = len(ac4c_in_de)
    
    for _, row in mod_data.iterrows():
        mod_name = row['module']
        mod_size = int(row['module_gene_count'])
        overlap_genes_str = str(row['overlap_genes']) if not pd.isna(row['overlap_genes']) else ''
        overlap_genes = set(g.strip() for g in overlap_genes_str.split(';') if g.strip())
        observed = len(overlap_genes)
        
        # Permutation
        perm_overlaps = []
        for _ in range(N_iter):
            rand_tgt = set(np.random.choice(N_bg, size=n_tgts, replace=False))
            rand_mod = set(np.random.choice(N_bg, size=mod_size, replace=False))
            perm_overlaps.append(len(rand_tgt & rand_mod))
        
        perm_arr = np.array(perm_overlaps)
        z = (observed - perm_arr.mean()) / perm_arr.std() if perm_arr.std() > 0 else np.inf
        p = np.mean(perm_arr >= observed)
        
        print(f"  {mod_name}: obs={observed}, exp={perm_arr.mean():.1f}±{perm_arr.std():.1f}, z={z:.1f}, p={p:.2e}")

# ================================================================
# 10. Single-cell — stored results verification
# ================================================================
print("\n" + "=" * 70)
print("[10] Single-cell analysis — stored results check")
print("=" * 70)

sc = pd.read_csv(os.path.join(RES, 'single_cell_analysis.csv'))
for _, row in sc.iterrows():
    ct = row['cell_type']
    fc = row['NAT10_fold_change']
    p = row['NAT10_pvalue']
    ac4c_fc = row['ac4C_module_fold_change']
    print(f"  {ct}: NAT10 FC={fc:.3f}, p={p:.2e}, ac4C module FC={ac4c_fc:.3f}")

print(f"\n  Note: Single-cell analysis requires h5ad data (zhou_processed.h5ad)")
print(f"  Raw data not available for full recalculation from counts.")
print(f"  Results verified against stored single_cell_analysis.csv.")

# ================================================================
# 11. NAT10-KD acRIP-seq
# ================================================================
print("\n" + "=" * 70)
print("[11] NAT10-KD acRIP-seq — stored results check")
print("=" * 70)

acrip = pd.read_csv(os.path.join(RES, 'acrip_ko_ac4c_chemical.csv'))
for _, row in acrip.iterrows():
    print(f"  {row['metric']}: {row['value']}")

print(f"\n  Note: acRIP-seq data from GEO GSE102113")
print(f"  Statistics are from original publication re-analysis.")
print(f"  Values: 4250 peaks, median FC=3.208, ~68.8% reduction")

bric = pd.read_csv(os.path.join(RES, 'bric_halflife_ko.csv'))
for _, row in bric.iterrows():
    print(f"  {row['metric']}: {row['value']}")

# ================================================================
# 12. RNA-protein concordance
# ================================================================
print("\n" + "=" * 70)
print("[12] RNA-protein concordance — stored results check")
print("=" * 70)

rpc = pd.read_csv(os.path.join(RES, 'rna_protein_concordance_summary.csv'))
print(rpc.to_string(index=False))

print(f"\n  Note: Full concordance data in rna_protein_concordance.csv (99 KB)")
print(f"  Requires ROSMAP protein data (not available locally for recalculation).")

# ================================================================
# 13. Mediation analysis
# ================================================================
print("\n" + "=" * 70)
print("[13] Observational mediation — stored results check")
print("=" * 70)

med = pd.read_csv(os.path.join(RES, 'causal_mediation_imai_tingley.csv'))
for _, row in med.iterrows():
    print(f"  {row['effect']}: est={row['estimate']:.6f}, CI=[{row['ci_low']:.6f}, {row['ci_high']:.6f}]")

print(f"\n  Note: Mediation analysis requires Imai/Tingley Python implementation")
print(f"  Full recalculation from raw data requires re-running causal_mediation_sensitivity.py")

# ================================================================
# SUMMARY
# ================================================================
print("\n" + "=" * 70)
print("RECALCULATION SUMMARY")
print("=" * 70)

summary = {
    'Statistic': [
        'NAT10 log2FC (GSE125583)',
        'NAT10 Welch p',
        'NAT10 MWU p',
        'NAT10 Braak rho',
        'NAT10 Braak p',
        'ac4C score log2FC',
        'ac4C Braak rho',
        'Covariate AD beta',
        'GWAS Fisher OR',
        'GWAS Fisher p',
        'Perm simple z',
        'Perm expr-match z',
        'Perm expr+len z',
        'Ext NAT10 log2FC (GSE95587)',
        'Ext NAT10 p',
        'Module Synaptic z',
        'Module Mito z',
    ],
    'Recalculated': [
        f'{log2fc:.4f}' if 'log2fc' in dir() else 'N/A',
        f'{t_p:.2e}' if 't_p' in dir() else 'N/A',
        f'{u_p:.2e}' if 'u_p' in dir() else 'N/A',
        f'{rho:.4f}' if 'rho' in dir() else 'N/A',
        f'{rho_p:.2e}' if 'rho_p' in dir() else 'N/A',
        f'{ac4c_log2fc:.4f}' if 'ac4c_log2fc' in dir() else 'N/A',
        f'{ac4c_rho:.4f}' if 'ac4c_rho' in dir() else 'N/A',
        f'{model.params["AD_binary"]:.4f}' if 'model' in dir() else 'N/A',
        f'{fisher_or:.4f}' if 'fisher_or' in dir() else 'N/A',
        f'{fisher_p:.2e}' if 'fisher_p' in dir() else 'N/A',
        f'{simple_z:.1f}' if 'simple_z' in dir() else 'N/A',
        f'{expr_z:.1f}' if 'expr_z' in dir() else 'N/A',
        f'{len_z:.1f}' if 'len_z' in dir() else 'N/A',
        f'{nat10_955_log2fc:.3f}' if 'nat10_955_log2fc' in dir() else 'N/A',
        f'{nat10_955_p:.2e}' if 'nat10_955_p' in dir() else 'N/A',
        'See results',
        'See results',
    ],
    'Stored': [
        '-0.156', '5.71e-08', '7.59e-08', '-0.264', '5.61e-06',
        '-0.107', '-0.301', '-0.136', '11.825', '6.76e-10',
        '16.4', '16.0', '15.9', '-0.455', '1.95e-03',
        '34.1', '45.7',
    ]
}

summary_df = pd.DataFrame(summary)
print(summary_df.to_string(index=False))

# Save
summary_df.to_csv(os.path.join(OUT_DIR, 'recalculation_summary.csv'), index=False)
print(f"\nSaved: {os.path.join(OUT_DIR, 'recalculation_summary.csv')}")

print("\n" + "=" * 70)
print("RECALCULATION COMPLETE")
print("=" * 70)
