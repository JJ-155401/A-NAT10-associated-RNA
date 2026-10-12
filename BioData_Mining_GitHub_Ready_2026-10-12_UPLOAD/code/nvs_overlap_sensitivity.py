"""
Problem 3 fix: NVS overlap sensitivity - remove ac4C∩NVS overlapping genes,
recalculate ac4C module vs NVS correlation in single-cell data.
"""
import scanpy as sc
import numpy as np
import pandas as pd
from scipy import stats
import os
sc.settings.verbosity = 0

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get('AD_NAT10_DATA_DIR', os.path.join(HERE, '..', 'data', 'external'))
results_dir = os.path.join(HERE, '..', 'results')
adata = sc.read_h5ad(os.path.join(BASE, 'scrna', 'zhou', 'zhou_merged.h5ad'))
print(f"Loaded: {adata.shape}")

# ============================================================
# 1. Cell type annotation (same as before)
# ============================================================
marker_genes = {
    'Excitatory neuron': ['SLC17A7', 'CAMK2A', 'SATB2', 'NRGN', 'SLC17A6'],
    'Inhibitory neuron': ['GAD1', 'GAD2', 'SLC32A1', 'SST', 'PVALB', 'VIP', 'LHX6'],
    'Astrocyte': ['AQP4', 'GFAP', 'GJA1', 'SLC1A2', 'SLC1A3', 'ALDH1L1'],
    'Oligodendrocyte': ['MBP', 'MOG', 'PLP1', 'OLIG2', 'SLC1A3', 'CNNM2'],
    'OPC': ['PDGFRA', 'CSPG4', 'OLIG2', 'SOX10'],
    'Microglia': ['CX3CR1', 'P2RY12', 'TMEM119', 'AIF1', 'CSF1R'],
    'Endothelial': ['PECAM1', 'FLT1', 'CLDN5', 'RGS5', 'ABCB1'],
}

for ct, genes in marker_genes.items():
    genes_in_data = [g for g in genes if g in adata.var_names]
    if len(genes_in_data) >= 2:
        sc.tl.score_genes(adata, gene_list=genes_in_data, score_name=f'score_{ct}', use_raw=False)

score_cols = [f'score_{ct}' for ct in marker_genes.keys() if f'score_{ct}' in adata.obs.columns]
score_matrix = adata.obs[score_cols].values
max_idx = np.argmax(score_matrix, axis=1)
ct_names = [ct for ct in marker_genes.keys() if f'score_{ct}' in adata.obs.columns]
adata.obs['cell_type'] = [ct_names[i] for i in max_idx]

# ============================================================
# 2. Load gene sets
# ============================================================
# ac4C targets
targets = pd.read_csv(os.path.join(HERE, '..', 'data', 'human_ac4c_targets_curated.csv'))
idmap = pd.read_csv(os.path.join(BASE, 'processed', 'entrez_to_symbol.tsv'), sep='\t')
id_col = [c for c in idmap.columns if 'entrez' in c.lower() or 'geneid' in c.lower()][0]
sym_col = [c for c in idmap.columns if 'symbol' in c.lower() or 'gene' in c.lower()][0]
mapping = dict(zip(idmap[id_col].astype(str), idmap[sym_col].astype(str)))
ac4c_genes_raw = targets['gene_symbol'].dropna().astype(str).tolist()
ac4c_genes = set([mapping.get(g, g) for g in ac4c_genes_raw])
ac4c_in_data = [g for g in ac4c_genes if g in adata.var_names]
print(f"ac4C targets in data: {len(ac4c_in_data)}/{len(ac4c_genes)}")

# NVS gene set (same as NVS script)
nvs_genes_all = [
    'SYN1','SYP','SNAP25','STX1A','VAMP1','VAMP2','SYT1','DLG4','SHANK1','SHANK3',
    'GRIA1','GRIA2','GRIN1','GRIN2A','GRIN2B','CAMK2A','ARC','BDNF','NTRK2',
    'APP','MAPT','PSEN1','PSEN2','APOE','BIN1','CLU','PICALM','SORL1','ABCA7',
    'TREM2','CD33','MEF2C','NFKB1','ACE','ADAM10','GRN','BACE1',
    'ATP5F1A','COX5A','NDUFS1','NDUFS3','SDHA','UQCRC1','CS','IDH2','FH',
    'SOD2','PRDX3','Pink1',
]
nvs_set = set(nvs_genes_all)
nvs_in_data = [g for g in nvs_genes_all if g in adata.var_names]
print(f"NVS genes in data: {len(nvs_in_data)}/{len(nvs_genes_all)}")

# ============================================================
# 3. Overlap analysis
# ============================================================
overlap_genes = ac4c_genes & nvs_set
overlap_in_data = [g for g in overlap_genes if g in adata.var_names]
print(f"\n=== Overlap: ac4C ∩ NVS ===")
print(f"Overlap (symbol space): {len(overlap_genes)}")
print(f"Overlap (in data): {len(overlap_in_data)}")
print(f"  {sorted(overlap_in_data)}")

# ============================================================
# 4. Compute scores: full vs non-overlapping
# ============================================================
# Full gene sets
sc.tl.score_genes(adata, gene_list=ac4c_in_data, score_name='ac4C_full', use_raw=False)
sc.tl.score_genes(adata, gene_list=nvs_in_data, score_name='NVS_full', use_raw=False)

# Non-overlapping sets
ac4c_nonoverlap = [g for g in ac4c_in_data if g not in overlap_in_data]
nvs_nonoverlap = [g for g in nvs_in_data if g not in overlap_in_data]
print(f"\nac4C non-overlap genes: {len(ac4c_nonoverlap)}")
print(f"NVS non-overlap genes: {len(nvs_nonoverlap)}")

sc.tl.score_genes(adata, gene_list=ac4c_nonoverlap, score_name='ac4C_no_overlap', use_raw=False)
sc.tl.score_genes(adata, gene_list=nvs_nonoverlap, score_name='NVS_no_overlap', use_raw=False)

# ============================================================
# 5. Correlation: full vs non-overlapping
# ============================================================
print("\n" + "="*60)
print("CORRELATION: ac4C module vs NVS")
print("="*60)

results = []
# Full gene sets
for ct in adata.obs['cell_type'].unique():
    mask = adata.obs['cell_type'] == ct
    sub = adata.obs[mask]
    if len(sub) < 100:
        continue
    
    r_full, p_full = stats.spearmanr(sub['ac4C_full'], sub['NVS_full'])
    r_nov, p_nov = stats.spearmanr(sub['ac4C_no_overlap'], sub['NVS_no_overlap'])
    
    results.append({
        'cell_type': ct,
        'n_cells': int(mask.sum()),
        'r_full': r_full,
        'p_full': p_full,
        'r_no_overlap': r_nov,
        'p_no_overlap': p_nov,
    })
    print(f"  {ct} (n={mask.sum()}):")
    print(f"    Full sets:        r={r_full:.3f}, p={p_full:.1e}")
    print(f"    Non-overlap only: r={r_nov:.3f}, p={p_nov:.1e}")

# Overall
r_full_all, p_full_all = stats.spearmanr(adata.obs['ac4C_full'], adata.obs['NVS_full'])
r_nov_all, p_nov_all = stats.spearmanr(adata.obs['ac4C_no_overlap'], adata.obs['NVS_no_overlap'])
print(f"\n  Overall (all cells):")
print(f"    Full sets:        r={r_full_all:.3f}, p={p_full_all:.1e}")
print(f"    Non-overlap only: r={r_nov_all:.3f}, p={p_nov_all:.1e}")

results_df = pd.DataFrame(results)
results_df.to_csv(os.path.join(results_dir, 'nvs_ac4c_sensitivity_overlap.csv'), index=False)

print("\n=== DONE ===")
