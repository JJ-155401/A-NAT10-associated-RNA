"""
Fig4 specificity control (Reviewer must-do #5):
Does ac4C-target enrichment with AD GWAS generalize to ANY gene set,
or is it specific to AD genetic architecture?
We compare the observed enrichment (OR=11.8, 13/303 vs 88 AD GWAS genes)
against:
  (A) 10,000 random 88-gene sets drawn from the expressed background (20,162 genes)
  (B) 10,000 random 88-gene sets drawn from a curated BRAIN-DISEASE gene pool
      (neurodegeneration / psychiatry / neurodevelopmental genes), to show the
      enrichment is NOT a generic property of brain-disease risk gene sets.
"""
import pandas as pd
import numpy as np
from scipy import stats
import os

np.random.seed(7)
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get('AD_NAT10_DATA_DIR', os.path.join(HERE, '..', 'data', 'external'))
results_dir = os.path.join(HERE, '..', 'results')

gwas = pd.read_csv(os.path.join(results_dir, 'ad_gwas_enrichment.csv')).iloc[0]
ac4c_n = int(gwas['ac4c_targets_n'])      # 303
gwas_n = int(gwas['ad_gwas_genes_n'])     # 88
bg = int(gwas['background_n'])            # 20162
overlap_genes = set(g.strip() for g in str(gwas['overlap_genes']).split(';') if g.strip())
ac4c_set = set(overlap_genes)  # placeholder; need full 303-set

# Load full ac4C target set (303) and the matching background universe (20,162)
targets = pd.read_csv(os.path.join(HERE, '..', 'data', 'human_ac4c_targets_curated.csv'))
idmap = pd.read_csv(os.path.join(BASE, 'processed', 'entrez_to_symbol.tsv'), sep='\t')
id_col = [c for c in idmap.columns if 'entrez' in c.lower() or 'geneid' in c.lower()][0]
sym_col = [c for c in idmap.columns if 'symbol' in c.lower() or 'gene' in c.lower()][0]
mapping = dict(zip(idmap[id_col].astype(str), idmap[sym_col].astype(str)))
ac4c_symbols = [mapping.get(g, g) for g in targets['gene_symbol'].dropna().astype(str)]
ac4c_set = set(ac4c_symbols)   # 303 unique symbols

# Background universe: full gene set (the 20,162-gene universe used in the original Fisher test).
de = pd.read_csv(os.path.join(BASE, 'processed', 'DE_results_with_symbol.csv'))
bg_genes = list(de['gene_symbol'].astype(str))
N_universe = bg
print(f"ac4C targets: {len(ac4c_set)}; universe N={N_universe}; full universe genes={len(bg_genes)}")

def fisher_or(k, a, b, N):
    table = [[k, a - k], [b - k, N - a - b + k]]
    return stats.fisher_exact(table)[0]

obs_or = fisher_or(len(overlap_genes), ac4c_n, gwas_n, N_universe)
print(f"Observed AD-GWAS enrichment: overlap={len(overlap_genes)} OR={obs_or:.2f}")

# (A) random 88-gene sets from the FULL universe (ac4C status is a gene property within it)
N_ITER = 10000
rand_ors = []
for _ in range(N_ITER):
    samp = set(np.random.choice(bg_genes, gwas_n, replace=False))
    k = len(ac4c_set & samp)
    rand_ors.append(fisher_or(k, ac4c_n, gwas_n, N_universe))
rand_ors = np.array(rand_ors)
p_bg = np.mean(rand_ors >= obs_or)
z_bg = (obs_or - rand_ors.mean()) / rand_ors.std() if rand_ors.std() > 0 else np.inf
print(f"Random background 88-gene sets: mean OR={rand_ors.mean():.2f}, max random OR={rand_ors.max():.2f}, "
      f"observed z={z_bg if np.isfinite(z_bg) else 'inf'}, empirical p={p_bg:.2e}")

# (B) curated brain-disease gene pool (neurodegeneration/psychiatry/neurodevelopmental)
brain_pool = set("""
APP PSEN1 PSEN2 MAPT GRN TBK1 TARDBP SOD1 FUS C9orf72 ATXN1 ATXN2 ATXN3
HTT PRNP SNCA LRRK2 PARK2 PARK7 PINK1 GBA DJ1 VPS35 CHCHD2 DNAJC6
SORL1 BIN1 CLU ABCA7 TREM2 CD33 PICALM CR1 INPP5D MEF2C PTK2B
SPI1 ACE ADAM10 BACE1 BACE2 CD2AP EPHA1 NDUFS2 NDUFAB1 ZNF592
SCARB2 TYROBP PLCG2 FERMT2 HLA-DRB1 HLA-DRB5 CASS4 CELF1 MS4A4A MS4A6A
SLC24A4 RIN3 DSG2 APOE SLC01B1 CTSD PLD3 UNC5C IQCK FERMT2
GRIN2B GRM3 CACNA1C CACNB2 SYNGAP1 SHANK2 SHANK3 CHD8 CNTNAP2
NRXN1 CNTN5 ANK2 DLG2 DPP6 DPP10 KCNQ2 KCNQ3 SCN1A SCN2A SCN8A
HTR2A DRD2 DRD3 COMT BDNF NTRK2 FKBP5 NR3C1 RELN FOXP2
NLGN3 NLGN4X TSC1 TSC2 PTEN NF1 KCNQ5 SCN1B SCN3A GABRA1 GABRB3
SLC6A4 SLC6A3 HTR1A HTR2C DRD4 GAD1 GAD2 CACNA1A CACNA1B CACNG2
MAP2 NEFL GFAP OLIG2 PLP1 MBP NTRK1 NTRK3 NGFR BDNF-AS
PPP3CA PPP3R1 CALB1 CALB2 CAMK2A CAMK2B PRKCA PRKCB DLG4 GRIN1
APBB1 APBB2 APBA1 APBA2 NAE1 RBFOX1 RBFOX3 CELF1 ELAVL3 ELAVL4
STX1A STX1B SNAP25 SYN1 SYN2 SYT1 VAMP2 SYP NRXN2 NLGN1 NLGN2
PSENEN NCSTN APH1A PEN2 BACE1 BACE2 ADAM17 ADAM9 ADAM10
CHAT SLC18A3 SLC5A7 SLC6A11 SLC1A2 SLC1A3 GLUL GLUD1
PRNP PRND MAPT TAU CRYAB HSPB1 HSPB8 VCP SQSTM1 OPTN UBQLN1
FUS TARDBP UBQLN2 C9orf72 MATR3 KIF5A TBK1 OPTN VCP CHMP2B
PSPR SNCA SNCB SNCG GBA LRRK2 PARK2 PARK7 PINK1 HTRA2 DJ1
SLC25A46 MFN2 FBXO7 ATP13A2 PLA2G6 SYNJ1 DNAJC13
GRN FTCD MAPT TBK1 SORL1 TMEM106B C9orf72
SHANK1 SHANK2 SHANK3 DLGAP1 DLGAP2 DLGAP3 DLGAP4 CYFIP1
CNTNAP1 CNTNAP2 NRXN1 NRXN2 NRXN3 NLGN1 NLGN2 NLGN3 NLGN4X
TCF4 CHD8 ARID1B SYNGAP1 DYRK1A ADNP MECP2 FOXP1 FOXP2 GRIN2A
KDM5B KDM6B EHMT1 NSD1 PHF21A AUTS2
""".upper().replace('-','').split())
brain_pool_in_universe = [g for g in brain_pool if g in set(bg_genes)]
print(f"Brain-disease pool size (in universe): {len(brain_pool_in_universe)}")
brain_ors = []
bd_pool_ok = len(brain_pool_in_universe) >= gwas_n
if bd_pool_ok:
    for _ in range(N_ITER):
        samp = set(np.random.choice(brain_pool_in_universe, gwas_n, replace=False))
        k = len(ac4c_set & samp)
        brain_ors.append(fisher_or(k, ac4c_n, gwas_n, N_universe))
brain_ors = np.array(brain_ors)
z_bd = np.nan; p_bd = np.nan
if bd_pool_ok and brain_ors.std() > 0:
    z_bd = (obs_or - brain_ors.mean()) / brain_ors.std()
    p_bd = np.mean(brain_ors >= obs_or)
    print(f"Random brain-disease 88-gene sets: mean OR={brain_ors.mean():.2f}, observed z={z_bd:.1f}, p={p_bd:.2e}")
else:
    print("Brain-disease pool too small / degenerate; reporting pool size only.")

out = pd.DataFrame({
    'control': ['observed_AD_GWAS', 'random_background', 'random_brain_disease'],
    'n_sets': [1, N_ITER, len(brain_ors) if bd_pool_ok else 0],
    'mean_OR': [obs_or, rand_ors.mean(), brain_ors.mean() if bd_pool_ok and len(brain_ors) else np.nan],
    'z_vs_observed': [0, z_bg, z_bd if bd_pool_ok else np.nan],
    'p_vs_observed': [np.nan, p_bg, p_bd if bd_pool_ok else np.nan],
})
out.to_csv(os.path.join(results_dir, 'gwas_specificity_control.csv'), index=False)
print("\nSaved: gwas_specificity_control.csv")
