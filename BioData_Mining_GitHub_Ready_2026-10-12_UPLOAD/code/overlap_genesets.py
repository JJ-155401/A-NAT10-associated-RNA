"""Orthogonal proteomic validation: PRIDE PXD055784 synaptosome proteomics
(5 AD + 5 HC) overlap with ac4C/AD-GWAS/NAT10-KD gene sets.
Uses local UniProt->Gene mapping (HUMAN_9606_idmapping.dat.gz).
"""
import pandas as pd, numpy as np, os, gzip, re
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('AD_NAT10_DATA_DIR', os.path.abspath(os.path.join(HERE, '..', 'data', 'external')))
DATA = os.path.join(ROOT, 'pride_pxd055784')
REF = os.path.join(ROOT, 'reference')
OUT = os.path.join(HERE, '..', 'results')
os.makedirs(OUT, exist_ok=True)

# ---------- UniProt AC -> gene symbol ----------
print('Loading UniProt mapping...')
ac2gene = {}
with gzip.open(os.path.join(REF, 'uniprot_human_idmapping.dat.gz'), 'rt', encoding='utf-8') as f:
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if len(parts) >= 3 and parts[1] == 'Gene_Name':
            ac2gene[parts[0]] = parts[2]
print(f'  {len(ac2gene)} mappings loaded')

# ---------- load gene lists ----------
ac4c_targets = set(pd.read_csv(os.path.join(HERE, '..', 'data', 'human_ac4c_targets_curated.csv'))['gene_symbol'])
ad_risk = set(pd.read_csv(os.path.join(ROOT,'derived','ad_gwas_risk_genes_real.csv'))['gene'])
kd_ad = pd.read_csv(os.path.join(ROOT,'derived','nat10_kd_targets_in_AD.csv'))
kd_ad_genes = set(kd_ad['symbol'])

# ---------- synaptosome proteins -> genes ----------
de = pd.read_csv(os.path.join(OUT, 'synaptosome_protein_DE.csv'))
prot_ids = list(de['Protein'])
gene2prot = {ac2gene.get(p): p for p in prot_ids if ac2gene.get(p)}
detected_genes = set(gene2prot.keys())
print(f'  synaptosome proteins mapped to genes: {len(detected_genes)}')

acetylated = set(l.strip() for l in open(os.path.join(OUT,'acetyl_protein_list.txt')) if l.strip())
acetylated_genes = set(ac2gene.get(p) for p in acetylated if ac2gene.get(p))
print(f'  acetylated proteins mapped to genes: {len(acetylated_genes)}')

# ---------- enrichment (Fisher exact) ----------
def fisher(geneset, background):
    universe = detected_genes
    gs = geneset & universe
    a = len(gs)
    in_set_bg = len(gs & background)
    in_set_notbg = a - in_set_bg
    notset_bg = len(background & universe) - in_set_bg
    notset_notbg = len(universe) - a - notset_bg - in_set_notbg
    _, p = stats.fisher_exact([[in_set_bg, notset_bg],[in_set_notbg, notset_notbg]])
    return a, in_set_bg, p

rows = []
for name, gs in [('ac4C_targets', ac4c_targets), ('AD_GWAS_risk', ad_risk), ('NAT10KD_AD', kd_ad_genes)]:
    a_det, obs_det, p_det = fisher(gs, detected_genes)
    a_ac, obs_ac, p_ac = fisher(gs, acetylated_genes)
    rows.append({'gene_set':name,'set_size':len(gs),
                 'detected_n':obs_det,'detected_total':a_det,'detected_enrich_p':p_det,
                 'acetylated_n':obs_ac,'acetylated_total':a_ac,'acetylated_enrich_p':p_ac})
    print(f'{name}: detected {obs_det}/{a_det} (p={p_det:.2e}); acetylated {obs_ac}/{a_ac} (p={p_ac:.2e})')

enr = pd.DataFrame(rows)
enr.to_csv(os.path.join(OUT,'synaptosome_gene_set_overlap.csv'), index=False)

# ---------- directional consistency for NAT10-KD AD targets ----------
kd_det = kd_ad[kd_ad['symbol'].isin(detected_genes)].copy().rename(columns={'log2FC':'log2FC_KD'})
kd_det['uniprot'] = kd_det['symbol'].map(gene2prot)
kd_det = kd_det.merge(de[['Protein','log2FC','p','q']].rename(columns={'log2FC':'log2FC_AD'}),
                      left_on='uniprot', right_on='Protein', how='left')
# NAT10 low in AD -> target expression follows KD direction.
# KD log2FC<0 (target lower in KD) & protein log2FC_AD<0 (lower in AD) => consistent
kd_det['protein_AD_dir'] = np.where(kd_det['log2FC_AD']>0,'up','down')
kd_det['consistent'] = ((kd_det['log2FC_KD']<0)&(kd_det['log2FC_AD']<0)) | ((kd_det['log2FC_KD']>0)&(kd_det['log2FC_AD']>0))
kd_det.to_csv(os.path.join(OUT,'nat10kd_targets_synaptosome_protein.csv'), index=False)
n_consistent = int(kd_det['consistent'].sum())
print(f'\nNAT10-KD AD targets detected in synaptosome: {len(kd_det)}; directionally consistent: {n_consistent}/{len(kd_det)}')
print(kd_det[['symbol','log2FC_KD','log2FC_AD','protein_AD_dir','consistent']].to_string())

print('\nSaved: synaptosome_gene_set_overlap.csv, nat10kd_targets_synaptosome_protein.csv')
print('DONE')
