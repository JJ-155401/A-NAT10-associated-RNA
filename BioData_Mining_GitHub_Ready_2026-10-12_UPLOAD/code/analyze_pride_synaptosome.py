"""Analysis of PRIDE PXD055784: human brain synaptosome proteomics (5 AD vs 5 HC).
Orthogonal proteomic validation of AD-ac4C hypothesis without wet lab.
- Protein-level differential abundance (peptides.txt)
- Acetyl-PTM landscape (modificationSpecificPeptides.txt)
- NAT10 (Q9H0A4) detection
"""
import pandas as pd, numpy as np, os, re
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get('AD_NAT10_DATA_DIR', os.path.join(HERE, '..', 'data', 'external'))
DATA = os.path.join(BASE, 'pride_pxd055784')
OUT = os.path.join(HERE, '..', 'results')
os.makedirs(OUT, exist_ok=True)

pep = pd.read_csv(os.path.join(DATA, 'peptides.txt'), sep='\t')
intensity_cols = [c for c in pep.columns if c.startswith('Intensity ')]
samples = [c.replace('Intensity ', '') for c in intensity_cols]
sample_group = {s: ('AD' if s.startswith('AD') else 'HC') for s in samples}
ad_cols = [c for c in intensity_cols if sample_group[c.replace('Intensity ', '')] == 'AD']
hc_cols = [c for c in intensity_cols if sample_group[c.replace('Intensity ', '')] == 'HC']

def uniprots(s):
    if not isinstance(s, str): return []
    return re.findall(r'([A-Z0-9]{6,10})', s.split(';')[0])

# ---- NAT10 peptide intensity ----
print('=== NAT10 (Q9H0A4 / NAT10) peptide-level ===')
nat10 = pep[pep['Proteins'].astype(str).str.contains('Q9H0A4')]
if len(nat10):
    for _, row in nat10.iterrows():
        ad = row[ad_cols].astype(float); hc = row[hc_cols].astype(float)
        print(f"  seq={row['Sequence']} AD_med={ad.median():.0f} HC_med={hc.median():.0f} log2FC={np.log2((ad.median()+1)/(hc.median()+1)):.3f}")
else:
    print('  NAT10 peptides not detected (expected: low-abundance writer enzyme; see MSBB/ROSMAP note)')

# ---- Protein-level differential abundance ----
pep['_prot'] = pep['Proteins'].apply(lambda s: uniprots(s)[0] if uniprots(s) else None)
pv = pep.dropna(subset=['_prot'])
prot_int = pv.groupby('_prot')[intensity_cols].sum(min_count=1)
prot_log = np.log2(prot_int.replace(0, np.nan))
ad_present = (prot_int[ad_cols] > 0).sum(axis=1)
hc_present = (prot_int[hc_cols] > 0).sum(axis=1)
keep = (ad_present >= 5) | (hc_present >= 5)
prot_log = prot_log[keep]

results = []
for prot, row in prot_log.iterrows():
    ad = row[[c for c in row.index if c in ad_cols]].dropna()
    hc = row[[c for c in row.index if c in hc_cols]].dropna()
    if len(ad) < 3 or len(hc) < 3: continue
    t, p = stats.ttest_ind(ad, hc, equal_var=False)
    lfc = ad.median() - hc.median()
    results.append((prot, lfc, p, len(ad), len(hc), ad.median(), hc.median()))

res = pd.DataFrame(results, columns=['Protein','log2FC','p','nAD','nHC','AD_med','HC_med'])
res['q'] = stats.false_discovery_control(res['p'], method='bh')
res = res.sort_values('p')
res.to_csv(os.path.join(OUT, 'synaptosome_protein_DE.csv'), index=False)
sig = res[res['q'] < 0.05]
print(f'\n=== Differential protein abundance: {prot_log.shape[0]} proteins tested ===')
print(f'  significant (q<0.05): {len(sig)}  (up in AD: {(sig.log2FC>0).sum()}, down: {(sig.log2FC<0).sum()})')
if 'Q9H0A4' in set(res.Protein):
    n = res[res.Protein=='Q9H0A4'].iloc[0]
    print(f'  NAT10 Q9H0A4: log2FC={n.log2FC:.3f} p={n.p:.3e} q={n.q:.3e}')

# ---- Acetyl-PTM landscape ----
print('\n=== Acetyl (Protein N-term) PTM landscape ===')
mod = pd.read_csv(os.path.join(DATA, 'modificationSpecificPeptides.txt'), sep='\t')
ac_col = 'Acetyl (Protein N-term)'
ac = mod[mod[ac_col].notna() & (mod[ac_col].astype(str) != '')].copy()
acetylated = set()
for s in ac['Proteins']: acetylated.update(uniprots(s))
print(f'  acetyl peptides: {len(ac)}, unique acetylated proteins: {len(acetylated)}')
with open(os.path.join(OUT, 'acetyl_protein_list.txt'), 'w') as f:
    for p in sorted(acetylated): f.write(p+'\n')

print('\nSaved: synaptosome_protein_DE.csv, acetyl_protein_list.txt')
print('DONE')
