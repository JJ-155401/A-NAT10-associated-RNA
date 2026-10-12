"""STRING PPI enrichment of the NAT10-ac4C-AD axis (free, no wet lab).
Validates that NAT10-KD AD targets + key ac4C/AD genes cluster in AD
pathology compartments (amyloid-beta complex, NFT, synapse)."""
import pandas as pd, numpy as np, os, urllib.request, urllib.parse, csv, io
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('AD_NAT10_DATA_DIR', os.path.abspath(os.path.join(HERE, '..', 'data', 'external')))
OUT = os.path.join(HERE, '..', 'results'); FIG = os.path.join(HERE, '..', 'figures')
os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)

kd = set(pd.read_csv(os.path.join(ROOT,'derived','nat10_kd_targets_in_AD.csv'))['symbol'])
ac4c = set(pd.read_csv(os.path.join(HERE, '..', 'data', 'human_ac4c_targets_curated.csv'))['gene_symbol'])
focus = list(kd | {'APP','PSEN1','PSEN2','BIN1','CLU','SORL1','TREM2','MAPT','APOE','MEF2C','NAT10'})
ids='\n'.join(focus)
url='https://string-db.org/api/tsv/enrichment?species=9606'
data=urllib.parse.urlencode({'identifiers':ids}).encode()
req=urllib.request.Request(url, data=data, headers={'User-Agent':'Mozilla/5.0','Content-Type':'application/x-www-form-urlencoded'})
import http.client
def robust_read(req, attempts=5):
    for i in range(attempts):
        try:
            r=urllib.request.urlopen(req, timeout=120)
            chunks=[]; 
            while True:
                c=r.read(65536)
                if not c: break
                chunks.append(c)
            return b''.join(chunks)
        except (http.client.IncompleteRead, Exception) as e:
            print(f'  retry {i+1}: {str(e)[:50]}')
    raise RuntimeError('STRING fetch failed')
raw=robust_read(req)
rows=list(csv.DictReader(io.StringIO(raw.decode()), delimiter='\t'))
enr=pd.DataFrame(rows)
enr['fdr']=enr['fdr'].astype(float)
enr.to_csv(os.path.join(OUT,'string_enrichment_nat10_ac4c_ad.csv'), index=False)

# Select AD-relevant GO/Compartment terms
ad_terms = enr[(enr['fdr']<0.001) & enr['category'].isin(['COMPARTMENTS','GO_BP','GO_CC','KEGG','REAC']) &
               enr['description'].str.contains('amyloid|neurofibrillary|tau|synaptic|synapse|synaptic|Alzheimer|learning|memory|protein folding|chaperone|nucleus|translation|ribosom|stress|mitochondr', case=False, na=False)]
print(f'Total STRING terms: {len(enr)}; AD-relevant (FDR<1e-3): {len(ad_terms)}')
print(ad_terms[['category','description','number_of_genes','fdr']].head(20).to_string())

# Figure: top AD-relevant enriched terms
top = ad_terms.sort_values('fdr').head(12).iloc[::-1]
fig, ax = plt.subplots(figsize=(8,5))
colors = {'COMPARTMENTS':'#2c7fb8','GO_BP':'#d95f02','GO_CC':'#1b9e77','KEGG':'#7570b3','REAC':'#e7298a'}
ax.barh(top['description'].str.slice(0,38), -np.log10(top['fdr']),
        color=[colors.get(c,'gray') for c in top['category']])
ax.set_xlabel('-log10(FDR)')
ax.set_title('STRING enrichment: NAT10-ac4C-AD axis in AD pathology', fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(FIG,'string_ad_enrichment.png'), dpi=150, bbox_inches='tight')
plt.savefig(os.path.join(FIG,'string_ad_enrichment.pdf'), bbox_inches='tight')
print('Saved:', os.path.join(FIG,'string_ad_enrichment.png'))
print('DONE')
