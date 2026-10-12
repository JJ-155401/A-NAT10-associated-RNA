"""NAT10 genetic & structural evidence (free public APIs) for NC mechanism pillar.
- AlphaFold structure availability (AF-Q9H0A4)
- ClinVar pathogenic variants in NAT10
- gnomAD constraint (pLI / missense z) for NAT10
- GWAS Catalog traits linked to NAT10 region (already computed: gwas_catalog_nat10_real.csv)
"""
import pandas as pd, os, urllib.request, urllib.parse, json
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('AD_NAT10_DATA_DIR', os.path.abspath(os.path.join(HERE, '..', 'data', 'external')))
OUT = os.path.join(HERE, '..', 'results'); FIG = os.path.join(HERE, '..', 'figures')
os.makedirs(OUT, exist_ok=True)

NAT10_UNIPROT='Q9H0A4'; NAT10_GENE='NAT10'; NAT10_ENSG='ENSG00000135372'

def get_json(url, params=None, headers=None):
    if params: url += '?'+urllib.parse.urlencode(params)
    req=urllib.request.Request(url, headers=headers or {'User-Agent':'Mozilla/5.0'})
    for _ in range(3):
        try:
            r=urllib.request.urlopen(req, timeout=60); return json.loads(r.read())
        except Exception as e:
            print('  retry', str(e)[:50]); 
    return None

rows=[]
# 1. AlphaFold (EBI files API)
af_url=f'https://alphafold.ebi.ac.uk/files/AF-{NAT10_UNIPROT}-F1-model_v4.cif'
try:
    req=urllib.request.Request(af_url, method='HEAD', headers={'User-Agent':'Mozilla/5.0'})
    r=urllib.request.urlopen(req, timeout=30); af_ok = (r.status==200)
except Exception as e:
    af_ok=False
rows.append({'evidence':'AlphaFold structure', 'value': 'available (model_v4)' if af_ok else 'absent',
             'detail': af_url if af_ok else ''})
print('AlphaFold Q9H0A4:', 'available' if af_ok else 'absent')

# 2. ClinVar (via EBI ClinVar REST)
cv=get_json(f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',
            {'db':'clinvar','term':f'{NAT10_GENE}[gene] AND (pathogenic OR likely_pathogenic)','retmode':'json'})
cv_n = cv['esearchresult']['count'] if cv else 'NA'
rows.append({'evidence':'ClinVar pathogenic/likely-pathogenic variants', 'value': cv_n, 'detail':''})
print('ClinVar pathogenic NAT10 variants:', cv_n)

# 3. gnomAD LoF constraint (graphql v4) via POST
import json as _json
gquery='{gene(gene_id:"%s", reference_genome:GRCh38){constraints{lof{pLI}missense{z}synonymous{z}}}}' % NAT10_ENSG
try:
    req=urllib.request.Request('https://gnomad.broadinstitute.org/api/',
        data=_json.dumps({'query':gquery}).encode(),
        headers={'Content-Type':'application/json','User-Agent':'Mozilla/5.0'})
    gr=urllib.request.urlopen(req, timeout=60)
    gj=_json.loads(gr.read())
    con=gj.get('data',{}).get('gene',{}).get('constraints',{})
    if con:
        rows.append({'evidence':'gnomAD pLI','value':con.get('lof',{}).get('pLI'),'detail':'LoF intolerance'})
        rows.append({'evidence':'gnomAD missense z','value':con.get('missense',{}).get('z'),'detail':''})
        rows.append({'evidence':'gnomAD syn z','value':con.get('synonymous',{}).get('z'),'detail':''})
        print('gnomAD pLI:', con.get('lof',{}).get('pLI'))
    else:
        rows.append({'evidence':'gnomAD constraint','value':'NA','detail':'no data'})
except Exception as e:
    rows.append({'evidence':'gnomAD constraint','value':'NA','detail':str(e)[:40]})

# 4. existing GWAS catalog traits
gw=pd.read_csv(os.path.join(ROOT,'derived','gwas_catalog_nat10_real.csv'))
rows.append({'evidence':'GWAS Catalog traits at NAT10 locus', 'value': len(gw),
             'detail': '; '.join(gw['trait'].head(5).astype(str))})

ev=pd.DataFrame(rows)
ev.to_csv(os.path.join(OUT,'nat10_genetic_structural_evidence.csv'), index=False)
print('\nSaved nat10_genetic_structural_evidence.csv')
print(ev.to_string())
print('DONE')
