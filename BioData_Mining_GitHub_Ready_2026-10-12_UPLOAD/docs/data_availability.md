# Data availability and upload scope

## Public resources

The study is based on public resources. The main accession identifiers are
GSE125583, GSE95587, GSE102113, GSE147528, GSE138852, and PXD055784, together
with the Bellenguez et al. 2022 AD GWAS, GTEx v10 eQTL data, STRING v11,
Human Protein Atlas, ROSMAP open proteomics, and the public 5xFAD ac4C
resource. The manuscript and README describe the role of each resource.

The repository does not mirror the original archives. In particular, raw
single-cell matrices, `.h5ad` files, bulk count matrices, PRIDE evidence files,
GWAS summary statistics, GTEx tables, UniProt mapping files, and large
annotation tables are intentionally excluded because of size, provenance, or
redistribution constraints.

## What is included

The uploadable package includes:

- analysis code in `code/`;
- the study-derived `data/human_ac4c_targets_curated.csv` target list;
- CSV source data and statistical summaries in `results/`;
- installation and reproduction documentation.

The result CSVs are not raw data. They are the tabular source data generated
for the manuscript figures, tables, and statistical statements.

## What should be downloaded separately

Download the public inputs from GEO, PRIDE/Synapse, GWAS Catalog/IGAP, GTEx,
UniProt, HPA, STRING, and the relevant public proteomics archive. Store them
outside Git history or under the ignored `data/external/` directory. Use
`AD_NAT10_DATA_DIR` to point the scripts to that location.

## Human data and privacy

The repository contains no direct identifiers, passwords, access tokens, or
private clinical records. Do not add controlled-access individual-level data
to the GitHub repository. If a public resource changes its redistribution
terms, link to the source accession and retain only derived aggregate outputs.
