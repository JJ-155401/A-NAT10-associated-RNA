# Data directory

`human_ac4c_targets_curated.csv` is the study-derived 303-gene ac4C target
list used in the final manuscript and is included in the repository.

Large public inputs are intentionally not committed. Place them under
`data/external/` using the layout in [`../docs/reproduction.md`](../docs/reproduction.md),
or set `AD_NAT10_DATA_DIR` to an external directory. The ignored external
directory may contain downloaded GEO, single-cell, PRIDE, GWAS, GTEx, UniProt,
and annotation files needed for full reruns.
