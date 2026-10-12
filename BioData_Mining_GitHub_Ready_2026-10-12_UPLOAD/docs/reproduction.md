# Reproduction guide

## 1. Environment

The scripts were tested with Python 3.10.11. Install the pinned dependencies
from the repository root:

```bash
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell
# .venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

## 2. Input directory convention

The repository contains the curated 303-gene target list at
`data/human_ac4c_targets_curated.csv`. Large public inputs are not bundled.
Place their downloaded and locally processed files under `data/external/`, or
set `AD_NAT10_DATA_DIR` to an equivalent directory. The scripts use paths
relative to the repository and write reproducible outputs to `results/`.

The expected external layout is:

```text
data/external/
  reference/
    uniprot_human_idmapping.dat.gz
  pride_pxd055784/
    peptides.txt
    modificationSpecificPeptides.txt
  processed/
    GSE125583_metadata.csv
    GSE125583_counts.parquet
    GSE95587_metadata.csv
    GSE95587_raw_counts.parquet
    DE_results_with_symbol.csv
    entrez_to_symbol.tsv
    gene_lengths.tsv                 # optional for length-matched nulls
  scrna/
    zhou/zhou_merged.h5ad
    GSE138852_counts.csv.gz
    GSE138852_covariates.csv.gz
  downloads/
    gwas/ad_gwas_nat10_region.csv
    ac4c_targets/arango2018_ac4c_targets.csv
    ac4c_targets/GSE102113_acRIP_peaks.csv
  gtex/gtex_v10_eqtl_real.csv
  derived/
    ad_gwas_risk_genes_real.csv
    gwas_catalog_nat10_real.csv
    nat10_kd_targets_in_AD.csv
    iv_details.csv
    real_ko_ad_merged_ac4c.csv
```

The precise file names above describe the interfaces expected by the scripts;
they do not imply that these files are included in the GitHub repository.

## 3. Public data retrieval

Download the source data from the accessions in
[`data_availability.md`](data_availability.md), then create the processed
tables expected above using the same gene identifiers and sample labels used
in the manuscript. `download_synapse.py` can assist with Synapse-hosted files,
but PRIDE/PXD055784 and other large resources may require their own access
procedure. Do not commit raw `.h5ad`, `.gz`, `.tar.gz`, `.zip`, or `.7z` files.

## 4. Recommended run order

1. Run the bulk discovery and replication scripts:
   `nat10_gse125583.py`, `ac4c_axis_gse125583.py`, and
   `nat10_gse95587_replication.py`.
2. Run the single-nucleus analysis and cross-cohort validation:
   `single_cell_cross_validation.py`.
3. Run the proteomics and RNA-protein analyses:
   `analyze_pride_synaptosome.py`, `overlap_genesets.py`, and
   `rna_protein_concordance.py`.
4. Run the GWAS, colocalization, neuronal-vulnerability, and permutation
   scripts, including `coloc_analysis.py`, `gwas_*.py`,
   `neuronal_vulnerability_score*.py`, and `permutation_test*.py`.
5. Run sensitivity and mediation analyses, then use
   `recalculate_all_stats.py` for the consolidated checks.

Some scripts are alternative or reviewer-sensitivity analyses and overwrite a
CSV with the same name. Preserve the committed `results/` files before
rerunning a sensitivity workflow if you need to compare versions.

## 5. Source data and deterministic checks

The committed CSV files in `results/` are the source data accompanying the
manuscript figures and reported summary statistics. They allow review of the
reported values without downloading every large public input. Randomized
permutation procedures use fixed seeds; exact reproduction still depends on
the public-data release, preprocessing, gene annotation release, and the
external input files matching the versions used in the study.

## 6. Important cell-count distinction

For GSE138852, the manuscript describes 13,214 total cells in the validation
cohort. The cross-validation output records 11,503 cells after the analysis QC
and cell-alignment steps. These are different denominators and should not be
substituted for one another.
