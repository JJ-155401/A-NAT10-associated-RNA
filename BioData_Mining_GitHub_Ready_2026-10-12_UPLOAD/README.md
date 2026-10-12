# A NAT10-associated RNA buffering layer links genetic susceptibility to neuronal vulnerability in Alzheimer's disease

Reproducibility package for the manuscript submitted to **BioData Mining**. The
repository contains the analysis scripts, the study-derived curated ac4C target
list, and CSV source data for the manuscript figures, tables, and reported
statistics.

This repository is a reproducibility and source-data package. It is not a copy
of every public archive used by the study: large GEO, single-cell, proteomics,
GWAS, eQTL, and annotation files remain at their original public accessions.
See [`docs/data_availability.md`](docs/data_availability.md) before running the
pipeline.

## Repository layout

```text
code/                 Python analysis scripts
data/                 Study-derived curated ac4C target list and data notes
docs/                 Reproduction, data, and manuscript mapping notes
results/              CSV source data and final statistical summaries
figures/              Generated figures (created when scripts are run)
requirements.txt      Tested Python dependencies
LICENSE               MIT license for this package
```

The bundled file `data/human_ac4c_targets_curated.csv` contains the 303-gene
curated ac4C target set used by the final manuscript. The CSV files in
`results/` are the available source data for the manuscript and should be
uploaded with the code.

## Analyses represented

| Module | Main scripts | Manuscript use |
|---|---|---|
| Bulk transcriptomics | `nat10_gse125583.py`, `ac4c_axis_gse125583.py`, `nat10_gse95587_replication.py` | NAT10/ac4C scores, AD versus control comparisons, Braak associations, and independent replication |
| Single-nucleus RNA-seq | `single_cell_cross_validation.py`, `check_counts.py` | Cell-type-resolved analysis in GSE147528 and validation in GSE138852 |
| Proteomics and RNA-protein concordance | `download_synapse.py`, `analyze_pride_synaptosome.py`, `rna_protein_concordance.py`, `string_enrichment.py` | PXD055784 synaptosome analysis, protein overlap, concordance, and STRING enrichment |
| Genetic convergence | `coloc_analysis.py`, `gwas_*.py`, `neuronal_vulnerability_score*.py`, `permutation_test*.py` | AD GWAS enrichment, matched null models, colocalization, and neuronal vulnerability analyses |
| Sensitivity and causal analyses | `causal_mediation_sensitivity.py`, `recalculate_all_stats.py`, `recalculate_or*.py`, `sensitivity_ac4c_targets.py` | Mediation sensitivity, odds-ratio checks, target-set sensitivity, and consolidated statistic checks |

## Quick start

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -r requirements.txt
python code/nat10_gse125583.py
```

Scripts write outputs to `results/` and generated figures to `figures/`.
Large public inputs are expected under `data/external/`. To keep those inputs
elsewhere, set the environment variable `AD_NAT10_DATA_DIR` to that directory
before running a script. For the complete sequence and expected input layout,
see [`docs/reproduction.md`](docs/reproduction.md).

The scripts use fixed random seeds where randomization is part of the analysis
(normally seed 42; the specificity-control script uses its documented seed).

## Public data and resources

| Resource | Accession or identifier | Role |
|---|---|---|
| AD bulk RNA-seq discovery cohort | GSE125583 | 219 AD and 70 control samples, n = 289 |
| AD bulk RNA-seq replication cohort | GSE95587 | 84 AD and 33 control samples, n = 117 |
| NAT10 perturbation and ac4C-related data | GSE102113 | NAT10 knockdown, ac4C-seq/acRIP and RNA-seq evidence |
| Human AD single-nucleus RNA-seq | GSE147528 | Zhou et al. cohort; 69,865 high-quality cells in the manuscript |
| Independent single-nucleus validation | GSE138852 | Mathys et al. cohort; 13,214 total cells in the manuscript and 11,503 cells after the analysis QC filter |
| Synaptosome proteomics | PXD055784 | Independent protein-level validation |
| AD GWAS | Bellenguez et al. 2022 / GCST90027158 | AD risk-gene convergence |
| Brain eQTL | GTEx v10 | NAT10-region colocalization input |
| Protein interaction and pathway annotation | STRING v11 | Functional enrichment |
| Protein and tissue annotation | Human Protein Atlas | NAT10 tissue and protein context |
| Additional public proteomics | ROSMAP open proteomics | External RNA-protein comparison where used |
| 5xFAD ac4C resource | Publicly archived study data | Disease-model ac4C comparison where used |

The exact download notes and expected local file names are listed in
[`docs/data_availability.md`](docs/data_availability.md). No passwords, tokens,
or private human-level data are included in this repository.

## Interpreting the source results

The final manuscript reports several analyses with different target/background
definitions. In particular, the main GWAS enrichment result and the target-list
sensitivity result must not be merged: the repository preserves both outputs.
Likewise, 13,214 is the reported total Mathys-cohort cell count, whereas
11,503 is the number remaining after the analysis QC and alignment steps used
for the cross-cohort analysis. See
[`docs/manuscript_mapping.md`](docs/manuscript_mapping.md).

## License

MIT. See [`LICENSE`](LICENSE).
