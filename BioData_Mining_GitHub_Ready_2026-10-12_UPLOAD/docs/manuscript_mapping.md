# Manuscript to repository mapping

## Manuscript title

**A NAT10-associated RNA buffering layer links genetic susceptibility to
neuronal vulnerability in Alzheimer's disease**

## Main analysis mapping

| Manuscript content | Repository location |
|---|---|
| Bulk NAT10/ac4C analysis and Braak associations | `code/nat10_gse125583.py`, `code/ac4c_axis_gse125583.py`, and matching CSVs in `results/` |
| Independent bulk replication | `code/nat10_gse95587_replication.py`, `results/external_replication_GSE95587.csv` |
| Zhou single-nucleus analysis | `results/single_cell_analysis.csv` |
| Mathys independent validation | `code/single_cell_cross_validation.py`, `results/single_cell_mathys2019.csv` |
| Cross-cohort summary | `results/cross_validation_summary.csv` |
| AD GWAS enrichment and sensitivity | `results/ad_gwas_enrichment.csv`, `gwas_enrichment_expression_matched.csv`, `gwas_enrichment_sensitivity.csv`, and leave-one-out/null outputs |
| Proteomic validation and RNA-protein concordance | `results/synaptosome_protein_DE.csv`, `results/rna_protein_concordance.csv`, and summary files |
| Mediation and robustness analyses | `results/causal_mediation_imai_tingley.csv`, `medsens_*.csv`, and sensitivity outputs |

## Denominator notes

The main GWAS enrichment file preserves the analysis with 303 curated targets,
88 AD GWAS genes, 13 overlapping genes, and the reported odds ratio of about
11.8. `gwas_enrichment_sensitivity.csv` contains a separate target-list and
background sensitivity analysis in which the current curated list is evaluated
with 10 overlaps and an odds ratio of about 17.72. These outputs use different
analysis definitions and must not be merged into one claim.

For GSE138852, 13,214 is the total cell count described for the validation
cohort, while 11,503 is the post-QC, aligned count used by the cross-validation
script. Both values may be correct when referring to their respective stages.
