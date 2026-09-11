# RBM3-HCC-MPC transcriptomics

## Study
**Context-dependent association of RBM3 with a macropinocytosis-associated transcriptional state in hepatocellular carcinoma**

This public repository accompanies the Scientific Reports submission and exposes the post-lock pipeline-sensitivity analysis, machine-readable sensitivity results, Supplementary Table S4, and v2 file-integrity provenance.

## Main interpretation
The transcriptomic MPC scores are literature-curated macropinocytosis-associated transcriptional modules, not direct measurements of macropinocytic flux. The analyses support a cohort-dependent RBM3-MPC transcriptional association in HCC and do not establish a direct RBM3-to-macropinocytosis regulatory mechanism.

## Reproducibility v2
The v2 package adds the post-lock Xena/GDC pipeline-sensitivity analysis. Xena and GDC contained the same 371 TCGA-LIHC sample keys. After exact matching and rescoring, RBM3-MPC Spearman correlations remained 0.361 and 0.636, respectively. Cross-pipeline concordance was high for RBM3 (rho = 0.886) and more moderate for the composite MPC score (rho = 0.655), indicating preprocessing sensitivity of the precise TCGA effect magnitude rather than a sample-composition artifact.

Publicly browsable files include:
- `00_code/postlock_sensitivity/28_TCGA_Xena_GDC_pipeline_sensitivity.R`
- `reference_tables/TCGA_Xena_GDC_pipeline_sensitivity_summary.csv`
- `reference_tables/TCGA_Xena_GDC_pipeline_sensitivity_genes.csv`
- `reference_tables/TCGA_Xena_GDC_matched_samples.csv`
- `reference_tables/Supplementary_Table_S4_TCGA_Xena_GDC_pipeline_sensitivity.xlsx`
- `provenance/reproducibility_manifest_v2.csv`
- `provenance/reproducibility_MD5_check_v2.csv`
- `RBM3_HCC_reproducibility_v2_SUBMISSION_CHECKSUMS.txt`

A complete submission reproducibility archive was generated and retained separately for journal submission. Its integrity record is included in this repository.

## Archive integrity
- MD5: `472e95e7fa5eedbd0c2354177c71877c`
- SHA-256: `3e2ea23aedcafaaee70df664e18d1b580393c5511cbe33a61cac6b94e44d100f`

The v2 payload contains 71 manifest-tracked files; all 71 passed the final MD5 verification. The original v1 local lock was retained separately and was not overwritten.

## Raw data
Raw public datasets are not redistributed. Dataset sources and accession identifiers are reported in the manuscript Methods section.

## R environment
The original analysis workstation exported a full `sessionInfo()` snapshot during v1 locking. That local file was not available to the repository-packaging environment, so package-version details were not reconstructed or invented.

## Repository
https://github.com/b2859186481-creator/RBM3-HCC-MPC-transcriptomics

## License
No explicit software license has been assigned at this stage.