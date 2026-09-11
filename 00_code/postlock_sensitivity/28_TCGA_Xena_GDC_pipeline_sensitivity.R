# ============================================================
# 28_TCGA_Xena_GDC_pipeline_sensitivity.R
#
# Goal:
# Explain the difference between the TCGA-LIHC RBM3-MPC
# association in:
#   UCSC Xena PanCanAtlas: rho ~ 0.361
#   GDC STAR-Counts/TPM:   rho ~ 0.636
#
# This is a post-lock sensitivity analysis.
# It does NOT overwrite the previously locked primary analyses.
# ============================================================

mpc_genes <- c(
  "RAC1", "CDC42", "RHOG", "ARF6", "ARF1", "RAB5A", "PAK1",
  "CTBP1", "PLD1", "DOCK1", "SNX5", "SLC9A1", "SDC1", "SLC4A7"
)

mean_z_score <- function(expr_mat, genes) {
  missing_genes <- setdiff(genes, rownames(expr_mat))
  if (length(missing_genes) > 0) {
    stop("Missing genes: ", paste(missing_genes, collapse = ", "))
  }
  x <- expr_mat[genes, , drop = FALSE]
  z <- t(scale(t(x)))
  colMeans(z, na.rm = TRUE)
}

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) {
    return(c(rho = NA_real_, p = NA_real_, n = sum(ok)))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

xena_scored <- readRDS("02_processed_data/PanCanAtlas_MPC_scored_primary_tumors.rds")
gdc_expr <- readRDS("02_processed_data/TCGA_LIHC_expr_log.rds")

xena_lihc <- xena_scored[xena_scored$cancer == "LIHC", , drop = FALSE]
cat("\nXena LIHC samples:", nrow(xena_lihc), "\n")
cat("GDC LIHC samples:", ncol(gdc_expr), "\n")

xena_lihc$sample_key <- substr(xena_lihc$sample, 1, 15)
gdc_sample_key <- substr(colnames(gdc_expr), 1, 15)

xena_dup <- unique(xena_lihc$sample_key[duplicated(xena_lihc$sample_key)])
gdc_dup <- unique(gdc_sample_key[duplicated(gdc_sample_key)])
cat("Duplicated Xena sample keys:", length(xena_dup), "\n")
cat("Duplicated GDC sample keys:", length(gdc_dup), "\n")
if (length(xena_dup) > 0) stop("Duplicate Xena sample keys detected.")
if (length(gdc_dup) > 0) stop("Duplicate GDC sample keys detected.")

common_keys <- intersect(xena_lihc$sample_key, gdc_sample_key)
xena_only <- setdiff(xena_lihc$sample_key, gdc_sample_key)
gdc_only <- setdiff(gdc_sample_key, xena_lihc$sample_key)

cat("\n============================================\n")
cat("SAMPLE MATCHING\n")
cat("============================================\n")
cat("Xena n:", nrow(xena_lihc), "\n")
cat("GDC n:", ncol(gdc_expr), "\n")
cat("Matched n:", length(common_keys), "\n")
cat("Xena-only n:", length(xena_only), "\n")
cat("GDC-only n:", length(gdc_only), "\n")
if (length(xena_only) > 0) print(xena_only)
if (length(gdc_only) > 0) print(gdc_only)

xena_original <- safe_spearman(xena_lihc$RBM3, xena_lihc$MPC_formation_score)
gdc_score_all <- mean_z_score(gdc_expr, mpc_genes)
gdc_original <- safe_spearman(as.numeric(gdc_expr["RBM3", , drop = TRUE]), gdc_score_all)

cat("\n============================================\n")
cat("ORIGINAL PIPELINE RESULTS\n")
cat("============================================\n")
cat("Xena original rho:", xena_original["rho"], "\n")
cat("Xena original P:", xena_original["p"], "\n")
cat("Xena original n:", xena_original["n"], "\n\n")
cat("GDC original rho:", gdc_original["rho"], "\n")
cat("GDC original P:", gdc_original["p"], "\n")
cat("GDC original n:", gdc_original["n"], "\n")

xena_idx <- match(common_keys, xena_lihc$sample_key)
gdc_idx <- match(common_keys, gdc_sample_key)
xena_gene_columns <- c("RBM3", mpc_genes)

missing_xena <- setdiff(xena_gene_columns, colnames(xena_lihc))
missing_gdc <- setdiff(xena_gene_columns, rownames(gdc_expr))
if (length(missing_xena) > 0) stop("Missing Xena genes: ", paste(missing_xena, collapse = ", "))
if (length(missing_gdc) > 0) stop("Missing GDC genes: ", paste(missing_gdc, collapse = ", "))

xena_match_expr <- t(as.matrix(xena_lihc[xena_idx, xena_gene_columns, drop = FALSE]))
colnames(xena_match_expr) <- common_keys
gdc_match_expr <- gdc_expr[xena_gene_columns, gdc_idx, drop = FALSE]
colnames(gdc_match_expr) <- common_keys
stopifnot(identical(colnames(xena_match_expr), colnames(gdc_match_expr)))

xena_match_score <- mean_z_score(xena_match_expr, mpc_genes)
gdc_match_score <- mean_z_score(gdc_match_expr, mpc_genes)
xena_matched <- safe_spearman(as.numeric(xena_match_expr["RBM3", , drop = TRUE]), xena_match_score)
gdc_matched <- safe_spearman(as.numeric(gdc_match_expr["RBM3", , drop = TRUE]), gdc_match_score)

cat("\n============================================\n")
cat("IDENTICAL-SAMPLE, RE-SCORED COMPARISON\n")
cat("============================================\n")
cat("Xena matched rho:", xena_matched["rho"], "\n")
cat("Xena matched P:", xena_matched["p"], "\n")
cat("GDC matched rho:", gdc_matched["rho"], "\n")
cat("GDC matched P:", gdc_matched["p"], "\n")

rbm3_cross <- safe_spearman(
  as.numeric(xena_match_expr["RBM3", , drop = TRUE]),
  as.numeric(gdc_match_expr["RBM3", , drop = TRUE])
)
mpc_cross <- safe_spearman(xena_match_score, gdc_match_score)

cat("\n============================================\n")
cat("CROSS-PIPELINE CONCORDANCE\n")
cat("============================================\n")
cat("RBM3 Xena-vs-GDC rho:", rbm3_cross["rho"], "\n")
cat("RBM3 concordance P:", rbm3_cross["p"], "\n\n")
cat("MPC score Xena-vs-GDC rho:", mpc_cross["rho"], "\n")
cat("MPC score concordance P:", mpc_cross["p"], "\n")

gene_diag <- do.call(
  rbind,
  lapply(mpc_genes, function(g) {
    cross <- safe_spearman(
      as.numeric(xena_match_expr[g, , drop = TRUE]),
      as.numeric(gdc_match_expr[g, , drop = TRUE])
    )
    xena_rbm3 <- safe_spearman(
      as.numeric(xena_match_expr["RBM3", , drop = TRUE]),
      as.numeric(xena_match_expr[g, , drop = TRUE])
    )
    gdc_rbm3 <- safe_spearman(
      as.numeric(gdc_match_expr["RBM3", , drop = TRUE]),
      as.numeric(gdc_match_expr[g, , drop = TRUE])
    )
    data.frame(
      Gene = g,
      Cross_pipeline_rho = as.numeric(cross["rho"]),
      Cross_pipeline_P = as.numeric(cross["p"]),
      Xena_RBM3_gene_rho = as.numeric(xena_rbm3["rho"]),
      GDC_RBM3_gene_rho = as.numeric(gdc_rbm3["rho"]),
      Delta_GDC_minus_Xena = as.numeric(gdc_rbm3["rho"]) - as.numeric(xena_rbm3["rho"]),
      stringsAsFactors = FALSE
    )
  })
)
gene_diag <- gene_diag[order(gene_diag$Cross_pipeline_rho), , drop = FALSE]

cat("\n============================================\n")
cat("GENE-LEVEL PIPELINE DIAGNOSTICS\n")
cat("============================================\n")
print(gene_diag, row.names = FALSE)

summary_results <- data.frame(
  Metric = c(
    "Xena_LIHC_n", "GDC_LIHC_n", "Matched_sample_n", "Xena_only_n", "GDC_only_n",
    "Xena_original_RBM3_MPC_rho", "GDC_original_RBM3_MPC_rho",
    "Xena_matched_rescored_RBM3_MPC_rho", "GDC_matched_rescored_RBM3_MPC_rho",
    "Cross_pipeline_RBM3_rho", "Cross_pipeline_MPC_score_rho"
  ),
  Value = c(
    nrow(xena_lihc), ncol(gdc_expr), length(common_keys), length(xena_only), length(gdc_only),
    as.numeric(xena_original["rho"]), as.numeric(gdc_original["rho"]),
    as.numeric(xena_matched["rho"]), as.numeric(gdc_matched["rho"]),
    as.numeric(rbm3_cross["rho"]), as.numeric(mpc_cross["rho"])
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_results, "05_tables/TCGA_Xena_GDC_pipeline_sensitivity_summary.csv", row.names = FALSE)
write.csv(gene_diag, "05_tables/TCGA_Xena_GDC_pipeline_sensitivity_genes.csv", row.names = FALSE)
write.csv(data.frame(sample_key = common_keys, stringsAsFactors = FALSE), "05_tables/TCGA_Xena_GDC_matched_samples.csv", row.names = FALSE)

cat("\n============================================\n")
cat("PIPELINE SENSITIVITY ANALYSIS COMPLETE\n")
cat("============================================\n")
print(summary_results, row.names = FALSE)
cat(
  "\nSaved:\n",
  "05_tables/TCGA_Xena_GDC_pipeline_sensitivity_summary.csv\n",
  "05_tables/TCGA_Xena_GDC_pipeline_sensitivity_genes.csv\n",
  "05_tables/TCGA_Xena_GDC_matched_samples.csv\n"
)
