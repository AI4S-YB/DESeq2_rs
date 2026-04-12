#!/usr/bin/env Rscript
# Export DESeq2 intermediate results on airway dataset for Rust validation.
# Usage: Rscript scripts/export_r_reference.R
# Requires: DESeq2, airway packages installed.

suppressPackageStartupMessages({
  library(DESeq2)
  library(airway)
})

# Load airway data
data("airway")
dds <- DESeqDataSet(airway, design = ~ dex)

# Set reference level (untreated)
dds$dex <- relevel(dds$dex, ref = "untrt")

# Output directory
outdir <- "tests/reference_data"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 1. Export raw counts (genes x samples)
counts_mat <- counts(dds)
write.table(counts_mat, file.path(outdir, "airway_counts.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

# 2. Export coldata
coldata <- as.data.frame(colData(dds)[, "dex", drop = FALSE])
write.table(coldata, file.path(outdir, "airway_coldata.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

# 3. Run DESeq step by step
dds <- estimateSizeFactors(dds)

# Export size factors
sf <- sizeFactors(dds)
write.table(data.frame(sample = names(sf), size_factor = sf),
            file.path(outdir, "r_size_factors.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Export base means (after size factor normalization)
base_means <- rowMeans(counts(dds, normalized = TRUE))
write.table(data.frame(gene = names(base_means), baseMean = base_means),
            file.path(outdir, "r_base_means.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 4. Estimate dispersions
dds <- estimateDispersions(dds)

# Export gene-wise dispersion estimates
write.table(
  data.frame(gene = rownames(dds), dispGeneEst = mcols(dds)$dispGeneEst),
  file.path(outdir, "r_disp_gene_estimates.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Export trend coefficients
disp_fn <- dispersionFunction(dds)
coeffs <- attr(disp_fn, "coefficients")
write.table(
  data.frame(param = names(coeffs), value = coeffs),
  file.path(outdir, "r_disp_trend_coeffs.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Export MAP dispersion estimates
write.table(
  data.frame(gene = rownames(dds), dispMAP = mcols(dds)$dispMAP),
  file.path(outdir, "r_disp_map_estimates.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Export final dispersions used
write.table(
  data.frame(gene = rownames(dds), dispersion = mcols(dds)$dispersion),
  file.path(outdir, "r_dispersions_final.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 5. Wald test
dds <- nbinomWaldTest(dds)

# Export beta coefficients (log2 fold changes for all terms)
beta_mat <- coef(dds)
write.table(cbind(gene = rownames(beta_mat), as.data.frame(beta_mat)),
            file.path(outdir, "r_beta_coefficients.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Export standard errors
se_mat <- mcols(dds)[, grep("^SE_", colnames(mcols(dds))), drop = FALSE]
write.table(cbind(gene = rownames(dds), as.data.frame(se_mat)),
            file.path(outdir, "r_beta_se.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 6. Results with independent filtering
res <- results(dds, contrast = c("dex", "trt", "untrt"), alpha = 0.1)
write.table(
  cbind(gene = rownames(res), as.data.frame(res)),
  file.path(outdir, "r_results.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Export design matrix for reference
mm <- model.matrix(~ dex, data = colData(dds))
write.table(mm, file.path(outdir, "r_design_matrix.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

cat("Exported all reference data to:", outdir, "\n")
cat("Genes:", nrow(dds), "Samples:", ncol(dds), "\n")
cat("Size factors:", paste(round(sf, 6), collapse = ", "), "\n")
cat("Trend coefficients:", paste(names(coeffs), round(coeffs, 6), sep = "=", collapse = ", "), "\n")

# Summary of significant genes at different thresholds
sig_01_1 <- sum(res$padj < 0.01 & abs(res$log2FoldChange) > 1, na.rm = TRUE)
cat("Genes with padj<0.01 & |log2FC|>1:", sig_01_1, "\n")
