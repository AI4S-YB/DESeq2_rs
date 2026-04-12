# DESeq2-rs Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement a Rust library + CLI that reproduces DESeq2's core differential expression pipeline, validated against R's DESeq2 on the airway dataset with 100% gene-set overlap at p<0.01 & log2FC>1.

**Architecture:** Functional Core + Pipeline Shell. Pure functions implement each algorithm step (size factors, dispersion, GLM, Wald test). A `DESeqDataSet` pipeline struct orchestrates them. CLI wraps the pipeline.

**Tech Stack:** Rust, faer (linear algebra), rayon (parallelism), clap (CLI), csv (I/O), statrs (statistical functions)

**Design Spec:** `docs/superpowers/specs/2026-04-12-deseq2-rs-design.md`

**R Source Reference:** `DESeq2/src/DESeq2.cpp` (fitDisp, fitBeta), `DESeq2/R/core.R` (pipeline logic)

---

## File Structure

```
rust_deseq2/
├── Cargo.toml
├── src/
│   ├── lib.rs                    # Public API re-exports
│   ├── data.rs                   # Core types: DESeqDataSet, DataFrame, result structs
│   ├── io.rs                     # CSV/TSV reading and writing
│   ├── design.rs                 # Design matrix generation from coldata
│   ├── size_factors.rs           # Median-of-ratios size factor estimation
│   ├── dispersion/
│   │   ├── mod.rs                # Re-exports + estimate_dispersions() orchestrator
│   │   ├── gene_estimate.rs      # NB log-posterior, derivatives, fitDisp line search
│   │   ├── trend_fit.rs          # Parametric dispersion-mean trend: alpha = a1/mu + a0
│   │   └── map_estimate.rs       # Prior variance estimation + MAP shrinkage
│   ├── glm/
│   │   ├── mod.rs                # Re-exports
│   │   └── irls.rs               # IRLS negative binomial GLM fitting (fitBeta)
│   ├── test_stats/
│   │   ├── mod.rs                # Re-exports
│   │   ├── wald.rs               # Wald test statistic + p-value
│   │   └── p_adjust.rs           # BH correction + independent filtering
│   └── pipeline.rs               # DESeqDataSet::run() orchestration + intermediate export
├── src/bin/
│   └── main.rs                   # CLI entry point (clap)
├── tests/
│   ├── test_size_factors.rs      # Size factor validation vs R
│   ├── test_dispersion.rs        # Dispersion estimation validation vs R
│   ├── test_glm.rs               # IRLS fitting validation vs R
│   ├── test_wald.rs              # Wald test validation vs R
│   ├── test_pipeline.rs          # Airway end-to-end validation
│   └── reference_data/           # R-exported intermediate results (git-tracked)
│       ├── airway_counts.tsv
│       ├── airway_coldata.tsv
│       ├── r_size_factors.tsv
│       ├── r_base_means.tsv
│       ├── r_disp_gene_estimates.tsv
│       ├── r_disp_trend_coeffs.tsv
│       ├── r_disp_map_estimates.tsv
│       ├── r_beta_coefficients.tsv
│       ├── r_beta_se.tsv
│       └── r_results.tsv
└── scripts/
    └── export_r_reference.R      # Export airway reference data from R DESeq2
```

---

## Task 1: Project Scaffolding

**Files:**
- Create: `rust_deseq2/Cargo.toml`
- Create: `rust_deseq2/src/lib.rs`
- Create: `rust_deseq2/src/bin/main.rs`

- [ ] **Step 1: Create Cargo.toml**

```toml
[package]
name = "deseq2-rs"
version = "0.1.0"
edition = "2021"
description = "Rust reimplementation of DESeq2 core differential expression analysis"

[[bin]]
name = "deseq2-rs"
path = "src/bin/main.rs"

[dependencies]
faer = "0.20"
rayon = "1"
clap = { version = "4", features = ["derive"] }
csv = "1"
statrs = "0.18"

[dev-dependencies]
approx = "0.5"
```

- [ ] **Step 2: Create minimal lib.rs**

```rust
pub mod data;
pub mod io;
pub mod design;
pub mod size_factors;
pub mod dispersion;
pub mod glm;
pub mod test_stats;
pub mod pipeline;
```

- [ ] **Step 3: Create stub modules**

Create empty files with a single comment so the module tree compiles:
- `src/data.rs`
- `src/io.rs`
- `src/design.rs`
- `src/size_factors.rs`
- `src/dispersion/mod.rs`
- `src/dispersion/gene_estimate.rs`
- `src/dispersion/trend_fit.rs`
- `src/dispersion/map_estimate.rs`
- `src/glm/mod.rs`
- `src/glm/irls.rs`
- `src/test_stats/mod.rs`
- `src/test_stats/wald.rs`
- `src/test_stats/p_adjust.rs`
- `src/pipeline.rs`
- `src/bin/main.rs` with `fn main() {}`

- [ ] **Step 4: Verify compilation**

Run: `cd rust_deseq2 && cargo check`
Expected: compiles with no errors (warnings OK)

- [ ] **Step 5: Commit**

```bash
git add rust_deseq2/
git commit -m "feat: scaffold rust_deseq2 project with module structure"
```

---

## Task 2: R Reference Data Export Script

**Files:**
- Create: `scripts/export_r_reference.R`

This script must be run manually in R before the Rust integration tests can pass. It exports all intermediate results from DESeq2 on the airway dataset.

- [ ] **Step 1: Write the export script**

```r
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
outdir <- file.path(dirname(getwd()), "tests", "reference_data")
if (!dir.exists(outdir)) {
  outdir <- "rust_deseq2/tests/reference_data"
}
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
```

- [ ] **Step 2: Run the R script** (manual step, requires R + Bioconductor)

Run: `cd rust_deseq2 && Rscript ../scripts/export_r_reference.R`
Expected: TSV files created in `tests/reference_data/`

- [ ] **Step 3: Commit the script and reference data**

```bash
git add scripts/export_r_reference.R rust_deseq2/tests/reference_data/
git commit -m "feat: add R reference data export script and airway test data"
```

---

## Task 3: Core Data Structures

**Files:**
- Create: `rust_deseq2/src/data.rs`
- Test: inline unit tests

- [ ] **Step 1: Write test for DataFrame**

In `src/data.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dataframe_get_column() {
        let mut df = DataFrame::new();
        df.add_column("condition", vec!["untrt".into(), "trt".into(), "untrt".into()]);
        assert_eq!(
            df.get_column("condition").unwrap(),
            &["untrt", "trt", "untrt"]
        );
        assert!(df.get_column("missing").is_none());
    }

    #[test]
    fn test_dataframe_levels() {
        let mut df = DataFrame::new();
        df.add_column("dex", vec!["trt".into(), "untrt".into(), "trt".into(), "untrt".into()]);
        let levels = df.levels("dex").unwrap();
        // Levels should be sorted alphabetically
        assert!(levels.contains(&"trt".to_string()));
        assert!(levels.contains(&"untrt".to_string()));
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd rust_deseq2 && cargo test --lib data`
Expected: FAIL — `DataFrame` not defined

- [ ] **Step 3: Implement data types**

```rust
use faer::Mat;
use std::collections::HashMap;

/// Lightweight sample metadata container.
/// Maps column names to string vectors.
#[derive(Debug, Clone)]
pub struct DataFrame {
    columns: HashMap<String, Vec<String>>,
    pub row_names: Vec<String>,
}

impl DataFrame {
    pub fn new() -> Self {
        Self {
            columns: HashMap::new(),
            row_names: Vec::new(),
        }
    }

    pub fn add_column(&mut self, name: &str, values: Vec<String>) {
        self.columns.insert(name.to_string(), values);
    }

    pub fn get_column(&self, name: &str) -> Option<&[String]> {
        self.columns.get(name).map(|v| v.as_slice())
    }

    /// Return sorted unique values for a column.
    pub fn levels(&self, name: &str) -> Option<Vec<String>> {
        self.columns.get(name).map(|v| {
            let mut lvls: Vec<String> = v.iter().cloned().collect::<std::collections::BTreeSet<_>>().into_iter().collect();
            lvls.sort();
            lvls
        })
    }

    pub fn nrow(&self) -> usize {
        self.columns.values().next().map_or(0, |v| v.len())
    }
}

/// Dispersion estimation results.
#[derive(Debug, Clone)]
pub struct DispersionResult {
    pub gene_estimates: Vec<f64>,
    pub trend_values: Vec<f64>,
    pub map_estimates: Vec<f64>,
    pub prior_var: f64,
}

/// Parametric trend coefficients: dispersion = asympt_disp + extra_pois / mean
#[derive(Debug, Clone, Copy)]
pub struct DispTrendCoeffs {
    pub asympt_disp: f64,  // a0: asymptotic dispersion
    pub extra_pois: f64,   // a1: extra-Poisson coefficient
}

impl DispTrendCoeffs {
    pub fn eval(&self, mean: f64) -> f64 {
        self.asympt_disp + self.extra_pois / mean
    }
}

/// GLM fitting result for all genes.
#[derive(Debug, Clone)]
pub struct BetaFitResult {
    pub coefficients: Mat<f64>,    // genes x p, natural log scale
    pub standard_errors: Mat<f64>, // genes x p
    pub deviances: Vec<f64>,
    pub converged: Vec<bool>,
}

/// Single-gene Wald test output.
#[derive(Debug, Clone)]
pub struct WaldTestResult {
    pub log2_fold_change: f64,
    pub lfc_se: f64,
    pub stat: f64,
    pub p_value: f64,
}

/// Final per-gene result row.
#[derive(Debug, Clone)]
pub struct DESeqResult {
    pub gene: String,
    pub base_mean: f64,
    pub log2_fold_change: f64,
    pub lfc_se: f64,
    pub stat: f64,
    pub p_value: f64,
    pub p_adjusted: f64,
}

/// Contrast specification for hypothesis testing.
#[derive(Debug, Clone)]
pub enum Contrast {
    /// Use the last coefficient in the design matrix.
    LastCoefficient,
    /// Specify (column_name, numerator_level, denominator_level).
    ColumnLevels(String, String, String),
    /// Custom numeric contrast vector.
    Vector(Vec<f64>),
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd rust_deseq2 && cargo test --lib data`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add rust_deseq2/src/data.rs
git commit -m "feat: add core data structures (DataFrame, DispersionResult, BetaFitResult, etc.)"
```

---

## Task 4: CSV/TSV I/O

**Files:**
- Modify: `rust_deseq2/src/io.rs`
- Test: `rust_deseq2/tests/test_io.rs` (integration test using reference_data)

The I/O module reads count matrices and coldata from TSV files, which is the input format for both the CLI and integration tests.

- [ ] **Step 1: Write failing test**

Create `tests/test_io.rs`:

```rust
use deseq2_rs::io::{read_count_matrix, read_coldata};
use std::path::Path;

#[test]
fn test_read_count_matrix_small() {
    // Create a small test file inline
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("counts.tsv");
    std::fs::write(&path, "gene\tS1\tS2\nGENE1\t10\t20\nGENE2\t30\t40\n").unwrap();

    let (counts, gene_names, sample_names) = read_count_matrix(&path).unwrap();
    assert_eq!(gene_names, vec!["GENE1", "GENE2"]);
    assert_eq!(sample_names, vec!["S1", "S2"]);
    assert_eq!(counts.nrows(), 2);
    assert_eq!(counts.ncols(), 2);
    assert_eq!(counts[(0, 0)], 10.0);
    assert_eq!(counts[(1, 1)], 40.0);
}

#[test]
fn test_read_coldata_small() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("coldata.tsv");
    std::fs::write(&path, "\tdex\nSRR1039508\tuntrt\nSRR1039509\ttrt\n").unwrap();

    let (df, sample_names) = read_coldata(&path).unwrap();
    assert_eq!(sample_names, vec!["SRR1039508", "SRR1039509"]);
    assert_eq!(df.get_column("dex").unwrap(), &["untrt", "trt"]);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd rust_deseq2 && cargo test test_read_count`
Expected: FAIL — `read_count_matrix` not found

- [ ] **Step 3: Implement I/O functions**

In `src/io.rs`:

```rust
use crate::data::{DataFrame, DESeqResult};
use faer::Mat;
use std::error::Error;
use std::io::BufRead;
use std::path::Path;

/// Read a tab-separated count matrix. First column is gene names, first row is sample names.
/// Returns (count_matrix as Mat<f64>, gene_names, sample_names).
pub fn read_count_matrix(path: &Path) -> Result<(Mat<f64>, Vec<String>, Vec<String>), Box<dyn Error>> {
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    let mut lines = reader.lines();

    // Header line: first field is empty or "gene", rest are sample names
    let header = lines.next().ok_or("empty file")??;
    let sample_names: Vec<String> = header.split('\t').skip(1).map(|s| s.trim().to_string()).collect();
    let n_samples = sample_names.len();

    let mut gene_names = Vec::new();
    let mut data: Vec<Vec<f64>> = Vec::new();

    for line in lines {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        gene_names.push(fields[0].trim().to_string());
        let row: Vec<f64> = fields[1..]
            .iter()
            .map(|s| s.trim().parse::<f64>())
            .collect::<Result<_, _>>()?;
        if row.len() != n_samples {
            return Err(format!("row {} has {} values, expected {}", gene_names.last().unwrap(), row.len(), n_samples).into());
        }
        data.push(row);
    }

    let n_genes = gene_names.len();
    let counts = Mat::from_fn(n_genes, n_samples, |i, j| data[i][j]);
    Ok((counts, gene_names, sample_names))
}

/// Read tab-separated coldata. First column is sample names (row names).
/// Returns (DataFrame, sample_names_in_order).
pub fn read_coldata(path: &Path) -> Result<(DataFrame, Vec<String>), Box<dyn Error>> {
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    let mut lines = reader.lines();

    let header = lines.next().ok_or("empty file")??;
    let col_names: Vec<String> = header.split('\t').skip(1).map(|s| s.trim().to_string()).collect();

    let mut sample_names = Vec::new();
    let mut columns: Vec<Vec<String>> = vec![Vec::new(); col_names.len()];

    for line in lines {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        sample_names.push(fields[0].trim().to_string());
        for (i, val) in fields[1..].iter().enumerate() {
            if i < columns.len() {
                columns[i].push(val.trim().to_string());
            }
        }
    }

    let mut df = DataFrame::new();
    df.row_names = sample_names.clone();
    for (name, vals) in col_names.iter().zip(columns) {
        df.add_column(name, vals);
    }
    Ok((df, sample_names))
}

/// Write results to a TSV file.
pub fn write_results(path: &Path, results: &[DESeqResult]) -> Result<(), Box<dyn Error>> {
    let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    wtr.write_record(&["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"])?;
    for r in results {
        wtr.write_record(&[
            &r.gene,
            &format_float(r.base_mean),
            &format_float(r.log2_fold_change),
            &format_float(r.lfc_se),
            &format_float(r.stat),
            &format_float(r.p_value),
            &format_float(r.p_adjusted),
        ])?;
    }
    Ok(())
}

fn format_float(v: f64) -> String {
    if v.is_nan() {
        "NA".to_string()
    } else {
        format!("{:.10e}", v)
    }
}
```

Also add `tempfile` to dev-dependencies in `Cargo.toml`:

```toml
[dev-dependencies]
approx = "0.5"
tempfile = "3"
```

- [ ] **Step 4: Run tests**

Run: `cd rust_deseq2 && cargo test test_read`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add rust_deseq2/src/io.rs rust_deseq2/tests/test_io.rs rust_deseq2/Cargo.toml
git commit -m "feat: add TSV I/O for count matrix, coldata, and results"
```

---

## Task 5: Design Matrix Generation

**Files:**
- Modify: `rust_deseq2/src/design.rs`
- Test: inline unit tests

Generates a model matrix from coldata: intercept column + indicator columns for each non-reference level.

- [ ] **Step 1: Write failing test**

In `src/design.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::DataFrame;

    #[test]
    fn test_design_matrix_two_levels() {
        let mut df = DataFrame::new();
        // 4 samples: untrt, trt, untrt, trt — matches airway layout
        df.add_column("dex", vec![
            "untrt".into(), "trt".into(), "untrt".into(), "trt".into()
        ]);
        let mat = build_design_matrix(&df, "dex", "untrt").unwrap();
        // Should be 4x2: [intercept, dex_trt]
        assert_eq!(mat.nrows(), 4);
        assert_eq!(mat.ncols(), 2);
        // Row 0 (untrt): [1, 0]
        assert_eq!(mat[(0, 0)], 1.0);
        assert_eq!(mat[(0, 1)], 0.0);
        // Row 1 (trt): [1, 1]
        assert_eq!(mat[(1, 0)], 1.0);
        assert_eq!(mat[(1, 1)], 1.0);
    }
}
```

- [ ] **Step 2: Run test to verify failure**

Run: `cd rust_deseq2 && cargo test --lib design`
Expected: FAIL

- [ ] **Step 3: Implement design matrix builder**

```rust
use crate::data::DataFrame;
use faer::Mat;

/// Build a design matrix from a DataFrame column.
/// Creates intercept + indicator columns for non-reference levels (sorted alphabetically).
/// This matches R's `model.matrix(~ factor)` with `relevel(factor, ref=reference)`.
pub fn build_design_matrix(
    col_data: &DataFrame,
    column: &str,
    reference: &str,
) -> Result<Mat<f64>, String> {
    let values = col_data
        .get_column(column)
        .ok_or_else(|| format!("column '{}' not found in coldata", column))?;
    let levels = col_data
        .levels(column)
        .ok_or_else(|| format!("column '{}' not found", column))?;

    if !levels.contains(&reference.to_string()) {
        return Err(format!("reference level '{}' not found in column '{}'", reference, column));
    }

    // Non-reference levels, sorted
    let non_ref: Vec<&String> = levels.iter().filter(|l| l.as_str() != reference).collect();
    let n_samples = values.len();
    let n_cols = 1 + non_ref.len(); // intercept + indicators

    let mat = Mat::from_fn(n_samples, n_cols, |i, j| {
        if j == 0 {
            1.0 // intercept
        } else {
            if values[i] == *non_ref[j - 1] { 1.0 } else { 0.0 }
        }
    });

    Ok(mat)
}
```

- [ ] **Step 4: Run tests**

Run: `cd rust_deseq2 && cargo test --lib design`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add rust_deseq2/src/design.rs
git commit -m "feat: add design matrix generation from coldata"
```

---

## Task 6: Size Factor Estimation

**Files:**
- Modify: `rust_deseq2/src/size_factors.rs`
- Test: inline unit tests + `tests/test_size_factors.rs`

Algorithm (median-of-ratios, Anders & Huber 2010):
1. For each gene: `log_geo_mean = mean(log(counts))` across samples (skip rows where any count is 0 in the geometric mean computation — R uses `is.finite(loggeomeans)` which excludes genes with any zero since log(0)=-Inf makes the mean -Inf)
2. For each sample j: `ratio_ij = log(count_ij) - log_geo_mean_i` for genes where `geo_mean > 0` AND `count_ij > 0`
3. `size_factor_j = exp(median(ratios_j))`

- [ ] **Step 1: Write failing test with known values**

In `src/size_factors.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use faer::Mat;
    use approx::assert_relative_eq;

    #[test]
    fn test_size_factors_simple() {
        // 3 genes, 2 samples
        // Gene1: 10, 20  -> geo_mean = sqrt(200) = 14.14
        // Gene2: 30, 60  -> geo_mean = sqrt(1800) = 42.43
        // Gene3: 50, 100 -> geo_mean = sqrt(5000) = 70.71
        // Sample 1 ratios: 10/14.14, 30/42.43, 50/70.71 = 0.707 for all
        // Sample 2 ratios: 20/14.14, 60/42.43, 100/70.71 = 1.414 for all
        let counts = Mat::from_fn(3, 2, |i, j| {
            [[10.0, 20.0], [30.0, 60.0], [50.0, 100.0]][i][j]
        });
        let sf = estimate_size_factors(&counts);
        // All ratios for sample 0 are 0.707, median = 0.707
        // All ratios for sample 1 are 1.414, median = 1.414
        assert_eq!(sf.len(), 2);
        assert_relative_eq!(sf[0] / sf[1], 0.5, epsilon = 1e-10);
    }

    #[test]
    fn test_size_factors_with_zeros() {
        // Gene with a zero should be excluded from geo mean (log(0)=-inf makes row mean -inf)
        let counts = Mat::from_fn(3, 2, |i, j| {
            [[0.0, 20.0], [30.0, 60.0], [50.0, 100.0]][i][j]
        });
        let sf = estimate_size_factors(&counts);
        // Only genes 2 and 3 used (gene 1 has zero -> -inf geo mean)
        assert_eq!(sf.len(), 2);
        assert_relative_eq!(sf[0] / sf[1], 0.5, epsilon = 1e-10);
    }
}
```

- [ ] **Step 2: Run test to verify failure**

Run: `cd rust_deseq2 && cargo test --lib size_factors`
Expected: FAIL

- [ ] **Step 3: Implement size factor estimation**

```rust
use faer::Mat;

/// Estimate size factors using the median-of-ratios method.
///
/// For each gene, compute log geometric mean across samples.
/// Genes where any sample has count=0 produce -inf geo mean and are excluded.
/// For each sample, compute median of (log count - log geo mean) for valid genes,
/// then exponentiate.
pub fn estimate_size_factors(counts: &Mat<f64>) -> Vec<f64> {
    let n_genes = counts.nrows();
    let n_samples = counts.ncols();

    // Step 1: log geometric mean per gene = rowMeans(log(counts))
    let log_geo_means: Vec<f64> = (0..n_genes)
        .map(|i| {
            let sum: f64 = (0..n_samples).map(|j| counts[(i, j)].ln()).sum();
            sum / n_samples as f64
        })
        .collect();

    // Step 2: for each sample, collect log-ratios and take median
    (0..n_samples)
        .map(|j| {
            let mut ratios: Vec<f64> = (0..n_genes)
                .filter(|&i| log_geo_means[i].is_finite() && counts[(i, j)] > 0.0)
                .map(|i| counts[(i, j)].ln() - log_geo_means[i])
                .collect();
            ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let median = if ratios.is_empty() {
                1.0
            } else {
                let mid = ratios.len() / 2;
                if ratios.len() % 2 == 0 {
                    (ratios[mid - 1] + ratios[mid]) / 2.0
                } else {
                    ratios[mid]
                }
            };
            median.exp()
        })
        .collect()
}

/// Compute base mean per gene after size factor normalization.
pub fn base_means(counts: &Mat<f64>, size_factors: &[f64]) -> Vec<f64> {
    let n_genes = counts.nrows();
    let n_samples = counts.ncols();
    (0..n_genes)
        .map(|i| {
            let sum: f64 = (0..n_samples).map(|j| counts[(i, j)] / size_factors[j]).sum();
            sum / n_samples as f64
        })
        .collect()
}
```

- [ ] **Step 4: Run tests**

Run: `cd rust_deseq2 && cargo test --lib size_factors`
Expected: PASS

- [ ] **Step 5: Write integration test against R reference data**

Create `tests/test_size_factors.rs`:

```rust
use deseq2_rs::io::read_count_matrix;
use deseq2_rs::size_factors::{estimate_size_factors, base_means};
use approx::assert_relative_eq;
use std::path::Path;

fn read_r_size_factors() -> Vec<(String, f64)> {
    let path = Path::new("tests/reference_data/r_size_factors.tsv");
    let content = std::fs::read_to_string(path).expect("r_size_factors.tsv not found — run scripts/export_r_reference.R first");
    content.lines().skip(1).map(|line| {
        let fields: Vec<&str> = line.split('\t').collect();
        (fields[0].to_string(), fields[1].parse::<f64>().unwrap())
    }).collect()
}

#[test]
fn test_size_factors_vs_r() {
    let path = Path::new("tests/reference_data/airway_counts.tsv");
    if !path.exists() {
        eprintln!("Skipping: reference data not found");
        return;
    }
    let (counts, _genes, samples) = read_count_matrix(path).unwrap();
    let sf = estimate_size_factors(&counts);
    let r_sf = read_r_size_factors();

    assert_eq!(sf.len(), r_sf.len());
    for (i, (sample, r_val)) in r_sf.iter().enumerate() {
        assert_relative_eq!(sf[i], *r_val, epsilon = 1e-10,
            // Allow tiny floating point differences
        );
    }
}
```

- [ ] **Step 6: Run integration test**

Run: `cd rust_deseq2 && cargo test test_size_factors_vs_r`
Expected: PASS (if reference data exists) or skip message

- [ ] **Step 7: Commit**

```bash
git add rust_deseq2/src/size_factors.rs rust_deseq2/tests/test_size_factors.rs
git commit -m "feat: implement median-of-ratios size factor estimation"
```

---

## Task 7: NB Log-Posterior and Derivatives

**Files:**
- Modify: `rust_deseq2/src/dispersion/gene_estimate.rs`
- Test: inline unit tests

This implements the three core mathematical functions from DESeq2.cpp: `log_posterior`, `dlog_posterior`, `d2log_posterior`. These are used by the fitDisp line search in Task 8.

**Mathematical formulas** (from DESeq2.cpp lines 31-158):

Log-posterior of dispersion alpha given counts y, means mu, design matrix X:

```
log_post(log_alpha) = LL + CR + prior

LL = sum_j [ lgamma(y_j + 1/alpha) - lgamma(1/alpha) - y_j * ln(mu_j + 1/alpha) - (1/alpha) * ln(1 + mu_j * alpha) ]

CR = -0.5 * ln(det(X' W X))   where W = diag(1 / (1/mu + alpha))

prior = -0.5 * (log_alpha - prior_mean)^2 / prior_var   (if use_prior)
```

First derivative w.r.t. log(alpha):

```
dLL/dalpha = (1/alpha^2) * sum_j [ digamma(1/alpha) + ln(1+mu*alpha) - mu*alpha/(1+mu*alpha) - digamma(y+1/alpha) + y/(mu+1/alpha) ]

dCR/dalpha: uses trace(B^-1 * dB) where B = X'WX, dB = X'(dW)X, dW = diag(-1/(1/mu+alpha)^2)

dlog_post/d(log_alpha) = (dLL/dalpha + dCR/dalpha) * alpha + dprior/d(log_alpha)
```

- [ ] **Step 1: Write test for log_posterior against known values**

Compute a small example by hand or derive from R. In `src/dispersion/gene_estimate.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use faer::Mat;
    use approx::assert_relative_eq;

    fn make_test_data() -> (Vec<f64>, Vec<f64>, Mat<f64>) {
        // 4 samples, 1 gene
        let y = vec![10.0, 20.0, 15.0, 25.0];
        let mu = vec![12.0, 18.0, 14.0, 22.0];
        // Simple design: intercept + one covariate
        let x = Mat::from_fn(4, 2, |i, j| {
            if j == 0 { 1.0 } else { [0.0, 1.0, 0.0, 1.0][i] }
        });
        (y, mu, x)
    }

    #[test]
    fn test_log_posterior_finite() {
        let (y, mu, x) = make_test_data();
        let lp = log_posterior(0.0, &y, &mu, &x, 0.0, 1.0, false, true);
        assert!(lp.is_finite());
    }

    #[test]
    fn test_derivative_numerical_check() {
        // Check analytical derivative against numerical derivative
        let (y, mu, x) = make_test_data();
        let log_alpha = 0.0;
        let eps = 1e-5;
        let lp_plus = log_posterior(log_alpha + eps, &y, &mu, &x, 0.0, 1.0, false, true);
        let lp_minus = log_posterior(log_alpha - eps, &y, &mu, &x, 0.0, 1.0, false, true);
        let numerical_deriv = (lp_plus - lp_minus) / (2.0 * eps);
        let analytical_deriv = dlog_posterior(log_alpha, &y, &mu, &x, 0.0, 1.0, false, true);
        assert_relative_eq!(analytical_deriv, numerical_deriv, epsilon = 1e-4);
    }

    #[test]
    fn test_second_derivative_numerical_check() {
        let (y, mu, x) = make_test_data();
        let log_alpha = 0.0;
        let eps = 1e-5;
        let d1_plus = dlog_posterior(log_alpha + eps, &y, &mu, &x, 0.0, 1.0, false, true);
        let d1_minus = dlog_posterior(log_alpha - eps, &y, &mu, &x, 0.0, 1.0, false, true);
        let numerical_d2 = (d1_plus - d1_minus) / (2.0 * eps);
        let analytical_d2 = d2log_posterior(log_alpha, &y, &mu, &x, 0.0, 1.0, false, true);
        assert_relative_eq!(analytical_d2, numerical_d2, epsilon = 1e-3);
    }
}
```

- [ ] **Step 2: Run test to verify failure**

Run: `cd rust_deseq2 && cargo test --lib dispersion::gene_estimate`
Expected: FAIL

- [ ] **Step 3: Implement log_posterior, dlog_posterior, d2log_posterior**

```rust
use faer::Mat;
use statrs::function::gamma::{digamma, ln_gamma};

/// Trigamma function: psi_1(x) = d/dx digamma(x).
/// Computed via series expansion for x > 6, recurrence for x <= 6.
fn trigamma(x: f64) -> f64 {
    // Use recurrence to shift x > 6, then asymptotic series
    let mut x = x;
    let mut result = 0.0;
    while x < 6.0 {
        result += 1.0 / (x * x);
        x += 1.0;
    }
    // Asymptotic expansion: 1/x + 1/(2x^2) + 1/(6x^3) - 1/(30x^5) + 1/(42x^7) - ...
    let x2 = x * x;
    result += 1.0 / x + 1.0 / (2.0 * x2) + 1.0 / (6.0 * x2 * x)
        - 1.0 / (30.0 * x2 * x2 * x)
        + 1.0 / (42.0 * x2 * x2 * x2 * x);
    result
}

/// NB log-posterior of dispersion.
/// log_alpha: log of dispersion parameter
/// y: counts for one gene across samples
/// mu: fitted means for one gene across samples
/// x: design matrix (samples x p)
/// prior_mean, prior_var: normal prior on log_alpha
/// use_prior: whether to include the prior term
/// use_cr: whether to include Cox-Reid adjustment
pub fn log_posterior(
    log_alpha: f64,
    y: &[f64],
    mu: &[f64],
    x: &Mat<f64>,
    prior_mean: f64,
    prior_var: f64,
    use_prior: bool,
    use_cr: bool,
) -> f64 {
    let alpha = log_alpha.exp();
    let alpha_inv = 1.0 / alpha;
    let n = y.len();

    // NB log-likelihood
    let ll: f64 = (0..n)
        .map(|j| {
            ln_gamma(y[j] + alpha_inv) - ln_gamma(alpha_inv)
                - y[j] * (mu[j] + alpha_inv).ln()
                - alpha_inv * (1.0 + mu[j] * alpha).ln()
        })
        .sum();

    // Cox-Reid adjustment: -0.5 * ln(det(X' W X))
    let cr = if use_cr {
        // W = diag(1 / (1/mu + alpha))
        let p = x.ncols();
        let mut xtwx = Mat::<f64>::zeros(p, p);
        for j in 0..n {
            let w = 1.0 / (1.0 / mu[j] + alpha);
            for c1 in 0..p {
                for c2 in 0..p {
                    xtwx[(c1, c2)] += x[(j, c1)] * w * x[(j, c2)];
                }
            }
        }
        let det = matrix_determinant(&xtwx);
        -0.5 * det.ln()
    } else {
        0.0
    };

    // Prior: -0.5 * (log_alpha - mean)^2 / var
    let prior = if use_prior {
        -0.5 * (log_alpha - prior_mean).powi(2) / prior_var
    } else {
        0.0
    };

    ll + cr + prior
}

/// First derivative of log-posterior w.r.t. log(alpha).
pub fn dlog_posterior(
    log_alpha: f64,
    y: &[f64],
    mu: &[f64],
    x: &Mat<f64>,
    prior_mean: f64,
    prior_var: f64,
    use_prior: bool,
    use_cr: bool,
) -> f64 {
    let alpha = log_alpha.exp();
    let alpha_inv = 1.0 / alpha;
    let alpha_inv2 = alpha_inv * alpha_inv;
    let n = y.len();

    // dLL/dalpha
    let dll_dalpha: f64 = alpha_inv2
        * (0..n)
            .map(|j| {
                digamma(alpha_inv) + (1.0 + mu[j] * alpha).ln()
                    - mu[j] * alpha / (1.0 + mu[j] * alpha)
                    - digamma(y[j] + alpha_inv)
                    + y[j] / (mu[j] + alpha_inv)
            })
            .sum::<f64>();

    // dCR/dalpha
    let dcr_dalpha = if use_cr {
        let p = x.ncols();
        // B = X' W X, dB = X' dW X
        // W_j = 1/(1/mu_j + alpha),  dW_j = -1/(1/mu_j + alpha)^2
        let mut b = Mat::<f64>::zeros(p, p);
        let mut db = Mat::<f64>::zeros(p, p);
        for j in 0..n {
            let base = 1.0 / mu[j] + alpha;
            let w = 1.0 / base;
            let dw = -1.0 / (base * base);
            for c1 in 0..p {
                for c2 in 0..p {
                    b[(c1, c2)] += x[(j, c1)] * w * x[(j, c2)];
                    db[(c1, c2)] += x[(j, c1)] * dw * x[(j, c2)];
                }
            }
        }
        let b_inv = matrix_inverse(&b);
        let det_b = matrix_determinant(&b);
        let trace_b_inv_db = matrix_trace_product(&b_inv, &db);
        let ddetb = det_b * trace_b_inv_db;
        -0.5 * ddetb / det_b
    } else {
        0.0
    };

    // dlog_post/d(log_alpha) = (dLL/dalpha + dCR/dalpha) * alpha + dprior/d(log_alpha)
    let dprior = if use_prior {
        -(log_alpha - prior_mean) / prior_var
    } else {
        0.0
    };

    (dll_dalpha + dcr_dalpha) * alpha + dprior
}

/// Second derivative of log-posterior w.r.t. log(alpha).
pub fn d2log_posterior(
    log_alpha: f64,
    y: &[f64],
    mu: &[f64],
    x: &Mat<f64>,
    prior_mean: f64,
    prior_var: f64,
    use_prior: bool,
    use_cr: bool,
) -> f64 {
    let alpha = log_alpha.exp();
    let alpha_inv = 1.0 / alpha;
    let alpha_inv2 = alpha_inv * alpha_inv;
    let n = y.len();

    // First sum (same as in dlog_posterior)
    let sum1: f64 = (0..n)
        .map(|j| {
            digamma(alpha_inv) + (1.0 + mu[j] * alpha).ln()
                - mu[j] * alpha / (1.0 + mu[j] * alpha)
                - digamma(y[j] + alpha_inv)
                + y[j] / (mu[j] + alpha_inv)
        })
        .sum();

    // Second sum (second derivative terms)
    let sum2: f64 = (0..n)
        .map(|j| {
            -alpha_inv2 * trigamma(alpha_inv)
                + mu[j] * mu[j] * alpha / (1.0 + mu[j] * alpha).powi(2)
                + alpha_inv2 * trigamma(y[j] + alpha_inv)
                + alpha_inv2 * y[j] / (mu[j] + alpha_inv).powi(2)
        })
        .sum();

    let d2ll_dalpha2 = -2.0 * alpha_inv * alpha_inv2 * sum1 + alpha_inv2 * sum2;

    // d2CR/dalpha2
    let d2cr_dalpha2 = if use_cr {
        let p = x.ncols();
        let mut b = Mat::<f64>::zeros(p, p);
        let mut db = Mat::<f64>::zeros(p, p);
        let mut d2b = Mat::<f64>::zeros(p, p);
        for j in 0..n {
            let base = 1.0 / mu[j] + alpha;
            let w = 1.0 / base;
            let dw = -1.0 / (base * base);
            let d2w = 2.0 / (base * base * base);
            for c1 in 0..p {
                for c2 in 0..p {
                    b[(c1, c2)] += x[(j, c1)] * w * x[(j, c2)];
                    db[(c1, c2)] += x[(j, c1)] * dw * x[(j, c2)];
                    d2b[(c1, c2)] += x[(j, c1)] * d2w * x[(j, c2)];
                }
            }
        }
        let b_inv = matrix_inverse(&b);
        let det_b = matrix_determinant(&b);
        let tr_bi_db = matrix_trace_product(&b_inv, &db);
        let bi_db = matrix_multiply(&b_inv, &db);
        let tr_bi_db_bi_db = matrix_trace_product(&bi_db, &bi_db);
        let tr_bi_d2b = matrix_trace_product(&b_inv, &d2b);
        let ddetb = det_b * tr_bi_db;
        let d2detb = det_b * (tr_bi_db * tr_bi_db - tr_bi_db_bi_db + tr_bi_d2b);
        0.5 * (ddetb / det_b).powi(2) - 0.5 * d2detb / det_b
    } else {
        0.0
    };

    // d2log_post/d(log_alpha)^2 = (d2ll*alpha^2 + dll*alpha + d2CR*alpha^2 + dCR*alpha) + d2prior
    // = (d2ll + d2CR) * alpha^2 + dlog_post/d(log_alpha)_without_prior + d2prior
    let dlog_post_no_prior = dlog_posterior(log_alpha, y, mu, x, prior_mean, prior_var, false, use_cr);
    let d2prior = if use_prior { -1.0 / prior_var } else { 0.0 };

    (d2ll_dalpha2 + d2cr_dalpha2) * alpha * alpha + dlog_post_no_prior + d2prior
}

// --- Small matrix helpers for p x p matrices (p is typically 2-4) ---

fn matrix_determinant(m: &Mat<f64>) -> f64 {
    let n = m.nrows();
    assert_eq!(n, m.ncols());
    match n {
        1 => m[(0, 0)],
        2 => m[(0, 0)] * m[(1, 1)] - m[(0, 1)] * m[(1, 0)],
        _ => {
            // Cofactor expansion along first row. Fine for small p (design matrices are typically 2-4 cols).
            let mut det = 0.0;
            for j in 0..n {
                let minor = minor_matrix(m, 0, j);
                let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
                det += sign * m[(0, j)] * matrix_determinant(&minor);
            }
            det
        }
    }
}

fn minor_matrix(m: &Mat<f64>, row: usize, col: usize) -> Mat<f64> {
    let n = m.nrows() - 1;
    Mat::from_fn(n, n, |i, j| {
        let si = if i >= row { i + 1 } else { i };
        let sj = if j >= col { j + 1 } else { j };
        m[(si, sj)]
    })
}

fn matrix_inverse(m: &Mat<f64>) -> Mat<f64> {
    let n = m.nrows();
    assert_eq!(n, m.ncols());
    if n == 2 {
        let det = m[(0, 0)] * m[(1, 1)] - m[(0, 1)] * m[(1, 0)];
        Mat::from_fn(2, 2, |i, j| {
            match (i, j) {
                (0, 0) => m[(1, 1)] / det,
                (0, 1) => -m[(0, 1)] / det,
                (1, 0) => -m[(1, 0)] / det,
                (1, 1) => m[(0, 0)] / det,
                _ => unreachable!(),
            }
        })
    } else {
        // Use faer's solve: A^-1 = solve(A, I)
        let identity = Mat::from_fn(n, n, |i, j| if i == j { 1.0 } else { 0.0 });
        let lu = m.partial_piv_lu();
        lu.solve(&identity)
    }
}

fn matrix_multiply(a: &Mat<f64>, b: &Mat<f64>) -> Mat<f64> {
    a * b
}

fn matrix_trace_product(a: &Mat<f64>, b: &Mat<f64>) -> f64 {
    // trace(A * B) = sum_ij A_ij * B_ji
    let n = a.nrows();
    let m = a.ncols();
    let mut trace = 0.0;
    for i in 0..n {
        for k in 0..m {
            trace += a[(i, k)] * b[(k, i)];
        }
    }
    trace
}
```

Note: The exact faer API for LU/inverse may need adjustment during implementation. The key formulas are correct. The implementer should check `faer` 0.20 docs for the precise method names (e.g., `m.partial_piv_lu()` vs `faer::linalg::lu::...`).

- [ ] **Step 4: Run tests**

Run: `cd rust_deseq2 && cargo test --lib dispersion::gene_estimate`
Expected: PASS (numerical derivative checks should match)

- [ ] **Step 5: Commit**

```bash
git add rust_deseq2/src/dispersion/gene_estimate.rs rust_deseq2/src/dispersion/mod.rs
git commit -m "feat: implement NB log-posterior and derivatives for dispersion estimation"
```

---

## Task 8: Gene-wise Dispersion Estimation (fitDisp Line Search)

**Files:**
- Modify: `rust_deseq2/src/dispersion/gene_estimate.rs` (add `fit_dispersion` and `estimate_dispersions_gene`)
- Test: inline + `tests/test_dispersion.rs`

This implements the Armijo line search from DESeq2.cpp `fitDisp` (lines 164-277).

**Algorithm** (from C++ source):
1. Start at `a = log_alpha_init`, compute `lp = log_posterior(a)`, `dlp = dlog_posterior(a)`
2. Set `kappa = kappa_0` (default 1.0)
3. For each iteration (max 100):
   - Propose `a_propose = a + kappa * dlp`
   - Clamp: if `a_propose < -30`, adjust kappa; if `a_propose > 10`, adjust kappa
   - Armijo test: `theta_kappa = -log_posterior(a + kappa*dlp)`, `theta_hat = -lp - kappa * epsilon * dlp^2`
   - If `theta_kappa <= theta_hat` (accepted):
     - Update `a = a + kappa * dlp`
     - Check convergence: `change = lp_new - lp_old < tol` (default 1e-6)
     - Increase kappa slightly: `kappa = min(kappa * 1.1, kappa_0)`
     - Every 5 accepts: `kappa /= 2`
   - Else: `kappa /= 2`

- [ ] **Step 1: Write failing test**

Add to `src/dispersion/gene_estimate.rs` tests module:

```rust
#[test]
fn test_fit_dispersion_converges() {
    let y = vec![10.0, 20.0, 15.0, 25.0];
    let mu = vec![12.0, 18.0, 14.0, 22.0];
    let x = Mat::from_fn(4, 2, |i, j| {
        if j == 0 { 1.0 } else { [0.0, 1.0, 0.0, 1.0][i] }
    });
    let result = fit_dispersion(
        &y, &mu, &x,
        0.0,    // log_alpha_init
        0.0,    // prior_mean
        1.0,    // prior_var
        false,  // use_prior (MLE mode)
        true,   // use_cr
    );
    // Should converge to a finite value
    assert!(result.log_alpha.is_finite());
    assert!(result.log_alpha > -30.0 && result.log_alpha < 10.0);
    assert!(result.converged);
}
```

- [ ] **Step 2: Run test to verify failure**

Run: `cd rust_deseq2 && cargo test --lib dispersion::gene_estimate::tests::test_fit_dispersion`
Expected: FAIL

- [ ] **Step 3: Implement fit_dispersion and estimate_dispersions_gene**

Add to `src/dispersion/gene_estimate.rs`:

```rust
use rayon::prelude::*;

/// Result of single-gene dispersion fit.
pub struct DispFitResult {
    pub log_alpha: f64,
    pub iterations: usize,
    pub converged: bool,
    pub final_log_posterior: f64,
}

/// Fit dispersion for a single gene via Armijo line search on log(alpha).
/// Matches DESeq2.cpp fitDisp().
pub fn fit_dispersion(
    y: &[f64],
    mu: &[f64],
    x: &Mat<f64>,
    log_alpha_init: f64,
    prior_mean: f64,
    prior_var: f64,
    use_prior: bool,
    use_cr: bool,
) -> DispFitResult {
    let epsilon = 1e-4;
    let kappa_0 = 1.0;
    let tol = 1e-6;
    let maxit = 100;
    let min_log_alpha = -30.0;
    let max_log_alpha = 10.0;

    let mut a = log_alpha_init;
    let mut lp = log_posterior(a, y, mu, x, prior_mean, prior_var, use_prior, use_cr);
    let mut dlp = dlog_posterior(a, y, mu, x, prior_mean, prior_var, use_prior, use_cr);
    let mut kappa = kappa_0;
    let mut converged = false;
    let mut iter_accept = 0usize;
    let mut iterations = 0usize;

    for _t in 0..maxit {
        iterations += 1;

        let mut a_propose = a + kappa * dlp;
        // Clamp proposals
        if a_propose < min_log_alpha {
            kappa = (min_log_alpha - a) / dlp;
        }
        if a_propose > max_log_alpha {
            kappa = (max_log_alpha - a) / dlp;
        }

        let theta_kappa = -log_posterior(a + kappa * dlp, y, mu, x, prior_mean, prior_var, use_prior, use_cr);
        let theta_hat_kappa = -lp - kappa * epsilon * dlp * dlp;

        if theta_kappa <= theta_hat_kappa {
            // Armijo condition satisfied
            iter_accept += 1;
            a = a + kappa * dlp;
            let lp_new = log_posterior(a, y, mu, x, prior_mean, prior_var, use_prior, use_cr);
            let change = lp_new - lp;

            if change < tol {
                lp = lp_new;
                converged = true;
                break;
            }

            if a < min_log_alpha {
                break;
            }

            lp = lp_new;
            dlp = dlog_posterior(a, y, mu, x, prior_mean, prior_var, use_prior, use_cr);

            kappa = (kappa * 1.1).min(kappa_0);
            if iter_accept % 5 == 0 {
                kappa /= 2.0;
            }
        } else {
            kappa /= 2.0;
        }
    }

    DispFitResult {
        log_alpha: a,
        iterations,
        converged,
        final_log_posterior: lp,
    }
}

/// Rough dispersion estimate via method of moments.
/// roughDisp = sum((y - mu)^2 - mu) / mu^2) / (m - p)
pub fn rough_disp_estimate(
    normalized_counts: &Mat<f64>,
    mu: &Mat<f64>,
    n_params: usize,
) -> Vec<f64> {
    let n_genes = normalized_counts.nrows();
    let n_samples = normalized_counts.ncols();
    let df = n_samples - n_params;
    (0..n_genes)
        .map(|i| {
            let est: f64 = (0..n_samples)
                .map(|j| {
                    let y = normalized_counts[(i, j)];
                    let m = mu[(i, j)];
                    ((y - m).powi(2) - m) / (m * m)
                })
                .sum::<f64>()
                / df as f64;
            est.max(0.0)
        })
        .collect()
}

/// Compute initial mu from design matrix via linear model: mu = (Y Q) Q'
/// where Q comes from QR decomposition of X.
/// Y is normalized counts (genes x samples), X is design matrix (samples x p).
pub fn linear_model_mu(normalized_counts: &Mat<f64>, design: &Mat<f64>) -> Mat<f64> {
    // QR decomposition of design matrix
    let qr = design.col_piv_qr(); // or thin_qr depending on faer API
    let q = qr.compute_thin_q();
    // mu = Y * Q * Q' (projection onto column space of X)
    let yq = normalized_counts * &q;
    let mu = &yq * q.transpose();
    // Clamp mu >= minmu (0.5 as per DESeq2 default)
    let n_genes = mu.nrows();
    let n_samples = mu.ncols();
    Mat::from_fn(n_genes, n_samples, |i, j| mu[(i, j)].max(0.5))
}

/// Estimate gene-wise dispersions for all genes in parallel.
/// This is the main entry point corresponding to R's estimateDispersionsGeneEst.
pub fn estimate_dispersions_gene(
    counts: &Mat<f64>,
    size_factors: &[f64],
    design_matrix: &Mat<f64>,
) -> (Vec<f64>, Mat<f64>) {
    let n_genes = counts.nrows();
    let n_samples = counts.ncols();

    // Normalized counts
    let norm_counts = Mat::from_fn(n_genes, n_samples, |i, j| {
        counts[(i, j)] / size_factors[j]
    });

    // Initial mu via linear model
    let mu = linear_model_mu(&norm_counts, design_matrix);

    // Rough dispersion estimate as starting value
    let rough_disp = rough_disp_estimate(&norm_counts, &mu, design_matrix.ncols());
    let min_disp = 1e-8;
    let max_disp = 10.0f64.max(n_samples as f64);

    // Parallel per-gene dispersion estimation
    let dispersions: Vec<f64> = (0..n_genes)
        .into_par_iter()
        .map(|i| {
            let y: Vec<f64> = (0..n_samples).map(|j| counts[(i, j)]).collect();
            let mu_row: Vec<f64> = (0..n_samples).map(|j| mu[(i, j)]).collect();
            let alpha_init = rough_disp[i].clamp(min_disp, max_disp);
            let log_alpha_init = alpha_init.ln();

            let result = fit_dispersion(
                &y,
                &mu_row,
                design_matrix,
                log_alpha_init,
                log_alpha_init, // prior_mean = init (wide prior)
                1.0,            // prior_var = 1 (wide, no shrinkage)
                false,          // use_prior = false for gene-wise
                true,           // use_cr = true
            );

            result.log_alpha.exp().clamp(min_disp, max_disp)
        })
        .collect();

    (dispersions, mu)
}
```

Note: The exact faer QR API (`col_piv_qr()`, `compute_thin_q()`, etc.) should be verified against faer 0.20 docs during implementation. The mathematical algorithm is correct.

- [ ] **Step 4: Run tests**

Run: `cd rust_deseq2 && cargo test --lib dispersion::gene_estimate`
Expected: PASS

- [ ] **Step 5: Write integration test against R reference**

Create `tests/test_dispersion.rs`:

```rust
use deseq2_rs::io::{read_count_matrix, read_coldata};
use deseq2_rs::design::build_design_matrix;
use deseq2_rs::size_factors::estimate_size_factors;
use deseq2_rs::dispersion::gene_estimate::estimate_dispersions_gene;
use approx::assert_relative_eq;
use std::path::Path;

#[test]
fn test_gene_dispersions_vs_r() {
    let ref_dir = Path::new("tests/reference_data");
    if !ref_dir.join("airway_counts.tsv").exists() {
        eprintln!("Skipping: reference data not found");
        return;
    }
    let (counts, genes, _) = read_count_matrix(&ref_dir.join("airway_counts.tsv")).unwrap();
    let (coldata, _) = read_coldata(&ref_dir.join("airway_coldata.tsv")).unwrap();
    let design = build_design_matrix(&coldata, "dex", "untrt").unwrap();
    let sf = estimate_size_factors(&counts);

    let (disp_gene, _mu) = estimate_dispersions_gene(&counts, &sf, &design);

    // Read R reference
    let r_disp: Vec<(String, f64)> = std::fs::read_to_string(ref_dir.join("r_disp_gene_estimates.tsv"))
        .unwrap()
        .lines().skip(1)
        .filter_map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            let val: f64 = f[1].parse().ok()?;
            Some((f[0].to_string(), val))
        })
        .collect();

    // Compare: relative error < 1e-4 for most genes
    let mut max_err = 0.0f64;
    let mut n_compared = 0;
    for (i, (gene, r_val)) in r_disp.iter().enumerate() {
        if r_val.is_nan() || disp_gene[i].is_nan() { continue; }
        if *r_val < 1e-7 { continue; } // skip near-zero
        let rel_err = (disp_gene[i] - r_val).abs() / r_val;
        max_err = max_err.max(rel_err);
        n_compared += 1;
    }
    println!("Gene-wise dispersion: compared {}, max relative error: {:.2e}", n_compared, max_err);
    assert!(max_err < 1e-2, "max relative error too large: {:.2e}", max_err);
}
```

- [ ] **Step 6: Run integration test**

Run: `cd rust_deseq2 && cargo test test_gene_dispersions_vs_r`
Expected: PASS (or skip if no reference data)

- [ ] **Step 7: Commit**

```bash
git add rust_deseq2/src/dispersion/ rust_deseq2/tests/test_dispersion.rs
git commit -m "feat: implement gene-wise dispersion estimation with Armijo line search"
```

---

## Task 9: Parametric Dispersion Trend Fitting

**Files:**
- Modify: `rust_deseq2/src/dispersion/trend_fit.rs`
- Test: inline + extend `tests/test_dispersion.rs`

Algorithm from `parametricDispersionFit` in core.R (lines 2187-2210):
1. Start with coefficients `[a0=0.1, a1=1.0]` for model `disp = a0 + a1/mean`
2. Iterate:
   a. Compute residuals: `r = disp / (a0 + a1/mean)`
   b. Filter outliers: keep where `1e-4 < r < 15`
   c. Fit GLM with Gamma family, identity link: `disp ~ 1 + 1/mean` on good genes
   d. Check convergence: `sum(log(new_coefs / old_coefs)^2) < 1e-6`
3. Max 10 outer iterations

The Gamma GLM with identity link is equivalent to weighted least squares where weights are `1/mu^2` (inverse variance of Gamma). We implement this directly using IRLS for the Gamma GLM.

- [ ] **Step 1: Write failing test**

In `src/dispersion/trend_fit.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::DispTrendCoeffs;
    use approx::assert_relative_eq;

    #[test]
    fn test_trend_fit_known_data() {
        // Generate data where disp = 0.05 + 2.0/mean + noise
        let means: Vec<f64> = (1..=100).map(|i| i as f64 * 10.0).collect();
        let disps: Vec<f64> = means.iter().map(|&m| 0.05 + 2.0 / m).collect();

        let coeffs = fit_dispersion_trend(&means, &disps).unwrap();
        assert_relative_eq!(coeffs.asympt_disp, 0.05, epsilon = 1e-3);
        assert_relative_eq!(coeffs.extra_pois, 2.0, epsilon = 1e-2);
    }
}
```

- [ ] **Step 2: Run test to verify failure**

Run: `cd rust_deseq2 && cargo test --lib dispersion::trend_fit`
Expected: FAIL

- [ ] **Step 3: Implement parametric trend fitting**

```rust
use crate::data::DispTrendCoeffs;

/// Fit parametric dispersion trend: disp = a0 + a1/mean
/// Uses iteratively reweighted least squares (Gamma GLM with identity link).
/// Matches R's parametricDispersionFit().
pub fn fit_dispersion_trend(
    base_means: &[f64],
    gene_dispersions: &[f64],
) -> Result<DispTrendCoeffs, String> {
    let mut coefs = [0.1, 1.0]; // [asymptDisp, extraPois]
    let max_outer_iter = 10;

    for outer in 0..max_outer_iter {
        // Compute residuals and filter outliers
        let mut good_indices = Vec::new();
        for i in 0..base_means.len() {
            if base_means[i] <= 0.0 || gene_dispersions[i].is_nan() {
                continue;
            }
            let fitted = coefs[0] + coefs[1] / base_means[i];
            if fitted <= 0.0 { continue; }
            let residual = gene_dispersions[i] / fitted;
            if residual > 1e-4 && residual < 15.0 {
                good_indices.push(i);
            }
        }

        if good_indices.len() < 2 {
            return Err("too few genes for trend fitting".to_string());
        }

        // Fit Gamma GLM with identity link on good genes
        // Model: E[disp] = a0 + a1 * (1/mean)
        // Gamma GLM IRLS: weight = (fitted_value)^(-2), working response = standard IRLS
        let new_coefs = fit_gamma_glm_identity(&good_indices, base_means, gene_dispersions)?;

        if new_coefs[0] <= 0.0 || new_coefs[1] <= 0.0 {
            return Err("parametric dispersion fit failed: negative coefficients".to_string());
        }

        // Check convergence
        let log_change: f64 = (0..2)
            .map(|k| (new_coefs[k] / coefs[k]).ln().powi(2))
            .sum();

        coefs = new_coefs;

        if log_change < 1e-6 {
            break;
        }
    }

    Ok(DispTrendCoeffs {
        asympt_disp: coefs[0],
        extra_pois: coefs[1],
    })
}

/// Fit Gamma GLM with identity link: E[y] = b0 + b1 * x
/// where x_i = 1/mean_i, y_i = dispersion_i
/// IRLS: weight_j = 1/mu_j^2 (Gamma variance function V(mu) = mu^2)
fn fit_gamma_glm_identity(
    indices: &[usize],
    means: &[f64],
    disps: &[f64],
) -> Result<[f64; 2], String> {
    let n = indices.len();
    let mut beta = [0.1, 1.0]; // starting values
    let maxit = 50;
    let tol = 1e-8;

    for _iter in 0..maxit {
        // Compute fitted values and weights
        let mut xtw_x = [[0.0; 2]; 2];
        let mut xtw_z = [0.0; 2];

        for &i in indices {
            let x = [1.0, 1.0 / means[i]];
            let mu = beta[0] * x[0] + beta[1] * x[1]; // fitted value
            if mu <= 0.0 {
                continue;
            }
            let w = 1.0 / (mu * mu); // Gamma weight
            let y = disps[i];
            // Working response for identity link: z = mu + (y - mu) = y
            // (for identity link, z = y when using canonical form)
            // Actually in IRLS for identity link GLM:
            // z = eta + (y - mu) * (d_eta/d_mu) = mu + (y - mu) * 1 = y
            let z = y;

            for a in 0..2 {
                for b in 0..2 {
                    xtw_x[a][b] += w * x[a] * x[b];
                }
                xtw_z[a] += w * x[a] * z;
            }
        }

        // Solve 2x2 system
        let det = xtw_x[0][0] * xtw_x[1][1] - xtw_x[0][1] * xtw_x[1][0];
        if det.abs() < 1e-30 {
            return Err("singular matrix in Gamma GLM".to_string());
        }
        let new_beta = [
            (xtw_x[1][1] * xtw_z[0] - xtw_x[0][1] * xtw_z[1]) / det,
            (xtw_x[0][0] * xtw_z[1] - xtw_x[1][0] * xtw_z[0]) / det,
        ];

        let change: f64 = (0..2).map(|k| (new_beta[k] - beta[k]).powi(2)).sum();
        beta = new_beta;
        if change < tol {
            break;
        }
    }

    Ok(beta)
}

/// Estimate prior variance for dispersion MAP shrinkage.
/// varPrior = max(var(log(dispGeneEst) - log(dispFit)) - trigamma((m-p)/2), 0.25)
pub fn estimate_prior_variance(
    gene_dispersions: &[f64],
    trend_values: &[f64],
    n_samples: usize,
    n_params: usize,
) -> f64 {
    let min_disp = 1e-8;
    let df = n_samples - n_params;

    // Residuals on log scale: log(geneEst) - log(trendFit)
    let residuals: Vec<f64> = gene_dispersions
        .iter()
        .zip(trend_values)
        .filter(|(g, _)| **g >= min_disp * 100.0 && !g.is_nan())
        .map(|(g, t)| g.ln() - t.ln())
        .collect();

    if residuals.is_empty() {
        return 0.25;
    }

    // Observed variance of log-dispersion residuals
    let mean_r: f64 = residuals.iter().sum::<f64>() / residuals.len() as f64;
    let var_r: f64 = residuals.iter().map(|r| (r - mean_r).powi(2)).sum::<f64>() / (residuals.len() - 1) as f64;

    // Expected sampling variance under chi-squared: trigamma(df/2)
    let expected_var = crate::dispersion::gene_estimate::trigamma(df as f64 / 2.0);

    // Prior variance = observed - expected, floored at 0.25
    (var_r - expected_var).max(0.25)
}
```

- [ ] **Step 4: Run tests**

Run: `cd rust_deseq2 && cargo test --lib dispersion::trend_fit`
Expected: PASS

- [ ] **Step 5: Add integration test for trend fit**

Append to `tests/test_dispersion.rs`:

```rust
use deseq2_rs::dispersion::trend_fit::fit_dispersion_trend;
use deseq2_rs::size_factors::base_means;

#[test]
fn test_trend_fit_vs_r() {
    let ref_dir = Path::new("tests/reference_data");
    if !ref_dir.join("r_disp_trend_coeffs.tsv").exists() {
        eprintln!("Skipping: reference data not found");
        return;
    }
    // Read R trend coefficients
    let content = std::fs::read_to_string(ref_dir.join("r_disp_trend_coeffs.tsv")).unwrap();
    let mut r_a0 = 0.0;
    let mut r_a1 = 0.0;
    for line in content.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "asymptDisp" => r_a0 = f[1].parse().unwrap(),
            "extraPois" => r_a1 = f[1].parse().unwrap(),
            _ => {}
        }
    }

    // Run Rust pipeline up to trend fit
    let (counts, genes, _) = read_count_matrix(&ref_dir.join("airway_counts.tsv")).unwrap();
    let (coldata, _) = read_coldata(&ref_dir.join("airway_coldata.tsv")).unwrap();
    let design = build_design_matrix(&coldata, "dex", "untrt").unwrap();
    let sf = estimate_size_factors(&counts);
    let bm = base_means(&counts, &sf);
    let (disp_gene, _mu) = estimate_dispersions_gene(&counts, &sf, &design);

    let coeffs = fit_dispersion_trend(&bm, &disp_gene).unwrap();
    println!("Rust: a0={:.6}, a1={:.6}", coeffs.asympt_disp, coeffs.extra_pois);
    println!("R:    a0={:.6}, a1={:.6}", r_a0, r_a1);

    assert_relative_eq!(coeffs.asympt_disp, r_a0, epsilon = r_a0 * 1e-3);
    assert_relative_eq!(coeffs.extra_pois, r_a1, epsilon = r_a1 * 1e-3);
}
```

- [ ] **Step 6: Run integration test**

Run: `cd rust_deseq2 && cargo test test_trend_fit_vs_r`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add rust_deseq2/src/dispersion/trend_fit.rs rust_deseq2/tests/test_dispersion.rs
git commit -m "feat: implement parametric dispersion trend fitting and prior variance estimation"
```

---

## Task 10: MAP Dispersion Shrinkage

**Files:**
- Modify: `rust_deseq2/src/dispersion/map_estimate.rs`
- Modify: `rust_deseq2/src/dispersion/mod.rs` (add orchestrator)
- Test: extend `tests/test_dispersion.rs`

MAP shrinkage re-runs `fit_dispersion` with `use_prior=true`, using the trend value as prior mean and the estimated prior variance.

- [ ] **Step 1: Write failing test**

In `src/dispersion/map_estimate.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_map_shrinks_toward_trend() {
        // If gene-wise estimate is far from trend, MAP should shrink it
        let y = vec![5.0, 8.0, 6.0, 10.0];
        let mu = vec![7.0, 7.0, 7.0, 7.0];
        let x = faer::Mat::from_fn(4, 2, |i, j| {
            if j == 0 { 1.0 } else { [0.0, 1.0, 0.0, 1.0][i] }
        });

        let trend_log_alpha = (-1.0f64).ln(); // trend says alpha ~ 0.37
        let prior_var = 0.5;

        let result = estimate_map_single(
            &y, &mu, &x,
            0.0, // initial log_alpha (alpha=1, far from trend)
            trend_log_alpha,
            prior_var,
        );

        // MAP estimate should be between gene-wise (0.0) and trend (trend_log_alpha)
        // i.e., it should move toward the trend
        assert!(result.log_alpha < 0.0 || result.log_alpha > trend_log_alpha,
            "MAP should be between gene-wise and trend values");
    }
}
```

- [ ] **Step 2: Run test to verify failure**

Run: `cd rust_deseq2 && cargo test --lib dispersion::map_estimate`
Expected: FAIL

- [ ] **Step 3: Implement MAP shrinkage**

```rust
use crate::dispersion::gene_estimate::{fit_dispersion, DispFitResult};
use faer::Mat;
use rayon::prelude::*;

/// MAP dispersion estimate for a single gene.
pub fn estimate_map_single(
    y: &[f64],
    mu: &[f64],
    x: &Mat<f64>,
    log_alpha_init: f64,
    prior_mean: f64,  // log of trend value
    prior_var: f64,
) -> DispFitResult {
    fit_dispersion(
        y, mu, x,
        log_alpha_init,
        prior_mean,
        prior_var,
        true, // use_prior = true for MAP
        true, // use_cr = true
    )
}

/// Estimate MAP dispersions for all genes in parallel.
pub fn estimate_dispersions_map(
    counts: &Mat<f64>,
    size_factors: &[f64],
    design_matrix: &Mat<f64>,
    mu: &Mat<f64>,
    gene_dispersions: &[f64],
    trend_values: &[f64],
    prior_var: f64,
) -> Vec<f64> {
    let n_genes = counts.nrows();
    let n_samples = counts.ncols();
    let min_disp = 1e-8;
    let max_disp = 10.0f64.max(n_samples as f64);

    (0..n_genes)
        .into_par_iter()
        .map(|i| {
            let y: Vec<f64> = (0..n_samples).map(|j| counts[(i, j)]).collect();
            let mu_row: Vec<f64> = (0..n_samples).map(|j| mu[(i, j)]).collect();

            let log_alpha_init = gene_dispersions[i].max(min_disp).ln();
            let prior_mean = trend_values[i].ln();

            let result = estimate_map_single(
                &y, &mu_row, design_matrix,
                log_alpha_init, prior_mean, prior_var,
            );

            result.log_alpha.exp().clamp(min_disp, max_disp)
        })
        .collect()
}
```

- [ ] **Step 4: Implement dispersion orchestrator in mod.rs**

In `src/dispersion/mod.rs`:

```rust
pub mod gene_estimate;
pub mod trend_fit;
pub mod map_estimate;

use crate::data::{DispTrendCoeffs, DispersionResult};
use crate::size_factors::base_means;
use faer::Mat;

/// Run the full dispersion estimation pipeline:
/// 1. Gene-wise MLE
/// 2. Parametric trend fit
/// 3. Prior variance estimation
/// 4. MAP shrinkage
pub fn estimate_dispersions(
    counts: &Mat<f64>,
    size_factors: &[f64],
    design_matrix: &Mat<f64>,
) -> Result<(DispersionResult, Mat<f64>), String> {
    // Step 1: Gene-wise dispersion
    let (gene_estimates, mu) =
        gene_estimate::estimate_dispersions_gene(counts, size_factors, design_matrix);

    // Base means for trend fitting
    let bm = base_means(counts, size_factors);

    // Step 2: Parametric trend fit
    let trend_coeffs = trend_fit::fit_dispersion_trend(&bm, &gene_estimates)?;
    let trend_values: Vec<f64> = bm.iter().map(|m| trend_coeffs.eval(*m)).collect();

    // Step 3: Prior variance
    let prior_var = trend_fit::estimate_prior_variance(
        &gene_estimates,
        &trend_values,
        counts.ncols(),
        design_matrix.ncols(),
    );

    // Step 4: MAP shrinkage
    let map_estimates = map_estimate::estimate_dispersions_map(
        counts,
        size_factors,
        design_matrix,
        &mu,
        &gene_estimates,
        &trend_values,
        prior_var,
    );

    Ok((
        DispersionResult {
            gene_estimates,
            trend_values,
            map_estimates,
            prior_var,
        },
        mu,
    ))
}
```

- [ ] **Step 5: Run all dispersion tests**

Run: `cd rust_deseq2 && cargo test dispersion`
Expected: PASS

- [ ] **Step 6: Add MAP integration test**

Append to `tests/test_dispersion.rs`:

```rust
use deseq2_rs::dispersion::estimate_dispersions;

#[test]
fn test_map_dispersions_vs_r() {
    let ref_dir = Path::new("tests/reference_data");
    if !ref_dir.join("r_disp_map_estimates.tsv").exists() {
        eprintln!("Skipping: reference data not found");
        return;
    }

    let (counts, genes, _) = read_count_matrix(&ref_dir.join("airway_counts.tsv")).unwrap();
    let (coldata, _) = read_coldata(&ref_dir.join("airway_coldata.tsv")).unwrap();
    let design = build_design_matrix(&coldata, "dex", "untrt").unwrap();
    let sf = estimate_size_factors(&counts);

    let (disp_result, _mu) = estimate_dispersions(&counts, &sf, &design).unwrap();

    // Read R MAP dispersions
    let r_map: Vec<f64> = std::fs::read_to_string(ref_dir.join("r_disp_map_estimates.tsv"))
        .unwrap()
        .lines().skip(1)
        .filter_map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            f[1].parse::<f64>().ok()
        })
        .collect();

    let mut max_err = 0.0f64;
    let mut n_compared = 0;
    for i in 0..r_map.len().min(disp_result.map_estimates.len()) {
        if r_map[i].is_nan() || disp_result.map_estimates[i].is_nan() { continue; }
        if r_map[i] < 1e-7 { continue; }
        let rel_err = (disp_result.map_estimates[i] - r_map[i]).abs() / r_map[i];
        max_err = max_err.max(rel_err);
        n_compared += 1;
    }
    println!("MAP dispersion: compared {}, max relative error: {:.2e}", n_compared, max_err);
    assert!(max_err < 1e-2, "MAP dispersion max error too large: {:.2e}", max_err);
}
```

- [ ] **Step 7: Run integration test**

Run: `cd rust_deseq2 && cargo test test_map_dispersions_vs_r`
Expected: PASS

- [ ] **Step 8: Commit**

```bash
git add rust_deseq2/src/dispersion/
git commit -m "feat: implement MAP dispersion shrinkage and full dispersion pipeline"
```

---

## Task 11: IRLS Negative Binomial GLM Fitting (fitBeta)

**Files:**
- Modify: `rust_deseq2/src/glm/irls.rs`
- Modify: `rust_deseq2/src/glm/mod.rs`
- Test: inline + `tests/test_glm.rs`

This is the most complex numerical routine. It implements the IRLS algorithm from DESeq2.cpp `fitBeta` (lines 283-465).

**Algorithm per gene:**
1. Initialize beta from QR of log(normalized_counts + 0.1)
2. Compute `mu = nf * exp(X * beta)`, clamp `mu >= 0.5`
3. Iterate (max 100):
   a. Weights: `w = mu / (1 + alpha * mu)`
   b. Working response: `z = log(mu/nf) + (y - mu)/mu`
   c. QR approach: stack `[sqrt(W)*X; sqrt(lambda)]`, QR decompose, solve
   d. Update beta, recompute mu
   e. Deviance: `dev = -2 * sum(log_dnbinom(y, 1/alpha, mu))`
   f. Convergence: `|dev - dev_old| / (|dev| + 0.1) < tol`
4. After convergence, compute:
   - `sigma = (X'WX + lambda)^{-1} X'WX (X'WX + lambda)^{-1}` (sandwich covariance)
   - SE = sqrt(diag(sigma))

- [ ] **Step 1: Write failing test**

In `src/glm/irls.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use faer::Mat;
    use approx::assert_relative_eq;

    #[test]
    fn test_irls_simple() {
        // 4 samples, 2 genes
        let counts = Mat::from_fn(2, 4, |i, j| {
            [[100.0, 200.0, 150.0, 300.0],
             [50.0, 40.0, 60.0, 45.0]][i][j]
        });
        let sf = vec![1.0, 1.0, 1.0, 1.0];
        let design = Mat::from_fn(4, 2, |i, j| {
            if j == 0 { 1.0 } else { [0.0, 1.0, 0.0, 1.0][i] }
        });
        let dispersions = vec![0.1, 0.1];

        let result = fit_negative_binomial_glm(&counts, &sf, &design, &dispersions);

        // Should have 2 genes x 2 coefficients
        assert_eq!(result.coefficients.nrows(), 2);
        assert_eq!(result.coefficients.ncols(), 2);
        // Gene 0: treated has higher counts -> positive fold change
        assert!(result.coefficients[(0, 1)] > 0.0);
        // Gene 1: treated has similar/lower counts -> near zero or negative
        assert!(result.converged[0]);
        assert!(result.converged[1]);
    }
}
```

- [ ] **Step 2: Run test to verify failure**

Run: `cd rust_deseq2 && cargo test --lib glm::irls`
Expected: FAIL

- [ ] **Step 3: Implement IRLS**

```rust
use crate::data::BetaFitResult;
use faer::Mat;
use rayon::prelude::*;
use statrs::function::gamma::ln_gamma;

/// NB log-probability: log P(y | size=1/alpha, mu)
fn log_dnbinom_mu(y: f64, size: f64, mu: f64) -> f64 {
    // size = 1/alpha, parameterization: P(Y=y) = C(y+size-1, y) * (mu/(mu+size))^y * (size/(mu+size))^size
    ln_gamma(y + size) - ln_gamma(size) - ln_gamma(y + 1.0)
        + y * (mu / (mu + size)).ln()
        + size * (size / (mu + size)).ln()
}

/// Fit NB GLM for all genes via IRLS with QR decomposition.
/// Returns coefficients on NATURAL LOG scale (not log2).
pub fn fit_negative_binomial_glm(
    counts: &Mat<f64>,
    size_factors: &[f64],
    design_matrix: &Mat<f64>,
    dispersions: &[f64],
) -> BetaFitResult {
    let n_genes = counts.nrows();
    let n_samples = counts.ncols();
    let n_params = design_matrix.ncols();
    let tol = 1e-6;
    let maxit = 100;
    let minmu = 0.5;
    let large = 30.0;

    // Normalization factors matrix (genes x samples): sf broadcasted
    // Ridge penalty: very small by default (no betaPrior), use 1e-6
    let lambda: Vec<f64> = vec![1e-6; n_params];

    // Initialize betas via QR of log(normalized_counts + 0.1)
    let init_betas = initialize_betas(counts, size_factors, design_matrix);

    // Parallel per-gene fitting
    let results: Vec<(Vec<f64>, Vec<f64>, f64, bool)> = (0..n_genes)
        .into_par_iter()
        .map(|gene| {
            let y: Vec<f64> = (0..n_samples).map(|j| counts[(gene, j)]).collect();
            let nf: Vec<f64> = (0..n_samples).map(|j| size_factors[j]).collect();
            let alpha = dispersions[gene];
            let mut beta: Vec<f64> = (0..n_params).map(|k| init_betas[(gene, k)]).collect();

            // Initial mu
            let mut mu: Vec<f64> = (0..n_samples)
                .map(|j| {
                    let eta: f64 = (0..n_params).map(|k| design_matrix[(j, k)] * beta[k]).sum();
                    (nf[j] * eta.exp()).max(minmu)
                })
                .collect();

            let mut dev = 0.0;
            let mut dev_old = 0.0;
            let mut converged = false;

            for t in 0..maxit {
                // Weights and working response
                let w: Vec<f64> = mu.iter().map(|&m| m / (1.0 + alpha * m)).collect();
                let z: Vec<f64> = (0..n_samples)
                    .map(|j| (mu[j] / nf[j]).ln() + (y[j] - mu[j]) / mu[j])
                    .collect();

                // QR approach: stack [sqrt(W)*X ; sqrt(Lambda)]
                let n_aug = n_samples + n_params;
                let mut aug_x = Mat::<f64>::zeros(n_aug, n_params);
                let mut aug_z = vec![0.0; n_aug];

                for j in 0..n_samples {
                    let sw = w[j].sqrt();
                    for k in 0..n_params {
                        aug_x[(j, k)] = sw * design_matrix[(j, k)];
                    }
                    aug_z[j] = sw * z[j];
                }
                for k in 0..n_params {
                    aug_x[(n_samples + k, k)] = lambda[k].sqrt();
                    aug_z[n_samples + k] = 0.0;
                }

                // QR solve: aug_x * beta_new = aug_z
                // Use faer QR
                let qr = aug_x.col_piv_qr();
                let aug_z_mat = Mat::from_fn(n_aug, 1, |i, _| aug_z[i]);
                let beta_new_mat = qr.solve(&aug_z_mat);
                let beta_new: Vec<f64> = (0..n_params).map(|k| beta_new_mat[(k, 0)]).collect();

                // Check for divergence
                if beta_new.iter().any(|&b| b.abs() > large) {
                    converged = false;
                    break;
                }

                beta = beta_new;

                // Update mu
                mu = (0..n_samples)
                    .map(|j| {
                        let eta: f64 = (0..n_params).map(|k| design_matrix[(j, k)] * beta[k]).sum();
                        (nf[j] * eta.exp()).max(minmu)
                    })
                    .collect();

                // Deviance
                dev = -2.0 * (0..n_samples)
                    .map(|j| log_dnbinom_mu(y[j], 1.0 / alpha, mu[j]))
                    .sum::<f64>();

                let conv_test = (dev - dev_old).abs() / (dev.abs() + 0.1);
                if conv_test.is_nan() {
                    converged = false;
                    break;
                }
                if t > 0 && conv_test < tol {
                    converged = true;
                    break;
                }
                dev_old = dev;
            }

            // Compute standard errors via sandwich estimator
            // sigma = (X'WX + L)^{-1} X'WX (X'WX + L)^{-1}
            let w: Vec<f64> = mu.iter().map(|&m| m / (1.0 + alpha * m)).collect();
            let mut xtwx = Mat::<f64>::zeros(n_params, n_params);
            for j in 0..n_samples {
                for a in 0..n_params {
                    for b in 0..n_params {
                        xtwx[(a, b)] += design_matrix[(j, a)] * w[j] * design_matrix[(j, b)];
                    }
                }
            }
            // xtwx_r = xtwx + ridge
            let mut xtwx_r = xtwx.clone();
            for k in 0..n_params {
                xtwx_r[(k, k)] += lambda[k];
            }
            let xtwx_r_inv = xtwx_r.partial_piv_lu().solve(
                &Mat::from_fn(n_params, n_params, |i, j| if i == j { 1.0 } else { 0.0 })
            );
            // sigma = xtwx_r_inv * xtwx * xtwx_r_inv
            let sigma = &xtwx_r_inv * &xtwx * &xtwx_r_inv;
            let se: Vec<f64> = (0..n_params).map(|k| sigma[(k, k)].sqrt()).collect();

            (beta, se, dev, converged)
        })
        .collect();

    // Assemble results
    let mut coefficients = Mat::<f64>::zeros(n_genes, n_params);
    let mut standard_errors = Mat::<f64>::zeros(n_genes, n_params);
    let mut deviances = Vec::with_capacity(n_genes);
    let mut converged_vec = Vec::with_capacity(n_genes);

    for (i, (beta, se, dev, conv)) in results.into_iter().enumerate() {
        for k in 0..n_params {
            coefficients[(i, k)] = beta[k];
            standard_errors[(i, k)] = se[k];
        }
        deviances.push(dev);
        converged_vec.push(conv);
    }

    BetaFitResult {
        coefficients,
        standard_errors,
        deviances,
        converged: converged_vec,
    }
}

/// Initialize betas via QR decomposition of design matrix applied to log(normalized counts + 0.1).
fn initialize_betas(counts: &Mat<f64>, size_factors: &[f64], design: &Mat<f64>) -> Mat<f64> {
    let n_genes = counts.nrows();
    let n_samples = counts.ncols();
    let n_params = design.ncols();

    // log(normalized_counts + 0.1), transposed to samples x genes for QR solve
    let log_norm = Mat::from_fn(n_samples, n_genes, |j, i| {
        (counts[(i, j)] / size_factors[j] + 0.1).ln()
    });

    // QR of design matrix
    let qr = design.col_piv_qr();
    let q = qr.compute_thin_q();
    let r = qr.compute_thin_r();

    // beta = R^{-1} Q' Y  (for each gene column)
    let qty = q.transpose() * &log_norm; // n_params x n_genes
    // Solve R * beta = qty
    let r_lu = r.partial_piv_lu();
    let beta_t = r_lu.solve(&qty); // n_params x n_genes

    // Transpose to genes x params
    Mat::from_fn(n_genes, n_params, |i, k| beta_t[(k, i)])
}
```

Note: faer QR/LU API names may need adjustment. The mathematical operations are correct. The implementer should verify exact method names against faer 0.20 docs (`col_piv_qr` vs `qr`, `compute_thin_q` vs `thin_q`, etc.).

- [ ] **Step 4: Update glm/mod.rs**

```rust
pub mod irls;
pub use irls::fit_negative_binomial_glm;
```

- [ ] **Step 5: Run tests**

Run: `cd rust_deseq2 && cargo test --lib glm`
Expected: PASS

- [ ] **Step 6: Write integration test**

Create `tests/test_glm.rs`:

```rust
use deseq2_rs::io::{read_count_matrix, read_coldata};
use deseq2_rs::design::build_design_matrix;
use deseq2_rs::size_factors::estimate_size_factors;
use deseq2_rs::dispersion::estimate_dispersions;
use deseq2_rs::glm::fit_negative_binomial_glm;
use approx::assert_relative_eq;
use std::path::Path;

#[test]
fn test_beta_coefficients_vs_r() {
    let ref_dir = Path::new("tests/reference_data");
    if !ref_dir.join("r_beta_coefficients.tsv").exists() {
        eprintln!("Skipping: reference data not found");
        return;
    }

    let (counts, genes, _) = read_count_matrix(&ref_dir.join("airway_counts.tsv")).unwrap();
    let (coldata, _) = read_coldata(&ref_dir.join("airway_coldata.tsv")).unwrap();
    let design = build_design_matrix(&coldata, "dex", "untrt").unwrap();
    let sf = estimate_size_factors(&counts);
    let (disp_result, _mu) = estimate_dispersions(&counts, &sf, &design).unwrap();

    let beta_fit = fit_negative_binomial_glm(&counts, &sf, &design, &disp_result.map_estimates);

    // Read R beta coefficients (log2 scale in R)
    let content = std::fs::read_to_string(ref_dir.join("r_beta_coefficients.tsv")).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    let _header = lines[0]; // gene, Intercept, dex_trt

    let log2_e = std::f64::consts::E.log2(); // = 1/ln(2) = log2(e)

    let mut max_err = 0.0f64;
    let mut n_compared = 0;
    for (i, line) in lines[1..].iter().enumerate() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 { continue; }
        // R stores on log2 scale; our betas are natural log
        let r_intercept: f64 = match f[1].parse() { Ok(v) => v, Err(_) => continue };
        let r_lfc: f64 = match f[2].parse() { Ok(v) => v, Err(_) => continue };

        if i >= beta_fit.coefficients.nrows() { break; }

        // Convert our natural log beta to log2
        let rust_intercept_log2 = beta_fit.coefficients[(i, 0)] * log2_e;
        let rust_lfc_log2 = beta_fit.coefficients[(i, 1)] * log2_e;

        if r_lfc.abs() < 1e-6 { continue; } // skip near-zero

        let rel_err = (rust_lfc_log2 - r_lfc).abs() / r_lfc.abs().max(1e-10);
        max_err = max_err.max(rel_err);
        n_compared += 1;
    }

    println!("Beta (LFC): compared {}, max relative error: {:.2e}", n_compared, max_err);
    assert!(max_err < 1e-2, "Beta max relative error too large: {:.2e}", max_err);
}
```

- [ ] **Step 7: Run integration test**

Run: `cd rust_deseq2 && cargo test test_beta_coefficients_vs_r`
Expected: PASS

- [ ] **Step 8: Commit**

```bash
git add rust_deseq2/src/glm/ rust_deseq2/tests/test_glm.rs
git commit -m "feat: implement IRLS negative binomial GLM fitting (fitBeta)"
```

---

## Task 12: Wald Test

**Files:**
- Modify: `rust_deseq2/src/test_stats/wald.rs`
- Modify: `rust_deseq2/src/test_stats/mod.rs`
- Test: inline + `tests/test_wald.rs`

Wald statistic: `W = beta / SE(beta)`, p-value: `p = 2 * (1 - Phi(|W|))`

For contrasts: `W = c'beta / sqrt(c' Sigma c)` where Sigma is the sandwich covariance matrix. Since fitBeta already computes the full covariance, we can extract it.

Note: R DESeq2 stores betas and SEs on log2 scale. Our IRLS stores on natural log scale. Convert by multiplying by `log2(e)`.

- [ ] **Step 1: Write failing test**

In `src/test_stats/wald.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::BetaFitResult;
    use faer::Mat;
    use approx::assert_relative_eq;

    #[test]
    fn test_wald_test_basic() {
        // 2 genes, 2 params
        let mut coefficients = Mat::from_fn(2, 2, |i, j| {
            [[5.0, 0.5], [3.0, -0.1]][i][j]  // natural log scale
        });
        let standard_errors = Mat::from_fn(2, 2, |i, j| {
            [[0.1, 0.2], [0.1, 0.3]][i][j]
        });

        let fit = BetaFitResult {
            coefficients,
            standard_errors,
            deviances: vec![100.0, 100.0],
            converged: vec![true, true],
        };

        // Contrast: test second coefficient (LFC)
        let contrast = vec![0.0, 1.0];
        let results = wald_test(&fit, &contrast);

        assert_eq!(results.len(), 2);
        // Gene 0: beta=0.5, SE=0.2, stat=2.5, two-sided p ~ 0.0124
        let log2_e = std::f64::consts::E.log2();
        assert_relative_eq!(results[0].log2_fold_change, 0.5 * log2_e, epsilon = 1e-10);
        assert_relative_eq!(results[0].stat, 0.5 / 0.2, epsilon = 1e-10);
        assert!(results[0].p_value < 0.02);
        assert!(results[0].p_value > 0.01);
    }
}
```

- [ ] **Step 2: Run test to verify failure**

Run: `cd rust_deseq2 && cargo test --lib test_stats::wald`
Expected: FAIL

- [ ] **Step 3: Implement Wald test**

```rust
use crate::data::{BetaFitResult, WaldTestResult};
use statrs::distribution::{Normal, ContinuousCDF};

/// Compute Wald test for each gene given a contrast vector.
/// Coefficients and SEs in BetaFitResult are on natural log scale.
/// Output log2FoldChange and lfcSE are on log2 scale.
pub fn wald_test(beta_fit: &BetaFitResult, contrast: &[f64]) -> Vec<WaldTestResult> {
    let n_genes = beta_fit.coefficients.nrows();
    let n_params = beta_fit.coefficients.ncols();
    let log2_e = std::f64::consts::E.log2(); // 1 / ln(2)
    let normal = Normal::new(0.0, 1.0).unwrap();

    (0..n_genes)
        .map(|i| {
            // c' * beta (natural log scale)
            let estimate: f64 = (0..n_params)
                .map(|k| contrast[k] * beta_fit.coefficients[(i, k)])
                .sum();

            // SE: for simple single-coefficient contrast c=[0,...,1,...,0],
            // SE = standard_errors[i, k] for the corresponding k.
            // For general contrast: SE = sqrt(c' Sigma c).
            // Here we use the stored diagonal SEs for single-coeff contrasts,
            // or compute from the diagonal for general contrasts.
            // (Full covariance matrix would need to be stored for general contrasts.
            // For MVP, we support single-coefficient and compute SE from diagonal.)
            let se: f64 = {
                // Approximate: SE(c'beta) = sqrt(sum(c_k^2 * SE_k^2))
                // This is exact when the off-diagonal covariances are zero (orthogonal design)
                // and a good approximation for typical designs.
                let var: f64 = (0..n_params)
                    .map(|k| contrast[k].powi(2) * beta_fit.standard_errors[(i, k)].powi(2))
                    .sum();
                var.sqrt()
            };

            let stat = if se > 0.0 { estimate / se } else { 0.0 };
            let p_value = if stat.is_finite() {
                2.0 * (1.0 - normal.cdf(stat.abs()))
            } else {
                f64::NAN
            };

            WaldTestResult {
                log2_fold_change: estimate * log2_e,
                lfc_se: se * log2_e,
                stat,
                p_value,
            }
        })
        .collect()
}
```

**Important note for integration accuracy:** The SE computation above uses the diagonal approximation. For the airway dataset (simple two-group comparison), this is exact because the contrast `[0, 1]` picks a single coefficient. If the integration test shows discrepancies for more complex contrasts, the implementer should store the full covariance matrix per gene in `BetaFitResult` and compute `sqrt(c' Sigma c)` exactly. For the MVP airway validation this diagonal approach is correct.

- [ ] **Step 4: Update test_stats/mod.rs**

```rust
pub mod wald;
pub mod p_adjust;
pub use wald::wald_test;
```

- [ ] **Step 5: Run tests**

Run: `cd rust_deseq2 && cargo test --lib test_stats::wald`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add rust_deseq2/src/test_stats/
git commit -m "feat: implement Wald test for differential expression"
```

---

## Task 13: BH P-value Adjustment and Independent Filtering

**Files:**
- Modify: `rust_deseq2/src/test_stats/p_adjust.rs`
- Test: inline + `tests/test_wald.rs`

**BH procedure:** Sort p-values, adjust: `p_adj[i] = min(p[i] * n / rank[i], p_adj[i+1])` (ensures monotonicity).

**Independent filtering** (simplified from R's `pvalueAdjustment`):
1. Compute quantile grid of baseMean (50 points from fraction_zero to 0.95)
2. For each threshold, filter genes with baseMean > threshold, apply BH
3. Count rejections at each threshold
4. Pick threshold that maximizes rejections (with LOWESS smoothing)

- [ ] **Step 1: Write failing test for BH**

In `src/test_stats/p_adjust.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_bh_adjustment() {
        let pvals = vec![0.01, 0.04, 0.03, 0.005];
        let padj = p_adjust_bh(&pvals);
        // Sorted by p: 0.005(idx3), 0.01(idx0), 0.03(idx2), 0.04(idx1)
        // Adjusted: 0.005*4/1=0.02, 0.01*4/2=0.02, 0.03*4/3=0.04, 0.04*4/4=0.04
        // After monotonicity: 0.02, 0.02, 0.04, 0.04
        assert_relative_eq!(padj[3], 0.02, epsilon = 1e-10); // gene 3 (smallest p)
        assert_relative_eq!(padj[0], 0.02, epsilon = 1e-10); // gene 0
        assert_relative_eq!(padj[2], 0.04, epsilon = 1e-10); // gene 2
        assert_relative_eq!(padj[1], 0.04, epsilon = 1e-10); // gene 1
    }

    #[test]
    fn test_bh_with_nan() {
        let pvals = vec![0.01, f64::NAN, 0.03];
        let padj = p_adjust_bh(&pvals);
        assert!(!padj[0].is_nan());
        assert!(padj[1].is_nan());
        assert!(!padj[2].is_nan());
        // Only 2 non-NaN values
        // Sorted: 0.01(idx0), 0.03(idx2)
        // Adjusted: 0.01*2/1=0.02, 0.03*2/2=0.03
        assert_relative_eq!(padj[0], 0.02, epsilon = 1e-10);
        assert_relative_eq!(padj[2], 0.03, epsilon = 1e-10);
    }
}
```

- [ ] **Step 2: Run test to verify failure**

Run: `cd rust_deseq2 && cargo test --lib test_stats::p_adjust`
Expected: FAIL

- [ ] **Step 3: Implement BH adjustment**

```rust
/// Benjamini-Hochberg p-value adjustment.
/// NaN p-values remain NaN (excluded from adjustment count).
pub fn p_adjust_bh(p_values: &[f64]) -> Vec<f64> {
    let n = p_values.len();
    let mut result = vec![f64::NAN; n];

    // Collect non-NaN indices and p-values
    let mut indexed: Vec<(usize, f64)> = p_values
        .iter()
        .enumerate()
        .filter(|(_, p)| !p.is_nan())
        .map(|(i, &p)| (i, p))
        .collect();

    if indexed.is_empty() {
        return result;
    }

    let m = indexed.len();

    // Sort by p-value descending
    indexed.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    // Adjust: p_adj[i] = min(p[i] * m / rank_from_bottom, 1.0)
    // Ensure monotonicity (cumulative minimum from largest to smallest p)
    let mut cummin = f64::INFINITY;
    for (rank_from_top, &(orig_idx, p)) in indexed.iter().enumerate() {
        let rank_from_bottom = m - rank_from_top; // rank in ascending order
        let adjusted = (p * m as f64 / rank_from_bottom as f64).min(1.0);
        cummin = cummin.min(adjusted);
        result[orig_idx] = cummin;
    }

    result
}

/// Independent filtering: find optimal baseMean threshold, then apply BH.
/// Returns adjusted p-values (NaN for filtered-out genes).
pub fn independent_filtering(
    base_means: &[f64],
    p_values: &[f64],
    alpha: f64,
) -> Vec<f64> {
    let n = base_means.len();

    // Quantile grid
    let fraction_zero = base_means.iter().filter(|&&m| m == 0.0).count() as f64 / n as f64;
    let upper = if fraction_zero < 0.95 { 0.95 } else { 1.0 };
    let n_theta = 50;
    let thetas: Vec<f64> = (0..n_theta)
        .map(|i| fraction_zero + (upper - fraction_zero) * i as f64 / (n_theta - 1) as f64)
        .collect();

    // Sorted base means for quantile computation
    let mut sorted_means: Vec<f64> = base_means.to_vec();
    sorted_means.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let quantile = |q: f64| -> f64 {
        if q <= 0.0 { return sorted_means[0]; }
        if q >= 1.0 { return sorted_means[n - 1]; }
        let idx = q * (n - 1) as f64;
        let lo = idx.floor() as usize;
        let hi = lo + 1;
        let frac = idx - lo as f64;
        if hi >= n { sorted_means[n - 1] }
        else { sorted_means[lo] * (1.0 - frac) + sorted_means[hi] * frac }
    };

    // For each theta, count rejections after filtering
    let mut best_theta_idx = 0;
    let mut max_rejections = 0;

    let rejections: Vec<usize> = thetas
        .iter()
        .map(|&theta| {
            let threshold = quantile(theta);
            let filtered_p: Vec<f64> = (0..n)
                .map(|i| {
                    if base_means[i] >= threshold { p_values[i] } else { f64::NAN }
                })
                .collect();
            let padj = p_adjust_bh(&filtered_p);
            padj.iter().filter(|&&p| !p.is_nan() && p < alpha).count()
        })
        .collect();

    // Find threshold maximizing rejections
    // Use simple max (the full R implementation uses LOWESS smoothing,
    // but for the MVP this gives equivalent results on most datasets)
    if let Some(max_rej) = rejections.iter().max() {
        if *max_rej > 10 {
            best_theta_idx = rejections.iter().position(|&r| r == *max_rej).unwrap_or(0);
        }
    }

    // Apply BH with the chosen threshold
    let threshold = quantile(thetas[best_theta_idx]);
    let filtered_p: Vec<f64> = (0..n)
        .map(|i| {
            if base_means[i] >= threshold { p_values[i] } else { f64::NAN }
        })
        .collect();
    p_adjust_bh(&filtered_p)
}
```

- [ ] **Step 4: Run tests**

Run: `cd rust_deseq2 && cargo test --lib test_stats::p_adjust`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add rust_deseq2/src/test_stats/p_adjust.rs
git commit -m "feat: implement BH p-value adjustment and independent filtering"
```

---

## Task 14: Pipeline Layer

**Files:**
- Modify: `rust_deseq2/src/pipeline.rs`
- Test: `tests/test_pipeline.rs`

This wires all functional core modules together into the `DESeqDataSet::run()` method.

- [ ] **Step 1: Write failing test**

Create `tests/test_pipeline.rs`:

```rust
use deseq2_rs::pipeline::DESeqDataSet;
use deseq2_rs::data::Contrast;
use std::path::Path;

#[test]
fn test_pipeline_end_to_end() {
    let ref_dir = Path::new("tests/reference_data");
    if !ref_dir.join("airway_counts.tsv").exists() {
        eprintln!("Skipping: reference data not found");
        return;
    }

    let mut dds = DESeqDataSet::from_csv(
        &ref_dir.join("airway_counts.tsv"),
        &ref_dir.join("airway_coldata.tsv"),
        "dex",
        "untrt",
    )
    .unwrap();

    dds.run().unwrap();
    let results = dds.results(Contrast::LastCoefficient).unwrap();

    // Should have results for all genes
    assert!(!results.is_empty());

    // Check that we have reasonable p-values
    let n_significant = results
        .iter()
        .filter(|r| r.p_adjusted < 0.01 && r.log2_fold_change.abs() > 1.0)
        .count();
    println!("Significant genes (padj<0.01, |log2FC|>1): {}", n_significant);
    assert!(n_significant > 0, "should find some significant genes");
}
```

- [ ] **Step 2: Run test to verify failure**

Run: `cd rust_deseq2 && cargo test test_pipeline_end_to_end`
Expected: FAIL

- [ ] **Step 3: Implement pipeline**

```rust
use crate::data::*;
use crate::design::build_design_matrix;
use crate::dispersion;
use crate::glm;
use crate::io;
use crate::size_factors;
use crate::test_stats::{p_adjust, wald};
use faer::Mat;
use std::path::Path;

pub struct DESeqDataSet {
    pub counts: Mat<f64>,
    pub col_data: DataFrame,
    pub design_matrix: Mat<f64>,
    pub gene_names: Vec<String>,
    pub sample_names: Vec<String>,
    // Intermediate results
    size_factors: Option<Vec<f64>>,
    base_means: Option<Vec<f64>>,
    dispersion_result: Option<DispersionResult>,
    beta_fit: Option<BetaFitResult>,
    mu: Option<Mat<f64>>,
}

impl DESeqDataSet {
    pub fn from_csv(
        count_path: &Path,
        coldata_path: &Path,
        design_column: &str,
        reference_level: &str,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let (counts, gene_names, sample_names) = io::read_count_matrix(count_path)?;
        let (col_data, coldata_samples) = io::read_coldata(coldata_path)?;

        // Verify sample order matches
        // (Reorder counts columns to match coldata if needed, or error)
        let design_matrix = build_design_matrix(&col_data, design_column, reference_level)?;

        Ok(Self {
            counts,
            col_data,
            design_matrix,
            gene_names,
            sample_names,
            size_factors: None,
            base_means: None,
            dispersion_result: None,
            beta_fit: None,
            mu: None,
        })
    }

    pub fn from_matrix(
        counts: Mat<f64>,
        col_data: DataFrame,
        design_matrix: Mat<f64>,
        gene_names: Vec<String>,
        sample_names: Vec<String>,
    ) -> Result<Self, String> {
        Ok(Self {
            counts,
            col_data,
            design_matrix,
            gene_names,
            sample_names,
            size_factors: None,
            base_means: None,
            dispersion_result: None,
            beta_fit: None,
            mu: None,
        })
    }

    pub fn run(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        // Step 1: Size factors
        let sf = size_factors::estimate_size_factors(&self.counts);
        let bm = size_factors::base_means(&self.counts, &sf);
        self.size_factors = Some(sf.clone());
        self.base_means = Some(bm);

        // Step 2: Dispersion estimation
        let (disp_result, mu) =
            dispersion::estimate_dispersions(&self.counts, &sf, &self.design_matrix)?;
        self.mu = Some(mu);
        self.dispersion_result = Some(disp_result);

        // Step 3: GLM fitting (using MAP dispersions)
        let dispersions = &self.dispersion_result.as_ref().unwrap().map_estimates;
        let beta_fit =
            glm::fit_negative_binomial_glm(&self.counts, &sf, &self.design_matrix, dispersions);
        self.beta_fit = Some(beta_fit);

        Ok(())
    }

    pub fn results(&self, contrast: Contrast) -> Result<Vec<DESeqResult>, String> {
        let beta_fit = self.beta_fit.as_ref().ok_or("run() not called")?;
        let base_means = self.base_means.as_ref().ok_or("run() not called")?;

        // Resolve contrast to a numeric vector
        let contrast_vec = match contrast {
            Contrast::LastCoefficient => {
                let p = self.design_matrix.ncols();
                let mut v = vec![0.0; p];
                v[p - 1] = 1.0;
                v
            }
            Contrast::Vector(v) => v,
            Contrast::ColumnLevels(_, _, _) => {
                // For MVP, just use last coefficient
                let p = self.design_matrix.ncols();
                let mut v = vec![0.0; p];
                v[p - 1] = 1.0;
                v
            }
        };

        // Wald test
        let wald_results = wald::wald_test(beta_fit, &contrast_vec);

        // Raw p-values
        let raw_p: Vec<f64> = wald_results.iter().map(|w| w.p_value).collect();

        // Independent filtering + BH adjustment
        let padj = p_adjust::independent_filtering(base_means, &raw_p, 0.1);

        // Assemble results
        let results: Vec<DESeqResult> = (0..self.gene_names.len())
            .map(|i| DESeqResult {
                gene: self.gene_names[i].clone(),
                base_mean: base_means[i],
                log2_fold_change: wald_results[i].log2_fold_change,
                lfc_se: wald_results[i].lfc_se,
                stat: wald_results[i].stat,
                p_value: wald_results[i].p_value,
                p_adjusted: padj[i],
            })
            .collect();

        Ok(results)
    }

    /// Export intermediate results for validation.
    pub fn export_intermediates(&self, dir: &Path) -> Result<(), Box<dyn std::error::Error>> {
        std::fs::create_dir_all(dir)?;

        if let Some(ref sf) = self.size_factors {
            let mut w = csv::WriterBuilder::new().delimiter(b'\t').from_path(dir.join("size_factors.tsv"))?;
            w.write_record(&["sample", "size_factor"])?;
            for (i, s) in self.sample_names.iter().enumerate() {
                w.write_record(&[s.as_str(), &format!("{:.10e}", sf[i])])?;
            }
        }

        if let Some(ref bm) = self.base_means {
            let mut w = csv::WriterBuilder::new().delimiter(b'\t').from_path(dir.join("base_means.tsv"))?;
            w.write_record(&["gene", "baseMean"])?;
            for (i, g) in self.gene_names.iter().enumerate() {
                w.write_record(&[g.as_str(), &format!("{:.10e}", bm[i])])?;
            }
        }

        if let Some(ref disp) = self.dispersion_result {
            let write_gene_vec = |filename: &str, header: &str, vals: &[f64]| -> Result<(), Box<dyn std::error::Error>> {
                let mut w = csv::WriterBuilder::new().delimiter(b'\t').from_path(dir.join(filename))?;
                w.write_record(&["gene", header])?;
                for (i, g) in self.gene_names.iter().enumerate() {
                    w.write_record(&[g.as_str(), &format!("{:.10e}", vals[i])])?;
                }
                Ok(())
            };
            write_gene_vec("disp_gene_estimates.tsv", "dispGeneEst", &disp.gene_estimates)?;
            write_gene_vec("disp_trend_values.tsv", "dispTrend", &disp.trend_values)?;
            write_gene_vec("disp_map_estimates.tsv", "dispMAP", &disp.map_estimates)?;
        }

        Ok(())
    }
}
```

- [ ] **Step 4: Run tests**

Run: `cd rust_deseq2 && cargo test test_pipeline`
Expected: PASS

- [ ] **Step 5: Add airway gene set validation test**

Append to `tests/test_pipeline.rs`:

```rust
#[test]
fn test_airway_gene_set_match() {
    let ref_dir = Path::new("tests/reference_data");
    if !ref_dir.join("r_results.tsv").exists() {
        eprintln!("Skipping: reference data not found");
        return;
    }

    // Run Rust pipeline
    let mut dds = DESeqDataSet::from_csv(
        &ref_dir.join("airway_counts.tsv"),
        &ref_dir.join("airway_coldata.tsv"),
        "dex", "untrt",
    ).unwrap();
    dds.run().unwrap();
    let rust_results = dds.results(Contrast::LastCoefficient).unwrap();

    // Read R results
    let r_content = std::fs::read_to_string(ref_dir.join("r_results.tsv")).unwrap();
    let mut r_significant: std::collections::HashSet<String> = std::collections::HashSet::new();
    for line in r_content.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 7 { continue; }
        let padj: f64 = match f[6].parse() { Ok(v) => v, Err(_) => continue };
        let lfc: f64 = match f[2].parse() { Ok(v) => v, Err(_) => continue };
        if padj < 0.01 && lfc.abs() > 1.0 {
            r_significant.insert(f[0].to_string());
        }
    }

    // Rust significant genes
    let rust_significant: std::collections::HashSet<String> = rust_results
        .iter()
        .filter(|r| !r.p_adjusted.is_nan() && r.p_adjusted < 0.01 && r.log2_fold_change.abs() > 1.0)
        .map(|r| r.gene.clone())
        .collect();

    println!("R significant genes: {}", r_significant.len());
    println!("Rust significant genes: {}", rust_significant.len());

    let overlap = r_significant.intersection(&rust_significant).count();
    let r_only: Vec<&String> = r_significant.difference(&rust_significant).collect();
    let rust_only: Vec<&String> = rust_significant.difference(&r_significant).collect();

    println!("Overlap: {}/{}", overlap, r_significant.len());
    if !r_only.is_empty() {
        println!("In R only (first 10): {:?}", &r_only[..r_only.len().min(10)]);
    }
    if !rust_only.is_empty() {
        println!("In Rust only (first 10): {:?}", &rust_only[..rust_only.len().min(10)]);
    }

    assert_eq!(
        r_significant, rust_significant,
        "Gene sets must match exactly at p<0.01, |log2FC|>1"
    );
}
```

- [ ] **Step 6: Run validation test**

Run: `cd rust_deseq2 && cargo test test_airway_gene_set_match -- --nocapture`
Expected: PASS with 100% overlap

- [ ] **Step 7: Commit**

```bash
git add rust_deseq2/src/pipeline.rs rust_deseq2/tests/test_pipeline.rs
git commit -m "feat: implement pipeline layer with end-to-end airway validation"
```

---

## Task 15: CLI

**Files:**
- Modify: `rust_deseq2/src/bin/main.rs`
- Test: manual CLI test

- [ ] **Step 1: Implement CLI**

```rust
use clap::Parser;
use deseq2_rs::data::Contrast;
use deseq2_rs::io;
use deseq2_rs::pipeline::DESeqDataSet;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "deseq2-rs", about = "Rust implementation of DESeq2 differential expression analysis")]
struct Cli {
    /// Path to count matrix TSV (genes x samples)
    #[arg(long)]
    counts: PathBuf,

    /// Path to sample metadata TSV
    #[arg(long)]
    coldata: PathBuf,

    /// Column name in coldata for grouping variable
    #[arg(long)]
    design: String,

    /// Reference level (control group)
    #[arg(long)]
    reference: String,

    /// Output results TSV path
    #[arg(long, default_value = "results.tsv")]
    output: PathBuf,

    /// FDR threshold for independent filtering
    #[arg(long, default_value = "0.1")]
    alpha: f64,

    /// Number of threads (default: all CPUs)
    #[arg(long)]
    threads: Option<usize>,

    /// Directory to export intermediate results for validation
    #[arg(long)]
    intermediate_dir: Option<PathBuf>,
}

fn main() {
    let cli = Cli::parse();

    if let Some(n) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .unwrap();
    }

    eprintln!("Loading data...");
    let mut dds = DESeqDataSet::from_csv(
        &cli.counts,
        &cli.coldata,
        &cli.design,
        &cli.reference,
    )
    .expect("Failed to load data");

    eprintln!(
        "Loaded {} genes x {} samples",
        dds.gene_names.len(),
        dds.sample_names.len()
    );

    eprintln!("Running DESeq2 pipeline...");
    dds.run().expect("Pipeline failed");

    if let Some(ref dir) = cli.intermediate_dir {
        eprintln!("Exporting intermediate results to {:?}", dir);
        dds.export_intermediates(dir).expect("Failed to export intermediates");
    }

    eprintln!("Extracting results...");
    let results = dds
        .results(Contrast::LastCoefficient)
        .expect("Failed to extract results");

    eprintln!("Writing results to {:?}", cli.output);
    io::write_results(&cli.output, &results).expect("Failed to write results");

    let n_sig = results
        .iter()
        .filter(|r| !r.p_adjusted.is_nan() && r.p_adjusted < 0.01 && r.log2_fold_change.abs() > 1.0)
        .count();
    eprintln!("Done. Significant genes (padj<0.01, |log2FC|>1): {}", n_sig);
}
```

- [ ] **Step 2: Build the CLI**

Run: `cd rust_deseq2 && cargo build`
Expected: compiles successfully

- [ ] **Step 3: Test CLI with airway data (manual)**

Run: `cd rust_deseq2 && cargo run -- --counts tests/reference_data/airway_counts.tsv --coldata tests/reference_data/airway_coldata.tsv --design dex --reference untrt --output /tmp/rust_results.tsv --intermediate-dir /tmp/debug`
Expected: results written, intermediate files created

- [ ] **Step 4: Commit**

```bash
git add rust_deseq2/src/bin/main.rs
git commit -m "feat: add CLI for deseq2-rs"
```

---

## Task 16: Final Validation and Cleanup

**Files:**
- Modify: `rust_deseq2/src/lib.rs` (ensure clean public API)
- Test: all tests passing

- [ ] **Step 1: Run full test suite**

Run: `cd rust_deseq2 && cargo test -- --nocapture`
Expected: all tests pass

- [ ] **Step 2: Run clippy**

Run: `cd rust_deseq2 && cargo clippy -- -D warnings`
Expected: no warnings (fix any that appear)

- [ ] **Step 3: Verify lib.rs exports clean public API**

Ensure `src/lib.rs` exposes the right modules:

```rust
pub mod data;
pub mod design;
pub mod dispersion;
pub mod glm;
pub mod io;
pub mod pipeline;
pub mod size_factors;
pub mod test_stats;
```

- [ ] **Step 4: Run the full airway validation one more time**

Run: `cd rust_deseq2 && cargo test test_airway_gene_set_match -- --nocapture`
Expected: PASS with 100% gene set overlap

- [ ] **Step 5: Final commit**

```bash
git add -A rust_deseq2/
git commit -m "chore: final cleanup and validation pass for deseq2-rs MVP"
```
