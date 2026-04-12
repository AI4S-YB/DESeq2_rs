# DESeq2-rs: Rust Reimplementation of DESeq2 Core Algorithm

## Overview

Rust reimplementation of the DESeq2 differential gene expression analysis pipeline, targeting numerical equivalence with the R/Bioconductor package on the airway dataset. This crate serves as the computation backend for a cross-platform desktop UI application.

## Goals

- **Minimum viable**: match R DESeq2 results on the airway dataset such that genes filtered at `p<0.01, log2FC>1` are identical
- **Dual interface**: library crate for programmatic use + CLI for standalone analysis
- **Cross-platform**: pure Rust dependencies, no system BLAS/LAPACK required
- **Performance**: per-gene parallelism from the start via rayon

## Scope

### In Scope (MVP)

- Size factor estimation (median-of-ratios)
- Dispersion estimation (gene-wise MLE → parametric trend fit → MAP shrinkage)
- Negative binomial GLM fitting (IRLS with QR decomposition)
- Wald test for differential expression
- Benjamini-Hochberg p-value adjustment with independent filtering
- CSV/TSV input (count matrix + sample metadata)
- TSV output (results table)
- Intermediate result export for validation
- Validation report comparing Rust vs R results

### Out of Scope (Future)

- LRT (likelihood ratio test)
- lfcShrink (apeglm, ashr, normal shrinkage)
- Formula parsing (R-style `~ condition + batch`)
- tximport-style transcript quantification input
- VST / rlog transformations
- Plotting
- Outlier replacement (Cook's distance based)

## Architecture: Functional Core + Pipeline Shell

### Design Rationale

- **Functional core**: pure functions with explicit inputs/outputs, no side effects. Enables independent testing, easy R-result comparison at every step, and natural rayon parallelism.
- **Pipeline shell**: `DESeqDataSet` struct that orchestrates the core functions, providing a clean API for both crate users and the CLI.

## Data Structures

### Core Types

Note on types: `Matrix<T>` refers to faer's `Mat<T>` (dense column-major matrix). `DataFrame` is a lightweight custom struct (not polars), holding column name-to-values mappings for sample metadata.

```rust
/// Input: count matrix + sample info + design matrix
pub struct DESeqDataSet {
    counts: Mat<f64>,           // genes x samples, raw counts (stored as f64 for faer compat)
    col_data: DataFrame,        // sample metadata — simple struct: HashMap<String, Vec<String>>
    design_matrix: Mat<f64>,    // model matrix (intercept + indicator columns)
    gene_names: Vec<String>,
    sample_names: Vec<String>,

    // Intermediate results (populated by pipeline steps)
    size_factors: Option<Vec<f64>>,
    normalized_counts: Option<Matrix<f64>>,
    dispersions: Option<DispersionResult>,
    beta_fit: Option<BetaFitResult>,
    results: Option<Vec<DESeqResult>>,
}

/// Full dispersion estimation results
pub struct DispersionResult {
    gene_estimates: Vec<f64>,   // gene-wise MLE
    trend_values: Vec<f64>,     // fitted trend alpha(mu)
    map_estimates: Vec<f64>,    // final MAP shrinkage values
    prior_var: f64,             // empirical Bayes prior variance
}

/// Per-gene GLM fitting results
pub struct BetaFitResult {
    coefficients: Matrix<f64>,    // genes x coefficients (log2 scale)
    standard_errors: Matrix<f64>, // genes x coefficients
    deviances: Vec<f64>,
    converged: Vec<bool>,
}

/// Final per-gene result
pub struct DESeqResult {
    gene: String,
    base_mean: f64,
    log2_fold_change: f64,
    lfc_se: f64,                // standard error of LFC
    stat: f64,                  // Wald statistic
    p_value: f64,
    p_adjusted: f64,            // BH adjusted
}
```

### Design Matrix Generation

No R formula parsing. Instead:
- CLI accepts `--design <column>` and `--reference <level>` to generate a 0/1 model matrix (intercept + level indicators vs reference)
- Crate API also accepts a pre-built custom design matrix

### Contrast Specification

```rust
pub enum Contrast {
    /// Default: last coefficient in the design matrix
    LastCoefficient,
    /// Specify column, numerator level, denominator level
    ColumnLevels(String, String, String),
    /// Custom numeric contrast vector
    Vector(Vec<f64>),
}
```

## Functional Core (Pure Functions)

### Size Factor Estimation

```rust
/// Median-of-ratios method (Anders & Huber 2010)
/// 1. Compute geometric mean per gene (skip zeros)
/// 2. Per sample: count / geometric_mean -> take median
pub fn estimate_size_factors(counts: &Matrix<u32>) -> Vec<f64>;
```

### Dispersion Estimation (Three Steps)

```rust
/// Step 1: Gene-wise dispersion MLE
/// Cox-Reid adjusted profile likelihood optimization
/// Armijo line search on log(alpha) scale
/// rayon parallel per-gene
pub fn estimate_dispersions_gene(
    counts: &Matrix<u32>,
    size_factors: &[f64],
    design_matrix: &Matrix<f64>,
    mu_hat: &Matrix<f64>,
) -> Vec<f64>;

/// Internal: single-gene dispersion optimization (line search)
/// Corresponds to R C++ fitDisp()
fn fit_dispersion_gene(
    counts: &[u32],
    mu: &[f64],
    design_matrix: &Matrix<f64>,
    log_alpha_init: f64,
    prior_mean: f64,
    prior_var: f64,
    use_prior: bool,
) -> DispFitResult;

/// Step 2: Parametric trend fitting
/// Fit alpha(mu) = a1/mu + a0 via nonlinear least squares (Gauss-Newton)
/// Returns (a0, a1) coefficients; trend for a given mean is: a1/mean + a0
pub fn fit_dispersion_trend(
    base_means: &[f64],
    gene_dispersions: &[f64],
) -> DispTrendCoeffs;  // struct { a0: f64, a1: f64 } with fn eval(&self, mean: f64) -> f64

/// Step 3: MAP shrinkage
/// Trend values as prior mean, empirical Bayes prior variance
/// Re-runs fit_dispersion_gene() with use_prior=true
/// rayon parallel per-gene
pub fn estimate_dispersions_map(
    counts: &Matrix<u32>,
    size_factors: &[f64],
    design_matrix: &Matrix<f64>,
    mu_hat: &Matrix<f64>,
    trend_values: &[f64],
    prior_var: f64,
) -> Vec<f64>;
```

### Beta Fitting (IRLS)

```rust
/// Negative binomial GLM coefficient fitting
/// IRLS + QR decomposition, corresponds to R C++ fitBeta()
/// rayon parallel per-gene
pub fn fit_negative_binomial_glm(
    counts: &Matrix<u32>,
    size_factors: &[f64],
    design_matrix: &Matrix<f64>,
    dispersions: &[f64],
) -> BetaFitResult;
```

IRLS algorithm per gene:
1. Initialize mu from normalized counts
2. Compute weights: `w_j = mu_j / (1 + alpha * mu_j)`
3. Compute working response: `z_j = log(mu/nf) + (y - mu) / mu`
4. Weighted least squares via QR: `beta = QR_solve(sqrt(W) * X, sqrt(W) * z)`
5. Update `mu = nf * exp(X * beta)`
6. Convergence: deviance relative change < 1e-6

### Wald Test & P-value Adjustment

```rust
/// Wald test: W = beta / SE(beta), p = 2 * (1 - Phi(|W|))
pub fn wald_test(beta_fit: &BetaFitResult, contrast: &[f64]) -> Vec<WaldTestResult>;

/// Benjamini-Hochberg p-value adjustment
pub fn p_adjust_bh(p_values: &[f64]) -> Vec<f64>;

/// Independent filtering: filter low-expression genes then BH correction
pub fn independent_filtering(
    base_means: &[f64],
    p_values: &[f64],
    alpha: f64,
) -> Vec<f64>;
```

### Key Numerical Details

- Dispersion optimization in **log(alpha)** space for stability
- Beta coefficients on **natural log** scale internally, converted to **log2** for output
- Cox-Reid adjustment: `-0.5 * log(det(X'WX))` reduces dispersion estimation bias
- IRLS convergence: deviance relative change < 1e-6
- Special functions (`lgamma`, `digamma`, `trigamma`) from `statrs` crate

## Pipeline Layer

```rust
impl DESeqDataSet {
    /// Load from files
    pub fn from_csv(
        count_path: &Path,
        coldata_path: &Path,
        design_column: &str,
        reference_level: &str,
    ) -> Result<Self>;

    /// Build from in-memory data (crate API)
    pub fn from_matrix(
        counts: Matrix<u32>,
        col_data: DataFrame,
        design_matrix: Matrix<f64>,
        gene_names: Vec<String>,
        sample_names: Vec<String>,
    ) -> Result<Self>;

    /// Run full pipeline (equivalent to R's DESeq())
    pub fn run(&mut self) -> Result<()> {
        self.estimate_size_factors()?;
        self.estimate_dispersions()?;
        self.fit_glm()?;
        self.wald_test()?;
        Ok(())
    }

    /// Extract results (equivalent to R's results())
    pub fn results(&self, contrast: Contrast) -> Result<Vec<DESeqResult>>;
}
```

## CLI

```
deseq2-rs \
  --counts airway_counts.tsv \
  --coldata airway_coldata.tsv \
  --design condition \
  --reference untreated \
  --output results.tsv
```

Output TSV columns: `gene`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`

Optional flags:
- `--alpha 0.1` — FDR threshold for independent filtering
- `--threads N` — rayon thread count (default: all CPUs)
- `--intermediate-dir ./debug/` — export intermediate results for validation

### Intermediate Output (--intermediate-dir)

```
debug/
├── size_factors.tsv
├── base_means.tsv
├── disp_gene_estimates.tsv
├── disp_trend_values.tsv
├── disp_map_estimates.tsv
├── beta_coefficients.tsv
├── wald_statistics.tsv
├── p_values_raw.tsv
└── p_values_adjusted.tsv
```

## Project Structure

```
rust_deseq2/
├── Cargo.toml
├── src/
│   ├── lib.rs                    # Public API re-exports
│   ├── data.rs                   # DESeqDataSet, DataFrame, result types
│   ├── io.rs                     # CSV/TSV reading and writing
│   ├── design.rs                 # Design matrix generation
│   ├── size_factors.rs           # Median-of-ratios
│   ├── dispersion/
│   │   ├── mod.rs
│   │   ├── gene_estimate.rs      # Gene-wise MLE (fitDisp line search)
│   │   ├── trend_fit.rs          # Parametric trend fitting (a1/mu + a0)
│   │   └── map_estimate.rs       # MAP shrinkage
│   ├── glm/
│   │   ├── mod.rs
│   │   └── irls.rs               # IRLS + QR (fitBeta)
│   ├── test_stats/
│   │   ├── mod.rs
│   │   ├── wald.rs               # Wald test
│   │   └── p_adjust.rs           # BH correction + independent filtering
│   ├── pipeline.rs               # DESeqDataSet::run() orchestration
│   └── special_functions.rs      # lgamma, digamma, trigamma wrappers
├── src/bin/
│   └── main.rs                   # CLI entry point (clap)
├── tests/
│   ├── test_size_factors.rs
│   ├── test_dispersion.rs
│   ├── test_glm.rs
│   ├── test_wald.rs
│   ├── test_pipeline.rs          # Airway end-to-end validation
│   └── reference_data/           # R-exported intermediate results
│       ├── airway_counts.tsv
│       ├── airway_coldata.tsv
│       ├── r_size_factors.tsv
│       ├── r_dispersions.tsv
│       ├── r_beta.tsv
│       └── r_results.tsv
├── scripts/
│   └── export_r_reference.R      # Export airway reference data from R DESeq2
```

## Dependencies

```toml
[dependencies]
faer = "0.20"               # Linear algebra (QR, matrix ops) - pure Rust
rayon = "1"                  # Data parallelism (per-gene)
clap = { version = "4", features = ["derive"] }  # CLI
csv = "1"                    # TSV/CSV I/O
statrs = "0.18"              # Normal CDF, digamma, trigamma

[dev-dependencies]
approx = "0.5"               # Float approximate comparison
```

## Validation Strategy

### Validation Flow

1. **R side**: `scripts/export_r_reference.R` runs DESeq2 on airway dataset, exports counts, coldata, and all intermediate results as TSV
2. **Rust side**: integration tests load the same counts + coldata, run the pipeline step by step, compare each step against R reference data

### Per-Step Accuracy Targets

| Step | Comparison Data | Tolerance |
|------|----------------|-----------|
| Size factors | Per-sample value | Relative error < 1e-10 (deterministic) |
| Base means | Per-gene value | Relative error < 1e-10 |
| Gene-wise dispersion | Per-gene value | Relative error < 1e-4 (iterative) |
| Trend fit parameters | a0, a1 | Relative error < 1e-3 (nonlinear fit) |
| MAP dispersion | Per-gene value | Relative error < 1e-4 |
| Beta coefficients | Genes x coeffs matrix | Relative error < 1e-4 |
| Wald stat & p-value | Per-gene pair | Relative error < 1e-4 |
| Adjusted p-value | Per-gene value | Relative error < 1e-4 |
| **Final gate** | **Gene set at p<0.01 & log2FC>1** | **100% overlap** |

### Validation Report

The CLI with `--intermediate-dir` generates a comparison summary showing per-step max relative error and the final gene set overlap, enabling quick identification of any numerical divergence.
