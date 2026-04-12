# DESeq2-rs TODO

Current state: MVP complete, 672/675 (99.6%) gene set overlap with R DESeq2 on airway dataset at p<0.01 & |log2FC|>1.

---

## Accuracy Improvements (Closing the 0.4% Gap)

### P0: Use NB GLM-fitted mu for dispersion estimation

R's `estimateDispersionsGeneEst` uses `fitNbinomGLMs` (not linear model) to compute mu when n_samples > 3. Our code uses the linear model projection, which causes ~5% error in the dispersion trend coefficient a1 (3.22 vs R's 3.38). This cascades through MAP shrinkage to p-values, producing 39 extra significant genes.

**Fix:** In `estimate_dispersions_gene`, after computing rough dispersions, run IRLS to get mu (reuse `glm::irls` logic), then use that mu for the `fitDisp` line search. A `compute_mu_from_glm` skeleton already exists in `gene_estimate.rs` (marked `#[allow(dead_code)]`).

### P0: Grid search fallback for non-converged genes

R uses `fitDispGrid` for genes where the line search didn't converge. A `fit_dispersion_grid` function is already implemented but disabled because the grid parameters didn't match R exactly, causing worse overlap. Needs tuning to match R's exact 15-point grid from `log(minDisp/10)` to `log(max(10, ncol))`.

### P1: Prior variance estimation refinement

Current prior variance is 15% higher than R's (1.51 vs 1.32), due to ~500 gene-wise outliers inflating the residual variance. Fixing the two items above should resolve this automatically.

### P1: Match R's dispersion-to-final assignment

R has additional logic when mapping `dispMAP` to `dispersion`:
- Outlier genes (high Cook's distance) may get their dispersion replaced
- Genes below `minDisp` threshold get special handling
- 329 genes differ between R's `dispMAP` and `dispersion` columns

---

## Feature Additions

### LRT (Likelihood Ratio Test)

Add `nbinomLRT` as an alternative to the Wald test. Requires fitting a reduced model and computing the likelihood ratio statistic. Useful for multi-factor designs.

### lfcShrink

Post-hoc log2 fold change shrinkage methods:
- `normal` (simplest, uses a normal prior on LFC)
- `apeglm` (adaptive prior, requires additional optimization)
- `ashr` (adaptive shrinkage, external dependency)

### Formula parsing

Support R-style formulas (`~ condition + batch`) instead of requiring `--design <column>` + `--reference <level>`. Would need a simple formula parser to generate multi-factor design matrices.

### Multi-factor designs

Currently only supports single-factor designs (intercept + one categorical variable). Extend to:
- Multiple categorical variables (`~ condition + batch`)
- Interaction terms (`~ condition * genotype`)
- Continuous covariates

### VST / rlog transformations

Variance-stabilizing transformation and regularized log transformation for visualization and clustering. These are commonly used downstream of DESeq2.

### Outlier handling

Cook's distance-based outlier detection and replacement, matching R's `cooksCutoff` logic in `results()`.

### tximport-style input

Support transcript-level quantification inputs (salmon, kallisto) with gene-level summarization.

---

## Engineering Improvements

### Performance

- **Benchmark suite**: Add criterion benchmarks for the hot path (dispersion estimation, IRLS)
- **SIMD**: Profile and optimize inner loops in `log_posterior` and IRLS with explicit SIMD
- **Memory**: Current implementation creates many temporary Vec allocations in the per-gene loops; consider pre-allocating per-thread buffers

### Testing

- **Integration tests against R**: Add `tests/test_pipeline.rs` that loads reference data and asserts per-step accuracy (test structure is in the plan, data files exist)
- **Property-based tests**: Use proptest for NB log-posterior derivatives (numerical vs analytical)
- **Edge cases**: Genes with single non-zero sample, extreme dispersion values, rank-deficient designs

### API

- **Error types**: Replace `Box<dyn Error>` and `String` errors with a proper error enum
- **Builder pattern**: `DESeqDataSet::builder()` for more ergonomic construction
- **Streaming I/O**: Support reading count matrices that don't fit in memory (chunked processing)
- **Serde**: Derive Serialize/Deserialize for result types to support JSON output

### CI/CD

- **GitHub Actions**: cargo test, clippy, rustfmt on push
- **R reference data**: Cache the airway reference data as a CI artifact
- **Cross-compilation**: Verify builds on Linux, macOS, Windows

### Desktop UI Integration

The crate is designed as a library for use in a cross-platform desktop UI (e.g., Tauri). Future work:
- Progress callbacks for long-running operations (`run()` with a progress reporter)
- Cancellation support (check a flag in the per-gene parallel loops)
- Async wrapper for non-blocking UI integration
