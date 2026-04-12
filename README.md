# DESeq2-rs

Rust reimplementation of [DESeq2](https://bioconductor.org/packages/DESeq2/) core differential expression analysis pipeline.

## Why This Project Exists

This project has a failed predecessor.

The first attempt burned roughly 200 million tokens trying to do a complete, bit-to-bit accurate reimplementation of DESeq2 in Rust. It used the full Superpowers agentic workflow — brainstorming, specs, plans, subagent-driven development, worktree isolation, oracle fixtures, field-level comparison tools. The architecture was clean, the documentation was thorough, the engineering looked impressive. But it never produced a single correct p-value.

The problem wasn't the tools. The problem was using tools to avoid the hard question: **what's the minimum thing that needs to work first?**

That attempt prioritized "design completeness" over "getting results". It built verification infrastructure before there was anything to verify. It pursued bit-to-bit accuracy before achieving approximate correctness. It kept adding layers of scaffolding around an uncertainty that should have been resolved by writing code and checking output.

This version took a different path.

## What Changed

**One goal, clearly defined:** Match R DESeq2's significant gene set on the airway dataset at `padj < 0.01` and `|log2FC| > 1`. Not bit-to-bit accuracy. Not feature completeness. Just: do these two implementations agree on which genes are differentially expressed?

**Core pipeline first, everything else later:** Size factors → dispersion estimation → NB GLM fitting → Wald test → BH adjustment. No LRT. No lfcShrink. No formula parsing. No VST. No outlier replacement. These are all important features, but none of them matter until the core pipeline produces correct results.

**Validate each step against R, then move on:** Export R's intermediate results (size factors, gene-wise dispersions, trend coefficients, MAP dispersions, beta coefficients, p-values). Compare at each step. Fix divergences. Don't build a comparison framework — just use a Python one-liner.

**Use the same tools, but with judgment:** This version also used Superpowers (brainstorming, writing-plans, subagent-driven-development). The difference was knowing when to follow the tool's suggestions and when to override them. When the plan said "implement grid search for non-converged genes" but the grid search made results worse, we reverted it. When the spec said "pure Rust linear algebra", we didn't fight the faer API for hours — we used simple 2×2 matrix formulas where p is small.

The entire implementation — from design spec to working CLI — was completed in a single session.

## Results

### Accuracy

Validated on the [airway](https://bioconductor.org/packages/airway/) dataset (64,102 genes × 8 samples):

| Metric | Value |
|--------|-------|
| Size factors | Perfect match (relative error < 1e-10) |
| Gene-wise dispersion | Median relative error 6.8e-12 |
| log2FoldChange | Median relative error 0.09% |
| **Gene set overlap (padj<0.01, \|log2FC\|>1)** | **672 / 675 = 99.6%** |

The 3 genes found by R but not Rust are all borderline (R padj = 0.008–0.009, Rust padj = 0.010–0.012). The 39 extra genes found by Rust have systematically lower adjusted p-values, caused by a ~5% difference in the dispersion trend coefficient.

### Speed

| | Wall clock | CPU time | Threads |
|---|---|---|---|
| R DESeq2 | 7.9s | 16.2s | 2 |
| **Rust (this project)** | **0.28s** | **1.35s** | **16** |
| Rust (single-threaded) | 0.94s | 0.94s | 1 |

**~28x faster** (wall clock, multi-threaded) or **~8x faster** (single-threaded) on the same dataset. Per-gene dispersion estimation and GLM fitting are parallelized with rayon.

## Architecture

Functional Core + Pipeline Shell:

- **Functional core**: Pure functions with explicit inputs/outputs. Each algorithm step (size factors, dispersion, GLM, Wald test) is independently testable and comparable against R's intermediate results.
- **Pipeline shell**: `DESeqDataSet` struct orchestrates the core functions, providing a clean API for both library users and the CLI.

```
src/
├── data.rs                  # Core types
├── io.rs                    # TSV I/O
├── design.rs                # Design matrix from coldata
├── size_factors.rs          # Median-of-ratios
├── dispersion/
│   ├── gene_estimate.rs     # NB log-posterior, derivatives, Armijo line search
│   ├── trend_fit.rs         # Parametric trend: α = a₀ + a₁/μ
│   └── map_estimate.rs      # Empirical Bayes MAP shrinkage
├── glm/irls.rs              # IRLS negative binomial GLM (fitBeta)
├── test_stats/
│   ├── wald.rs              # Wald statistic + p-value
│   └── p_adjust.rs          # BH correction + independent filtering
├── pipeline.rs              # DESeqDataSet::run()
└── bin/main.rs              # CLI
```

Dependencies: `faer` (linear algebra, pure Rust), `rayon` (parallelism), `statrs` (statistical functions), `clap` (CLI), `csv` (I/O). No system BLAS/LAPACK required.

## Usage

### CLI

```bash
deseq2-rs \
  --counts counts.tsv \
  --coldata sample_info.tsv \
  --design dex \
  --reference untrt \
  --output results.tsv
```

### As a library

```rust
use deseq2_rs::pipeline::DESeqDataSet;
use deseq2_rs::data::Contrast;

let mut dds = DESeqDataSet::from_csv(
    "counts.tsv", "coldata.tsv", "dex", "untrt"
)?;
dds.run()?;
let results = dds.results(Contrast::LastCoefficient)?;
```

## The Remaining 0.4% Gap

The 39 extra significant genes in Rust come from a specific chain:

1. R uses NB GLM-fitted μ for dispersion estimation; we use linear model μ
2. This causes ~500 gene-wise dispersion outliers (>1% error)
3. Which shifts the trend coefficient a₁ by ~5% (3.22 vs R's 3.38)
4. Which shifts the MAP prior, slightly underestimating dispersions
5. Which produces smaller standard errors → more significant p-values

Fixing this requires implementing R's full mu iteration loop (`fitNbinomGLMs` inside `estimateDispersionsGeneEst`). See [docs/TODO.md](docs/TODO.md) for details.

## Lessons Learned

The first attempt taught an expensive lesson: **tools amplify judgment, they don't replace it.**

A framework that helps you brainstorm, plan, and execute is powerful — if you already know what the most important next step is. If you don't, it will faithfully help you build elaborate scaffolding around your confusion.

The second attempt worked because the goal was specific ("match this gene set"), the path was incremental ("validate each step against R"), and the judgment calls were made by a human who understood the tradeoffs ("99.6% overlap is good enough for MVP, don't chase the last 0.4% today").

No tool gives you abilities you don't have. But if you know what you're doing, the right tools let you do it 28 times faster.

## License

MIT
