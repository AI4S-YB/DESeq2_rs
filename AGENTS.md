# Repository Guidelines

## Project Structure & Module Organization

This is a Rust 2021 crate named `deseq2-rs`, with a library plus a CLI binary. Core library modules live in `src/`: `data.rs` defines shared types, `io.rs` handles TSV input/output, `design.rs` builds design matrices, `size_factors.rs` implements normalization, `dispersion/` contains dispersion estimation, `glm/` contains IRLS fitting, `test_stats/` contains Wald and p-adjust logic, and `pipeline.rs` orchestrates the full workflow. The CLI entry point is `src/bin/main.rs`.

Tests are split between inline unit tests in `src/**` and integration tests in `tests/`. R-derived fixtures live in `tests/reference_data/`. Documentation and implementation notes are in `docs/`, and `scripts/export_r_reference.R` regenerates reference outputs from R/DESeq2.

## Build, Test, and Development Commands

- `cargo build`: compile the library and `deseq2-rs` binary.
- `cargo test`: run unit and integration tests.
- `cargo test dispersion`: run tests whose names or module paths include `dispersion`.
- `cargo fmt`: format Rust code with rustfmt before submitting changes.
- `cargo clippy -- -D warnings`: run lint checks and fail on warnings.
- `cargo run -- --counts tests/reference_data/airway_counts.tsv --coldata tests/reference_data/airway_coldata.tsv --design dex --reference untrt --output /tmp/results.tsv`: run the CLI against the bundled airway fixture.
- `Rscript scripts/export_r_reference.R`: refresh R reference TSVs when DESeq2 comparison data must change.

## Coding Style & Naming Conventions

Follow idiomatic Rust formatted by rustfmt. Use four-space indentation, `snake_case` for functions, variables, files, and modules, `CamelCase` for types, and `SCREAMING_SNAKE_CASE` for constants. Keep algorithmic code in pure functions where practical, and let `pipeline.rs` handle orchestration. Prefer explicit `Result` returns for recoverable errors; reserve `unwrap`/`expect` for tests and CLI boundary failures.

## Testing Guidelines

Add focused unit tests beside the module being changed, and integration tests under `tests/` for I/O or end-to-end behavior. Name tests descriptively, usually `test_<behavior>` or `test_<step>_vs_r`. For numerical work, use tolerant comparisons with `approx` rather than exact floating-point equality. When changing statistical algorithms, compare against relevant files in `tests/reference_data/` and document any intentional divergence.

## Commit & Pull Request Guidelines

Recent history uses Conventional Commit-style prefixes such as `feat:`, `fix:`, `refactor:`, and `docs:`. Keep commits scoped and imperative, for example `fix: improve dispersion convergence check`.

Pull requests should describe the algorithmic or CLI behavior changed, list validation commands run, and mention any fixture updates. Link related issues or TODO items, and include before/after accuracy or performance notes when touching pipeline, dispersion, GLM, or statistical test code.
