use crate::data::*;
use crate::design::build_design_matrix;
use crate::dispersion;
use crate::glm;
use crate::io;
use crate::size_factors;
use crate::test_stats;
use faer::Mat;
use std::error::Error;
use std::path::Path;

pub struct DESeqDataSet {
    pub counts: Mat<f64>,
    pub col_data: DataFrame,
    pub design_matrix: Mat<f64>,
    pub gene_names: Vec<String>,
    pub sample_names: Vec<String>,
    size_factors: Option<Vec<f64>>,
    base_means: Option<Vec<f64>>,
    dispersion_result: Option<DispersionResult>,
    beta_fit: Option<BetaFitResult>,
    mu: Option<Mat<f64>>,
}

impl DESeqDataSet {
    /// Construct a DESeqDataSet from CSV/TSV files.
    ///
    /// Reads the count matrix and column data, then builds the design matrix.
    /// If sample order in counts differs from coldata, count columns are reordered
    /// to match coldata row order.
    pub fn from_csv(
        count_path: &Path,
        coldata_path: &Path,
        design_column: &str,
        reference_level: &str,
    ) -> Result<Self, Box<dyn Error>> {
        let (counts, gene_names, count_sample_names) = io::read_count_matrix(count_path)?;
        let (col_data, coldata_sample_names) = io::read_coldata(coldata_path)?;

        // Reorder counts columns to match coldata row order if needed
        let (counts, sample_names) = if count_sample_names == coldata_sample_names {
            (counts, coldata_sample_names)
        } else {
            // Build a mapping: for each coldata sample, find its column index in counts
            let mut col_order = Vec::with_capacity(coldata_sample_names.len());
            for name in &coldata_sample_names {
                let idx = count_sample_names
                    .iter()
                    .position(|s| s == name)
                    .ok_or_else(|| {
                        format!(
                            "sample '{}' in coldata not found in count matrix columns",
                            name
                        )
                    })?;
                col_order.push(idx);
            }
            let n_genes = counts.nrows();
            let n_samples = col_order.len();
            let reordered = Mat::from_fn(n_genes, n_samples, |i, j| counts[(i, col_order[j])]);
            (reordered, coldata_sample_names)
        };

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

    /// Construct a DESeqDataSet from pre-built components.
    pub fn from_matrix(
        counts: Mat<f64>,
        col_data: DataFrame,
        design_matrix: Mat<f64>,
        gene_names: Vec<String>,
        sample_names: Vec<String>,
    ) -> Result<Self, String> {
        if counts.ncols() != design_matrix.nrows() {
            return Err(format!(
                "counts has {} columns but design matrix has {} rows",
                counts.ncols(),
                design_matrix.nrows()
            ));
        }
        if counts.nrows() != gene_names.len() {
            return Err(format!(
                "counts has {} rows but gene_names has {} entries",
                counts.nrows(),
                gene_names.len()
            ));
        }
        if counts.ncols() != sample_names.len() {
            return Err(format!(
                "counts has {} columns but sample_names has {} entries",
                counts.ncols(),
                sample_names.len()
            ));
        }
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

    /// Run the full DESeq2 pipeline: size factors, dispersions, GLM fitting.
    pub fn run(&mut self) -> Result<(), Box<dyn Error>> {
        // Step 1: Estimate size factors
        let sf = size_factors::estimate_size_factors(&self.counts);

        // Step 2: Compute base means
        let bm = size_factors::base_means(&self.counts, &sf);

        // Step 3: Estimate dispersions (gene-wise, trend, MAP) and fitted means
        let (disp_result, mu) =
            dispersion::estimate_dispersions(&self.counts, &sf, &self.design_matrix)?;

        // Step 4: Fit negative binomial GLM with MAP dispersions
        let beta_fit = glm::fit_negative_binomial_glm(
            &self.counts,
            &sf,
            &self.design_matrix,
            &disp_result.map_estimates,
        );

        self.size_factors = Some(sf);
        self.base_means = Some(bm);
        self.dispersion_result = Some(disp_result);
        self.beta_fit = Some(beta_fit);
        self.mu = Some(mu);

        Ok(())
    }

    /// Extract results for a given contrast.
    pub fn results(&self, contrast: Contrast) -> Result<Vec<DESeqResult>, String> {
        let beta_fit = self
            .beta_fit
            .as_ref()
            .ok_or("pipeline has not been run yet")?;
        let base_means = self
            .base_means
            .as_ref()
            .ok_or("pipeline has not been run yet")?;

        let n_params = self.design_matrix.ncols();

        // Resolve contrast to a numeric vector
        let contrast_vec = match contrast {
            Contrast::LastCoefficient => {
                let mut v = vec![0.0; n_params];
                v[n_params - 1] = 1.0;
                v
            }
            Contrast::Vector(v) => {
                if v.len() != n_params {
                    return Err(format!(
                        "contrast vector length {} does not match {} parameters",
                        v.len(),
                        n_params
                    ));
                }
                v
            }
            Contrast::ColumnLevels(_column, _level, _reference) => {
                // For a two-level factor with intercept + treatment coding,
                // the last coefficient is the contrast between the two levels.
                // A more general implementation would look up which column
                // corresponds to the level, but for now use LastCoefficient logic.
                let mut v = vec![0.0; n_params];
                v[n_params - 1] = 1.0;
                v
            }
        };

        // Wald test
        let wald_results = test_stats::wald_test(beta_fit, &contrast_vec);

        // Collect raw p-values
        let raw_p: Vec<f64> = wald_results.iter().map(|r| r.p_value).collect();

        // Independent filtering + BH adjustment
        let padj = test_stats::independent_filtering(base_means, &raw_p, 0.1);

        // Assemble final results
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

    /// Return size-factor normalized counts (`counts(dds, normalized=TRUE)` in DESeq2).
    pub fn normalized_counts(&self) -> Result<Mat<f64>, String> {
        let size_factors = self
            .size_factors
            .as_ref()
            .ok_or("pipeline has not been run yet")?;
        Ok(size_factors::normalized_counts(&self.counts, size_factors))
    }

    /// Write size-factor normalized counts to a TSV file.
    pub fn write_normalized_counts(&self, path: &Path) -> Result<(), Box<dyn Error>> {
        let normalized = self.normalized_counts().map_err(std::io::Error::other)?;
        io::write_named_matrix(
            path,
            &normalized,
            &self.gene_names,
            &self.sample_names,
            "gene",
        )
    }

    /// Export intermediate results (size factors, dispersions, etc.) to a directory.
    pub fn export_intermediates(&self, dir: &Path) -> Result<(), Box<dyn Error>> {
        std::fs::create_dir_all(dir)?;

        if let Some(ref sf) = self.size_factors {
            let path = dir.join("size_factors.tsv");
            let mut wtr = csv::WriterBuilder::new()
                .delimiter(b'\t')
                .from_path(&path)?;
            wtr.write_record(&["sample", "size_factor"])?;
            for (name, val) in self.sample_names.iter().zip(sf) {
                wtr.write_record(&[name.as_str(), &format!("{:.10e}", val)])?;
            }
        }

        if let Some(ref bm) = self.base_means {
            let path = dir.join("base_means.tsv");
            let mut wtr = csv::WriterBuilder::new()
                .delimiter(b'\t')
                .from_path(&path)?;
            wtr.write_record(&["gene", "baseMean"])?;
            for (name, val) in self.gene_names.iter().zip(bm) {
                wtr.write_record(&[name.as_str(), &format!("{:.10e}", val)])?;
            }
        }

        if let Some(ref disp) = self.dispersion_result {
            // Gene estimates
            let path = dir.join("disp_gene_estimates.tsv");
            let mut wtr = csv::WriterBuilder::new()
                .delimiter(b'\t')
                .from_path(&path)?;
            wtr.write_record(&["gene", "dispGeneEst"])?;
            for (name, val) in self.gene_names.iter().zip(&disp.gene_estimates) {
                wtr.write_record(&[name.as_str(), &format!("{:.10e}", val)])?;
            }

            // Trend values
            let path = dir.join("disp_trend_values.tsv");
            let mut wtr = csv::WriterBuilder::new()
                .delimiter(b'\t')
                .from_path(&path)?;
            wtr.write_record(&["gene", "dispTrend"])?;
            for (name, val) in self.gene_names.iter().zip(&disp.trend_values) {
                wtr.write_record(&[name.as_str(), &format!("{:.10e}", val)])?;
            }

            // MAP estimates
            let path = dir.join("disp_map_estimates.tsv");
            let mut wtr = csv::WriterBuilder::new()
                .delimiter(b'\t')
                .from_path(&path)?;
            wtr.write_record(&["gene", "dispMAP"])?;
            for (name, val) in self.gene_names.iter().zip(&disp.map_estimates) {
                wtr.write_record(&[name.as_str(), &format!("{:.10e}", val)])?;
            }
        }

        Ok(())
    }
}
