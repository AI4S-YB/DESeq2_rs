pub mod gene_estimate;
pub mod map_estimate;
pub mod trend_fit;

use crate::data::DispersionResult;
use crate::size_factors::base_means;
use faer::Mat;

/// Run the full dispersion estimation pipeline.
///
/// Steps:
/// 1. Gene-wise MLE dispersion estimation
/// 2. Compute per-gene base means
/// 3. Fit parametric trend dispersion = a0 + a1/mean
/// 4. Estimate empirical Bayes prior variance
/// 5. MAP shrinkage of dispersions toward the fitted trend
///
/// # Returns
/// `(DispersionResult, mu_matrix)` where `mu_matrix` is the fitted means
/// (genes x samples) from the gene-wise GLM fit.
pub fn estimate_dispersions(
    counts: &Mat<f64>,
    size_factors: &[f64],
    design_matrix: &Mat<f64>,
) -> Result<(DispersionResult, Mat<f64>), String> {
    // Step 1: Gene-wise MLE
    let (gene_estimates, mu) =
        gene_estimate::estimate_dispersions_gene(counts, size_factors, design_matrix);

    // Step 2: Base means for trend fitting
    let bm = base_means(counts, size_factors);

    // Step 3: Parametric trend fit
    let trend_coeffs = trend_fit::fit_dispersion_trend(&bm, &gene_estimates)?;
    let trend_values: Vec<f64> = bm.iter().map(|m| trend_coeffs.eval(*m)).collect();

    // Step 4: Prior variance
    let prior_var = trend_fit::estimate_prior_variance(
        &gene_estimates,
        &trend_values,
        counts.ncols(),
        design_matrix.ncols(),
    );

    // Step 5: MAP shrinkage
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
