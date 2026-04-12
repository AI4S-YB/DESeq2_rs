// MAP dispersion shrinkage for DESeq2
//
// Re-runs fit_dispersion with use_prior=true to shrink gene-wise dispersions
// toward the fitted trend using an empirical Bayes normal prior.

use crate::dispersion::gene_estimate::fit_dispersion;
use faer::Mat;
use rayon::prelude::*;

/// Compute MAP (maximum a posteriori) dispersion estimates by shrinking
/// gene-wise MLEs toward the parametric trend.
///
/// For each gene, re-runs `fit_dispersion` with `use_prior=true` using:
/// - `prior_mean = log(trend_values[i])`
/// - `prior_var` = empirical Bayes prior variance from `estimate_prior_variance`
///
/// # Arguments
/// * `counts` - raw count matrix (genes x samples)
/// * `size_factors` - per-sample size factors (unused here; mu is precomputed)
/// * `design_matrix` - model design matrix (samples x params)
/// * `mu` - fitted mean matrix from gene-wise estimation (genes x samples)
/// * `gene_dispersions` - gene-wise MLE dispersions
/// * `trend_values` - trend-fitted dispersion values
/// * `prior_var` - empirical Bayes prior variance
///
/// # Returns
/// Vector of MAP dispersion estimates, one per gene.
pub fn estimate_dispersions_map(
    counts: &Mat<f64>,
    _size_factors: &[f64],
    design_matrix: &Mat<f64>,
    mu: &Mat<f64>,
    gene_dispersions: &[f64],
    trend_values: &[f64],
    prior_var: f64,
) -> Vec<f64> {
    let n_genes = counts.nrows();
    let n_samples = counts.ncols();
    let min_disp = 1e-8;
    let max_disp = 10.0_f64.max(n_samples as f64);

    (0..n_genes)
        .into_par_iter()
        .map(|i| {
            let y: Vec<f64> = (0..n_samples).map(|j| counts[(i, j)]).collect();
            let mu_row: Vec<f64> = (0..n_samples).map(|j| mu[(i, j)]).collect();
            let log_alpha_init = gene_dispersions[i].max(min_disp).ln();
            let prior_mean = trend_values[i].max(min_disp).ln();
            let result = fit_dispersion(
                &y,
                &mu_row,
                design_matrix,
                log_alpha_init,
                prior_mean,
                prior_var,
                true,
                true,
            );
            result.log_alpha.exp().clamp(min_disp, max_disp)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use faer::Mat;

    #[test]
    fn test_map_shrinks_toward_trend() {
        let y = vec![5.0, 8.0, 6.0, 10.0];
        let mu = vec![7.0, 7.0, 7.0, 7.0];
        let x = Mat::from_fn(4, 2, |i, j| {
            if j == 0 { 1.0 } else { [0.0, 1.0, 0.0, 1.0][i] }
        });
        // Gene-wise without prior gives some value
        let gene_result = crate::dispersion::gene_estimate::fit_dispersion(
            &y, &mu, &x, 0.0, 0.0, 1.0, false, true,
        );
        // MAP with prior should shift toward prior_mean
        let map_result = crate::dispersion::gene_estimate::fit_dispersion(
            &y, &mu, &x, gene_result.log_alpha, -1.0, 0.5, true, true,
        );
        // MAP log_alpha should be closer to -1.0 (prior) than gene_result
        let gene_dist = (gene_result.log_alpha - (-1.0)).abs();
        let map_dist = (map_result.log_alpha - (-1.0)).abs();
        assert!(
            map_dist <= gene_dist,
            "MAP should shrink toward prior: gene_dist={}, map_dist={}",
            gene_dist,
            map_dist
        );
    }
}
