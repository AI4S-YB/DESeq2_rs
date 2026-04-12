// Parametric dispersion trend fitting for DESeq2
//
// Implements fit_dispersion_trend and estimate_prior_variance.
// Matches R's parametricDispersionFit() using iterated Gamma GLM with identity link.

use crate::data::DispTrendCoeffs;
use crate::dispersion::gene_estimate::trigamma;

/// Fit parametric dispersion trend: dispersion = a0 + a1/mean
///
/// Uses iterated Gamma GLM with identity link (IRLS) matching R's
/// `parametricDispersionFit()`. Outlier genes are excluded each iteration
/// via a residual filter.
///
/// # Arguments
/// * `means` - per-gene base means
/// * `disps` - per-gene dispersion estimates (gene-wise MLEs)
///
/// # Returns
/// `DispTrendCoeffs { asympt_disp: a0, extra_pois: a1 }` or an error string.
pub fn fit_dispersion_trend(
    means: &[f64],
    disps: &[f64],
) -> Result<DispTrendCoeffs, String> {
    let n = means.len();
    if n != disps.len() {
        return Err(format!(
            "means ({}) and disps ({}) must have the same length",
            n,
            disps.len()
        ));
    }
    if n < 3 {
        return Err("Need at least 3 data points to fit dispersion trend".to_string());
    }

    // Initial coefficients: a0=0.1 (asymptotic dispersion), a1=1.0 (extra Poisson)
    let mut a0: f64 = 0.1;
    let mut a1: f64 = 1.0;

    let max_outer = 10;

    for _outer in 0..max_outer {
        // Step a: compute residuals r_i = disp_i / (a0 + a1/mean_i)
        // Filter to genes where 1e-4 < r_i < 15
        let mut good_means: Vec<f64> = Vec::new();
        let mut good_disps: Vec<f64> = Vec::new();

        for i in 0..n {
            let m = means[i];
            if m <= 0.0 || !disps[i].is_finite() || disps[i] <= 0.0 {
                continue;
            }
            let fitted = a0 + a1 / m;
            if fitted <= 0.0 {
                continue;
            }
            let r = disps[i] / fitted;
            if r > 1e-4 && r < 15.0 {
                good_means.push(m);
                good_disps.push(disps[i]);
            }
        }

        if good_means.len() < 3 {
            return Err("Too few genes pass the outlier filter for trend fitting".to_string());
        }

        // Step c: Fit Gamma GLM with identity link using IRLS
        // Model: E[Y] = mu_i = a0 * 1 + a1 * (1/mean_i)
        // Design matrix X: [1, 1/mean_i]
        // Variance function: V(mu) = mu^2 (Gamma)
        // Weight: w_i = 1 / mu_i^2
        // Working response: z_i = y_i  (identity link, so z = y)
        //
        // WLS normal equations: (X^T W X) beta = X^T W z
        // where W = diag(w_i)

        let ng = good_means.len();
        // Build design matrix rows: [1, 1/mean_i]
        let x0: Vec<f64> = vec![1.0; ng];
        let x1: Vec<f64> = good_means.iter().map(|&m| 1.0 / m).collect();

        let new_coeffs = irls_gamma_identity(&x0, &x1, &good_disps, a0, a1, 50)?;

        // Step d: check convergence
        let new_a0 = new_coeffs.0.max(1e-10);
        let new_a1 = new_coeffs.1.max(1e-10);

        let change = (new_a0.ln() - a0.ln()).powi(2) + (new_a1.ln() - a1.ln()).powi(2);

        a0 = new_a0;
        a1 = new_a1;

        if change < 1e-6 {
            break;
        }
    }

    Ok(DispTrendCoeffs {
        asympt_disp: a0,
        extra_pois: a1,
    })
}

/// IRLS for Gamma GLM with identity link.
///
/// Iterates WLS updates:
///   mu_i = a0 + a1 * x1_i   (identity link)
///   weight_i = 1 / mu_i^2   (Gamma variance)
///   working_response_i = y_i (identity link, z = y)
///
/// Returns (a0, a1).
fn irls_gamma_identity(
    _x0: &[f64],
    x1: &[f64],
    y: &[f64],
    a0_init: f64,
    a1_init: f64,
    max_iter: usize,
) -> Result<(f64, f64), String> {
    let n = y.len();
    let mut a0 = a0_init.max(1e-10);
    let mut a1 = a1_init.max(1e-10);

    for _iter in 0..max_iter {
        // Compute fitted values and weights
        // mu_i = a0 + a1/mean_i = a0*1 + a1*x1_i
        let mut w00 = 0.0_f64;
        let mut w01 = 0.0_f64;
        let mut w11 = 0.0_f64;
        let mut rhs0 = 0.0_f64;
        let mut rhs1 = 0.0_f64;

        let mut all_ok = true;
        for i in 0..n {
            let mu_i = a0 + a1 * x1[i];
            if mu_i <= 0.0 {
                all_ok = false;
                break;
            }
            let w_i = 1.0 / (mu_i * mu_i);
            // WLS: X^T W X * beta = X^T W z, z_i = y_i
            w00 += w_i * 1.0 * 1.0;
            w01 += w_i * 1.0 * x1[i];
            w11 += w_i * x1[i] * x1[i];
            rhs0 += w_i * 1.0 * y[i];
            rhs1 += w_i * x1[i] * y[i];
        }

        if !all_ok {
            break;
        }

        // Solve 2x2 system: [[w00, w01],[w01, w11]] * [a0, a1] = [rhs0, rhs1]
        let det = w00 * w11 - w01 * w01;
        if det.abs() < 1e-30 {
            return Err("Singular matrix in IRLS for dispersion trend".to_string());
        }

        let new_a0 = (w11 * rhs0 - w01 * rhs1) / det;
        let new_a1 = (w00 * rhs1 - w01 * rhs0) / det;

        // Guard against non-positive coefficients
        let new_a0 = new_a0.max(1e-10);
        let new_a1 = new_a1.max(1e-10);

        let change = (new_a0 - a0).powi(2) + (new_a1 - a1).powi(2);
        a0 = new_a0;
        a1 = new_a1;

        if change < 1e-12 {
            break;
        }
    }

    Ok((a0, a1))
}

/// Estimate empirical Bayes prior variance for MAP dispersion shrinkage.
///
/// Matches DESeq2's `estimateDispersionsPriorVar`:
///
/// - residuals = log(geneDisp) - log(trendDisp)   (only genes with geneDisp >= 1e-6)
/// - observedVar = var(residuals)
/// - expectedVar = trigamma((n_samples - n_params) / 2)
/// - priorVar = max(observedVar - expectedVar, 0.25)
///
/// # Arguments
/// * `gene_disps` - per-gene dispersion MLEs
/// * `trend_disps` - per-gene trend-fitted dispersions
/// * `n_samples` - number of samples
/// * `n_params` - number of model parameters (design matrix columns)
pub fn estimate_prior_variance(
    gene_disps: &[f64],
    trend_disps: &[f64],
    n_samples: usize,
    n_params: usize,
) -> f64 {
    // Compute log residuals for genes with geneDisp >= 1e-6
    let residuals: Vec<f64> = gene_disps
        .iter()
        .zip(trend_disps.iter())
        .filter_map(|(&g, &t)| {
            if g >= 1e-6 && t > 0.0 && g.is_finite() && t.is_finite() {
                Some(g.ln() - t.ln())
            } else {
                None
            }
        })
        .collect();

    if residuals.is_empty() {
        return 0.25;
    }

    // Compute variance of residuals
    let n = residuals.len() as f64;
    let mean = residuals.iter().sum::<f64>() / n;
    let observed_var = if n > 1.0 {
        residuals.iter().map(|&r| (r - mean).powi(2)).sum::<f64>() / (n - 1.0)
    } else {
        0.0
    };

    // Expected variance from trigamma (sampling variation in log-dispersion MLE)
    let df = (n_samples as f64 - n_params as f64) / 2.0;
    let expected_var = if df > 0.0 { trigamma(df) } else { 0.0 };

    let prior_var = (observed_var - expected_var).max(0.25);
    prior_var
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_trend_fit_known_data() {
        let means: Vec<f64> = (1..=100).map(|i| i as f64 * 10.0).collect();
        let disps: Vec<f64> = means.iter().map(|&m| 0.05 + 2.0 / m).collect();
        let coeffs = fit_dispersion_trend(&means, &disps).unwrap();
        assert_relative_eq!(coeffs.asympt_disp, 0.05, epsilon = 1e-3);
        assert_relative_eq!(coeffs.extra_pois, 2.0, epsilon = 1e-2);
    }

    #[test]
    fn test_prior_variance_basic() {
        let gene = vec![0.1, 0.5, 0.01, 1.0, 0.2];
        let trend = vec![0.08, 0.4, 0.012, 0.9, 0.15];
        let pv = estimate_prior_variance(&gene, &trend, 8, 2);
        assert!(pv >= 0.25);
        assert!(pv.is_finite());
    }
}
