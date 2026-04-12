use crate::data::{BetaFitResult, WaldTestResult};
use statrs::distribution::{ContinuousCDF, Normal};

/// Compute Wald test for each gene given a contrast vector.
/// Coefficients in BetaFitResult are on natural log scale.
/// Output log2FoldChange and lfcSE are on log2 scale.
pub fn wald_test(beta_fit: &BetaFitResult, contrast: &[f64]) -> Vec<WaldTestResult> {
    let n_genes = beta_fit.coefficients.nrows();
    let n_params = beta_fit.coefficients.ncols();
    let log2_e = std::f64::consts::E.log2(); // 1/ln(2) ≈ 1.4427
    let normal = Normal::new(0.0, 1.0).unwrap();

    (0..n_genes)
        .map(|i| {
            // c' * beta (natural log scale)
            let estimate: f64 = (0..n_params)
                .map(|k| contrast[k] * beta_fit.coefficients[(i, k)])
                .sum();

            // SE: sqrt(sum(c_k^2 * SE_k^2)) — exact for single-coefficient contrasts
            let se: f64 = {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::BetaFitResult;
    use faer::Mat;
    use approx::assert_relative_eq;

    #[test]
    fn test_wald_test_basic() {
        let coefficients = Mat::from_fn(2, 2, |i, j| {
            [[5.0, 0.5], [3.0, -0.1]][i][j]
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

        let contrast = vec![0.0, 1.0];
        let results = wald_test(&fit, &contrast);
        assert_eq!(results.len(), 2);

        let log2_e = std::f64::consts::E.log2();
        // Gene 0: beta=0.5, SE=0.2 -> stat=2.5
        assert_relative_eq!(results[0].log2_fold_change, 0.5 * log2_e, epsilon = 1e-10);
        assert_relative_eq!(results[0].stat, 2.5, epsilon = 1e-10);
        assert!(results[0].p_value > 0.01 && results[0].p_value < 0.02);
    }
}
