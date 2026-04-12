// IRLS Negative Binomial GLM Fitting (fitBeta)
//
// Implements the iteratively reweighted least squares algorithm for fitting
// negative binomial generalized linear models, matching DESeq2's C++ fitBeta
// function (DESeq2.cpp lines 283-465).

use crate::data::BetaFitResult;
use faer::prelude::*;
use faer::Mat;
use rayon::prelude::*;
use statrs::function::gamma::ln_gamma;

// ---------------------------------------------------------------------------
// Constants matching DESeq2.cpp
// ---------------------------------------------------------------------------

const MINMU: f64 = 0.5;
const LARGE: f64 = 30.0;
const TOL: f64 = 1e-6;
const MAXIT: usize = 100;

// ---------------------------------------------------------------------------
// NB log-PMF
// ---------------------------------------------------------------------------

/// Negative binomial log-PMF parameterized by mean `mu`.
///
/// `log_dnbinom_mu(y, size, mu)` where `size = 1/alpha`.
///
/// Formula: `ln_gamma(y+size) - ln_gamma(size) - ln_gamma(y+1)
///           + y*ln(mu/(mu+size)) + size*ln(size/(mu+size))`
pub fn log_dnbinom_mu(y: f64, size: f64, mu: f64) -> f64 {
    // Handle edge case: if mu is essentially zero
    if mu < 1e-300 {
        if y == 0.0 {
            return 0.0; // P(Y=0 | mu=0) = 1
        } else {
            return f64::NEG_INFINITY;
        }
    }

    ln_gamma(y + size) - ln_gamma(size) - ln_gamma(y + 1.0)
        + y * (mu / (mu + size)).ln()
        + size * (size / (mu + size)).ln()
}

// ---------------------------------------------------------------------------
// Beta initialization
// ---------------------------------------------------------------------------

/// Initialize beta coefficients via normal equations on log-transformed
/// normalized counts.
///
/// Computes `Y = log(counts / sf + 0.1)` then solves `X'X * beta = X'Y`
/// for each gene.
///
/// Returns a genes x params matrix.
pub fn initialize_betas(
    counts: &Mat<f64>,
    size_factors: &[f64],
    design: &Mat<f64>,
) -> Mat<f64> {
    let n_genes = counts.nrows();
    let n_samples = counts.ncols();
    let n_params = design.ncols();

    // Build Y matrix: n_samples x n_genes, where Y[j][i] = log(counts[i][j] / sf[j] + 0.1)
    let y_mat = Mat::<f64>::from_fn(n_samples, n_genes, |j, i| {
        (counts[(i, j)] / size_factors[j] + 0.1).ln()
    });

    // Solve X'X * beta = X'Y via LU decomposition
    // X'X is (n_params x n_params), X'Y is (n_params x n_genes)
    let xtx: Mat<f64> = design.transpose() * design;
    let xty: Mat<f64> = design.transpose() * &y_mat;

    let lu = xtx.partial_piv_lu();
    let beta_t = lu.solve(&xty); // n_params x n_genes

    // Transpose to genes x params
    Mat::from_fn(n_genes, n_params, |i, j| beta_t[(j, i)])
}

// ---------------------------------------------------------------------------
// Small matrix utilities (for p x p systems, typically 2x2)
// ---------------------------------------------------------------------------

/// Solve A * x = b for a small p x p system using LU decomposition.
/// Returns the solution vector as a Vec<f64>.
fn solve_small_system(a: &Mat<f64>, b: &[f64]) -> Vec<f64> {
    let p = a.nrows();
    let b_mat = Mat::from_fn(p, 1, |i, _| b[i]);
    let lu = a.partial_piv_lu();
    let x = lu.solve(&b_mat);
    (0..p).map(|i| x[(i, 0)]).collect()
}

/// Compute the inverse of a small p x p matrix using LU decomposition.
fn invert_small(a: &Mat<f64>) -> Mat<f64> {
    let p = a.nrows();
    let eye = Mat::<f64>::from_fn(p, p, |i, j| if i == j { 1.0 } else { 0.0 });
    let lu = a.partial_piv_lu();
    lu.solve(&eye)
}

// ---------------------------------------------------------------------------
// Per-gene IRLS fitting
// ---------------------------------------------------------------------------

/// Result from fitting a single gene.
struct SingleGeneResult {
    beta: Vec<f64>,
    se: Vec<f64>,
    deviance: f64,
    converged: bool,
}

/// Fit a single gene using IRLS.
///
/// Follows DESeq2.cpp fitBeta non-QR path (lines 386-425) using
/// the normal equations approach: `(X'WX + L) * beta = X'Wz`.
fn fit_single_gene(
    y: &[f64],         // counts for this gene (n_samples)
    sf: &[f64],        // size factors (n_samples)
    x: &Mat<f64>,      // design matrix (n_samples x n_params)
    alpha: f64,        // dispersion for this gene
    beta_init: &[f64], // initial beta (n_params)
    lambda: &[f64],    // ridge penalty (n_params)
) -> SingleGeneResult {
    let n_samples = y.len();
    let n_params = x.ncols();

    let mut beta = beta_init.to_vec();
    let mut mu = vec![0.0; n_samples];
    let size = 1.0 / alpha;

    // Compute initial mu = sf * exp(X * beta), clamped >= MINMU
    for j in 0..n_samples {
        let mut xb = 0.0;
        for k in 0..n_params {
            xb += x[(j, k)] * beta[k];
        }
        mu[j] = (sf[j] * xb.exp()).max(MINMU);
    }

    let mut dev = 0.0;
    let mut dev_old = 0.0;
    let mut converged = false;

    // Ridge penalty matrix (diagonal)
    let ridge = Mat::from_fn(n_params, n_params, |i, j| {
        if i == j { lambda[i] } else { 0.0 }
    });

    for t in 0..MAXIT {
        // Compute weights: w_j = mu_j / (1 + alpha * mu_j)
        let w: Vec<f64> = mu.iter().map(|&m| m / (1.0 + alpha * m)).collect();

        // Working response: z_j = ln(mu_j / sf_j) + (y_j - mu_j) / mu_j
        let z: Vec<f64> = (0..n_samples)
            .map(|j| (mu[j] / sf[j]).ln() + (y[j] - mu[j]) / mu[j])
            .collect();

        // Build X'WX + ridge: (p x p)
        let mut xtwx = Mat::from_fn(n_params, n_params, |r, c| {
            let mut sum = 0.0;
            for j in 0..n_samples {
                sum += x[(j, r)] * w[j] * x[(j, c)];
            }
            sum
        });
        // Add ridge penalty
        for k in 0..n_params {
            xtwx[(k, k)] += ridge[(k, k)];
        }

        // Build X'Wz: (p,)
        let xtwz: Vec<f64> = (0..n_params)
            .map(|k| {
                let mut sum = 0.0;
                for j in 0..n_samples {
                    sum += x[(j, k)] * w[j] * z[j];
                }
                sum
            })
            .collect();

        // Solve (X'WX + L) * beta_new = X'Wz
        beta = solve_small_system(&xtwx, &xtwz);

        // Check for divergence
        if beta.iter().any(|&b| b.abs() > LARGE) {
            converged = false;
            break;
        }

        // Update mu
        for j in 0..n_samples {
            let mut xb = 0.0;
            for k in 0..n_params {
                xb += x[(j, k)] * beta[k];
            }
            mu[j] = (sf[j] * xb.exp()).max(MINMU);
        }

        // Compute deviance
        dev = 0.0;
        for j in 0..n_samples {
            dev += -2.0 * log_dnbinom_mu(y[j], size, mu[j]);
        }

        // Convergence test
        let conv_test = (dev - dev_old).abs() / (dev.abs() + 0.1);
        if conv_test.is_nan() {
            converged = false;
            break;
        }
        if t > 0 && conv_test < TOL {
            converged = true;
            break;
        }
        dev_old = dev;

        // If we reach the last iteration without converging, mark as not converged
        if t == MAXIT - 1 {
            converged = false;
        }
    }

    // Compute sandwich covariance for standard errors
    // sigma = (X'WX + L)^{-1} * X'WX * (X'WX + L)^{-1}
    // Recalculate W with final mu
    let w_final: Vec<f64> = mu.iter().map(|&m| m / (1.0 + alpha * m)).collect();

    let mut xtwx_final = Mat::from_fn(n_params, n_params, |r, c| {
        let mut sum = 0.0;
        for j in 0..n_samples {
            sum += x[(j, r)] * w_final[j] * x[(j, c)];
        }
        sum
    });

    let xtwx_no_ridge = xtwx_final.clone();

    // Add ridge
    for k in 0..n_params {
        xtwx_final[(k, k)] += ridge[(k, k)];
    }

    let xtwxr_inv = invert_small(&xtwx_final);

    // sigma = xtwxr_inv * xtwx_no_ridge * xtwxr_inv
    let sigma: Mat<f64> = &xtwxr_inv * &xtwx_no_ridge * &xtwxr_inv;

    let se: Vec<f64> = (0..n_params)
        .map(|k| {
            let v = sigma[(k, k)];
            if v > 0.0 {
                v.sqrt()
            } else {
                f64::NAN
            }
        })
        .collect();

    SingleGeneResult {
        beta,
        se,
        deviance: dev,
        converged,
    }
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Fit negative binomial GLM for all genes using IRLS.
///
/// This is the Rust equivalent of DESeq2's fitBeta C++ function.
///
/// # Arguments
/// * `counts` - Gene expression count matrix (genes x samples)
/// * `size_factors` - Size factors for each sample
/// * `design_matrix` - Design matrix (samples x params)
/// * `dispersions` - Per-gene dispersion estimates
///
/// # Returns
/// `BetaFitResult` with coefficients and SEs on the natural log scale.
pub fn fit_negative_binomial_glm(
    counts: &Mat<f64>,
    size_factors: &[f64],
    design_matrix: &Mat<f64>,
    dispersions: &[f64],
) -> BetaFitResult {
    let n_genes = counts.nrows();
    let n_params = design_matrix.ncols();

    // Initialize betas
    let beta_init = initialize_betas(counts, size_factors, design_matrix);

    // Ridge penalty: small constant for all coefficients
    let lambda = vec![1e-6; n_params];

    // Fit each gene in parallel
    let results: Vec<SingleGeneResult> = (0..n_genes)
        .into_par_iter()
        .map(|i| {
            // Extract count row for gene i
            let y: Vec<f64> = (0..counts.ncols()).map(|j| counts[(i, j)]).collect();
            let beta_i: Vec<f64> = (0..n_params).map(|k| beta_init[(i, k)]).collect();

            fit_single_gene(&y, size_factors, design_matrix, dispersions[i], &beta_i, &lambda)
        })
        .collect();

    // Assemble results into matrices
    let coefficients = Mat::from_fn(n_genes, n_params, |i, k| results[i].beta[k]);
    let standard_errors = Mat::from_fn(n_genes, n_params, |i, k| results[i].se[k]);
    let deviances: Vec<f64> = results.iter().map(|r| r.deviance).collect();
    let converged: Vec<bool> = results.iter().map(|r| r.converged).collect();

    BetaFitResult {
        coefficients,
        standard_errors,
        deviances,
        converged,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use faer::Mat;

    #[test]
    fn test_log_dnbinom_mu_basic() {
        // Sanity: log-density should be negative for reasonable values
        let ld = log_dnbinom_mu(10.0, 10.0, 10.0);
        assert!(ld < 0.0);
        assert!(ld.is_finite());

        // y=0 should give a valid (finite) log-density
        let ld0 = log_dnbinom_mu(0.0, 10.0, 10.0);
        assert!(ld0.is_finite());
        assert!(ld0 < 0.0);
    }

    #[test]
    fn test_log_dnbinom_mu_matches_r() {
        // Compare with R: dnbinom(5, size=10, mu=5, log=TRUE)
        // Computed via lgamma formula: ~ -1.9458
        let ld = log_dnbinom_mu(5.0, 10.0, 5.0);
        assert_relative_eq!(ld, -1.9458, epsilon = 0.01);
    }

    #[test]
    fn test_initialize_betas() {
        let counts = Mat::from_fn(2, 4, |i, j| {
            [[100.0, 200.0, 150.0, 300.0], [50.0, 40.0, 60.0, 45.0]][i][j]
        });
        let sf = vec![1.0, 1.0, 1.0, 1.0];
        let design = Mat::from_fn(4, 2, |i, j| {
            if j == 0 {
                1.0
            } else {
                [0.0, 1.0, 0.0, 1.0][i]
            }
        });

        let betas = initialize_betas(&counts, &sf, &design);
        assert_eq!(betas.nrows(), 2);
        assert_eq!(betas.ncols(), 2);

        // Gene 0: treated has higher counts => positive beta[1]
        assert!(betas[(0, 1)] > 0.0);
        // Intercept should be roughly log(125) for gene 0 untreated mean
        assert!(betas[(0, 0)] > 4.0 && betas[(0, 0)] < 6.0);
    }

    #[test]
    fn test_irls_simple() {
        let counts = Mat::from_fn(2, 4, |i, j| {
            [[100.0, 200.0, 150.0, 300.0], [50.0, 40.0, 60.0, 45.0]][i][j]
        });
        let sf = vec![1.0, 1.0, 1.0, 1.0];
        let design = Mat::from_fn(4, 2, |i, j| {
            if j == 0 {
                1.0
            } else {
                [0.0, 1.0, 0.0, 1.0][i]
            }
        });
        let dispersions = vec![0.1, 0.1];

        let result = fit_negative_binomial_glm(&counts, &sf, &design, &dispersions);
        assert_eq!(result.coefficients.nrows(), 2);
        assert_eq!(result.coefficients.ncols(), 2);
        // Gene 0: treated (cols 1,3) has higher counts -> positive beta[1]
        assert!(result.coefficients[(0, 1)] > 0.0);
        assert!(result.converged[0]);
        assert!(result.converged[1]);
        // SEs should be positive and finite
        for i in 0..2 {
            for j in 0..2 {
                assert!(result.standard_errors[(i, j)] > 0.0);
                assert!(result.standard_errors[(i, j)].is_finite());
            }
        }
    }

    #[test]
    fn test_irls_known_fold_change() {
        // Gene with exactly 2x counts in treated vs untreated
        // Expected log fold change ~ ln(2) = 0.693
        let counts = Mat::from_fn(1, 4, |_, j| [100.0, 200.0, 100.0, 200.0][j]);
        let sf = vec![1.0, 1.0, 1.0, 1.0];
        let design = Mat::from_fn(4, 2, |i, j| {
            if j == 0 {
                1.0
            } else {
                [0.0, 1.0, 0.0, 1.0][i]
            }
        });
        let dispersions = vec![0.1];
        let result = fit_negative_binomial_glm(&counts, &sf, &design, &dispersions);
        // Beta[1] should be close to ln(2)
        assert_relative_eq!(result.coefficients[(0, 1)], 2.0f64.ln(), epsilon = 0.01);
    }

    #[test]
    fn test_irls_with_size_factors() {
        // Gene with same normalized expression but different size factors
        // Counts = sf * base_expression, so log fold change should be ~0
        let counts = Mat::from_fn(1, 4, |_, j| [50.0, 200.0, 50.0, 200.0][j]);
        let sf = vec![0.5, 2.0, 0.5, 2.0];
        let design = Mat::from_fn(4, 2, |i, j| {
            if j == 0 {
                1.0
            } else {
                [0.0, 1.0, 0.0, 1.0][i]
            }
        });
        let dispersions = vec![0.1];
        let result = fit_negative_binomial_glm(&counts, &sf, &design, &dispersions);
        // After normalization, all samples have ~100 counts, so LFC should be ~0
        assert_relative_eq!(result.coefficients[(0, 1)], 0.0, epsilon = 0.1);
    }

    #[test]
    fn test_irls_deviance_finite() {
        let counts = Mat::from_fn(1, 4, |_, j| [100.0, 100.0, 100.0, 100.0][j]);
        let sf = vec![1.0, 1.0, 1.0, 1.0];
        let design = Mat::from_fn(4, 2, |i, j| {
            if j == 0 {
                1.0
            } else {
                [0.0, 1.0, 0.0, 1.0][i]
            }
        });
        let dispersions = vec![0.1];
        let result = fit_negative_binomial_glm(&counts, &sf, &design, &dispersions);
        assert!(result.deviances[0].is_finite());
        assert!(result.deviances[0] >= 0.0);
    }
}
