// Gene-wise dispersion estimation for DESeq2
//
// Implements the NB log-posterior and its derivatives, plus the Armijo line search
// for fitting per-gene dispersion parameters. Matches DESeq2's C++ code in DESeq2.cpp.

use faer::Mat;
use rayon::prelude::*;
use statrs::function::gamma::{digamma, ln_gamma};

// ---------------------------------------------------------------------------
// Trigamma function (not available in statrs)
// ---------------------------------------------------------------------------

/// Trigamma function: the second derivative of ln(Gamma(x)).
/// Uses recurrence relation for small x and asymptotic expansion for large x.
pub fn trigamma(mut x: f64) -> f64 {
    // For negative integers or zero, trigamma is undefined
    if x <= 0.0 && x == x.floor() {
        return f64::NAN;
    }

    // Reflection formula for x < 0: trigamma(1-x) + pi^2/sin^2(pi*x) = trigamma(x)
    // (We avoid this and just use recurrence to push x positive.)
    let mut result = 0.0;

    // Use recurrence: trigamma(x) = trigamma(x+1) + 1/x^2
    // to push x into region where asymptotic expansion is accurate
    while x < 6.0 {
        result += 1.0 / (x * x);
        x += 1.0;
    }

    // Asymptotic expansion for x >= 6:
    // trigamma(x) ~ 1/x + 1/(2x^2) + 1/(6x^3) - 1/(30x^5) + 1/(42x^7) - 1/(30x^9) + ...
    let x2 = x * x;
    let x_inv = 1.0 / x;
    result += x_inv
        + 1.0 / (2.0 * x2)
        + 1.0 / (6.0 * x2 * x)
        - 1.0 / (30.0 * x2 * x2 * x)
        + 1.0 / (42.0 * x2 * x2 * x2 * x)
        - 1.0 / (30.0 * x2 * x2 * x2 * x2 * x)
        + 5.0 / (66.0 * x2 * x2 * x2 * x2 * x2 * x);

    result
}

// ---------------------------------------------------------------------------
// Small dense matrix helpers (for p x p matrices, typically p <= 4)
// ---------------------------------------------------------------------------

/// Compute determinant of a small dense matrix using cofactor expansion.
fn matrix_determinant(m: &[Vec<f64>]) -> f64 {
    let n = m.len();
    if n == 0 {
        return 1.0;
    }
    if n == 1 {
        return m[0][0];
    }
    if n == 2 {
        return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    }
    // Cofactor expansion along first row
    let mut det = 0.0;
    for j in 0..n {
        let minor: Vec<Vec<f64>> = (1..n)
            .map(|i| {
                (0..n)
                    .filter(|&k| k != j)
                    .map(|k| m[i][k])
                    .collect()
            })
            .collect();
        let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
        det += sign * m[0][j] * matrix_determinant(&minor);
    }
    det
}

/// Compute inverse of a small dense matrix.
fn matrix_inverse(m: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = m.len();
    if n == 2 {
        let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        let inv_det = 1.0 / det;
        return vec![
            vec![m[1][1] * inv_det, -m[0][1] * inv_det],
            vec![-m[1][0] * inv_det, m[0][0] * inv_det],
        ];
    }
    // General case: adjugate / determinant
    let det = matrix_determinant(m);
    let inv_det = 1.0 / det;
    let mut result = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            let minor: Vec<Vec<f64>> = (0..n)
                .filter(|&r| r != j)
                .map(|r| (0..n).filter(|&c| c != i).map(|c| m[r][c]).collect())
                .collect();
            let sign = if (i + j) % 2 == 0 { 1.0 } else { -1.0 };
            result[i][j] = sign * matrix_determinant(&minor) * inv_det;
        }
    }
    result
}

/// Multiply two small dense matrices A (m x k) and B (k x n) -> C (m x n).
fn matrix_multiply(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let m = a.len();
    let k = a[0].len();
    let n = b[0].len();
    let mut c = vec![vec![0.0; n]; m];
    for i in 0..m {
        for j in 0..n {
            let mut sum = 0.0;
            for l in 0..k {
                sum += a[i][l] * b[l][j];
            }
            c[i][j] = sum;
        }
    }
    c
}

/// Compute trace(A * B) = sum_{i,j} A_{ij} * B_{ji} without forming the product.
fn matrix_trace_product(a: &[Vec<f64>], b: &[Vec<f64>]) -> f64 {
    let n = a.len();
    let m = a[0].len();
    let mut tr = 0.0;
    for i in 0..n {
        for j in 0..m {
            tr += a[i][j] * b[j][i];
        }
    }
    tr
}

/// Build X^T * diag(w) * X from a faer Mat and a weight vector.
/// Returns a p x p matrix as Vec<Vec<f64>>.
fn xt_diag_w_x(x: &Mat<f64>, w: &[f64]) -> Vec<Vec<f64>> {
    let n = x.nrows();
    let p = x.ncols();
    let mut result = vec![vec![0.0; p]; p];
    for i in 0..p {
        for j in 0..p {
            let mut sum = 0.0;
            for k in 0..n {
                sum += x[(k, i)] * w[k] * x[(k, j)];
            }
            result[i][j] = sum;
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Log-posterior and derivatives — matching DESeq2.cpp exactly
// ---------------------------------------------------------------------------

/// Negative binomial log-posterior of the dispersion parameter.
///
/// Matches DESeq2.cpp `log_posterior` (lines 31-64), without the weights path.
///
/// # Arguments
/// * `log_alpha` - log of the dispersion parameter
/// * `y` - observed counts for one gene (length m, one per sample)
/// * `mu` - fitted mean values (length m)
/// * `x` - design matrix (m x p)
/// * `prior_mean` - prior mean for log(alpha)
/// * `prior_var` - prior variance for log(alpha)
/// * `use_prior` - whether to include the normal prior term
/// * `use_cr` - whether to include the Cox-Reid adjustment
pub fn log_posterior(
    log_alpha: f64,
    y: &[f64],
    mu: &[f64],
    x: &Mat<f64>,
    prior_mean: f64,
    prior_var: f64,
    use_prior: bool,
    use_cr: bool,
) -> f64 {
    let alpha = log_alpha.exp();
    let alpha_neg1 = 1.0 / alpha;

    // Log-likelihood part
    let ll_part: f64 = y
        .iter()
        .zip(mu.iter())
        .map(|(&yi, &mui)| {
            ln_gamma(yi + alpha_neg1)
                - ln_gamma(alpha_neg1)
                - yi * (mui + alpha_neg1).ln()
                - alpha_neg1 * (1.0 + mui * alpha).ln()
        })
        .sum();

    // Cox-Reid adjustment
    let cr_term = if use_cr {
        let w_diag: Vec<f64> = mu
            .iter()
            .map(|&mui| 1.0 / (1.0 / mui + alpha))
            .collect();
        let b = xt_diag_w_x(x, &w_diag);
        let det_b = matrix_determinant(&b);
        -0.5 * det_b.ln()
    } else {
        0.0
    };

    // Prior part
    let prior_part = if use_prior {
        -0.5 * (log_alpha - prior_mean).powi(2) / prior_var
    } else {
        0.0
    };

    ll_part + prior_part + cr_term
}

/// First derivative of the log-posterior w.r.t. log(alpha).
///
/// Matches DESeq2.cpp `dlog_posterior` (lines 68-107).
pub fn dlog_posterior(
    log_alpha: f64,
    y: &[f64],
    mu: &[f64],
    x: &Mat<f64>,
    prior_mean: f64,
    prior_var: f64,
    use_prior: bool,
    use_cr: bool,
) -> f64 {
    let alpha = log_alpha.exp();
    let alpha_neg1 = 1.0 / alpha;
    let alpha_neg2 = 1.0 / (alpha * alpha);

    // dLL/dalpha
    let ll_part: f64 = alpha_neg2
        * y.iter()
            .zip(mu.iter())
            .map(|(&yi, &mui)| {
                digamma(alpha_neg1) + (1.0 + mui * alpha).ln()
                    - mui * alpha / (1.0 + mui * alpha)
                    - digamma(yi + alpha_neg1)
                    + yi / (mui + alpha_neg1)
            })
            .sum::<f64>();

    // dCR/dalpha
    let cr_term = if use_cr {
        let w_diag: Vec<f64> = mu
            .iter()
            .map(|&mui| 1.0 / (1.0 / mui + alpha))
            .collect();
        let dw_diag: Vec<f64> = mu
            .iter()
            .map(|&mui| -1.0 / (1.0 / mui + alpha).powi(2))
            .collect();
        let b = xt_diag_w_x(x, &w_diag);
        let db = xt_diag_w_x(x, &dw_diag);
        let b_inv = matrix_inverse(&b);
        let det_b = matrix_determinant(&b);
        let ddetb = det_b * matrix_trace_product(&b_inv, &db);
        -0.5 * ddetb / det_b
    } else {
        0.0
    };

    // Prior part (w.r.t. log_alpha)
    let prior_part = if use_prior {
        -1.0 * (log_alpha - prior_mean) / prior_var
    } else {
        0.0
    };

    // Return dlog_post/dalpha * alpha + prior_part
    // because we take derivatives w.r.t. log(alpha)
    (ll_part + cr_term) * alpha + prior_part
}

/// Second derivative of the log-posterior w.r.t. log(alpha).
///
/// Matches DESeq2.cpp `d2log_posterior` (lines 111-158).
pub fn d2log_posterior(
    log_alpha: f64,
    y: &[f64],
    mu: &[f64],
    x: &Mat<f64>,
    prior_mean: f64,
    prior_var: f64,
    use_prior: bool,
    use_cr: bool,
) -> f64 {
    let alpha = log_alpha.exp();
    let alpha_neg1 = 1.0 / alpha;
    let alpha_neg2 = 1.0 / (alpha * alpha);

    // First sum: same as in dlog_posterior (the part inside ll_part before multiplication)
    let sum1: f64 = y
        .iter()
        .zip(mu.iter())
        .map(|(&yi, &mui)| {
            digamma(alpha_neg1) + (1.0 + mui * alpha).ln()
                - mui * alpha / (1.0 + mui * alpha)
                - digamma(yi + alpha_neg1)
                + yi / (mui + alpha_neg1)
        })
        .sum();

    // Second sum: the derivative of the first sum w.r.t alpha
    let sum2: f64 = y
        .iter()
        .zip(mu.iter())
        .map(|(&yi, &mui)| {
            -alpha_neg2 * trigamma(alpha_neg1)
                + mui * mui * alpha / (1.0 + mui * alpha).powi(2)
                + alpha_neg2 * trigamma(yi + alpha_neg1)
                + alpha_neg2 * yi / (mui + alpha_neg1).powi(2)
        })
        .sum();

    let ll_part = -2.0 / (alpha * alpha * alpha) * sum1 + alpha_neg2 * sum2;

    // Cox-Reid second derivative
    let cr_term = if use_cr {
        let w_diag: Vec<f64> = mu
            .iter()
            .map(|&mui| 1.0 / (1.0 / mui + alpha))
            .collect();
        let dw_diag: Vec<f64> = mu
            .iter()
            .map(|&mui| -1.0 / (1.0 / mui + alpha).powi(2))
            .collect();
        let d2w_diag: Vec<f64> = mu
            .iter()
            .map(|&mui| 2.0 / (1.0 / mui + alpha).powi(3))
            .collect();
        let b = xt_diag_w_x(x, &w_diag);
        let b_inv = matrix_inverse(&b);
        let db = xt_diag_w_x(x, &dw_diag);
        let d2b = xt_diag_w_x(x, &d2w_diag);
        let det_b = matrix_determinant(&b);
        let tr_binv_db = matrix_trace_product(&b_inv, &db);
        let ddetb = det_b * tr_binv_db;
        let binv_db = matrix_multiply(&b_inv, &db);
        let d2detb = det_b
            * (tr_binv_db.powi(2) - matrix_trace_product(&binv_db, &binv_db)
                + matrix_trace_product(&b_inv, &d2b));
        0.5 * (ddetb / det_b).powi(2) - 0.5 * d2detb / det_b
    } else {
        0.0
    };

    // Prior part (w.r.t. log_alpha): second derivative of -0.5*(log_alpha-mean)^2/var
    let prior_part = if use_prior {
        -1.0 / prior_var
    } else {
        0.0
    };

    // d2log_post/dlogalpha2 = (d2log_post/dalpha2 * alpha^2 + dlog_post/dlogalpha) + prior
    // where dlog_post/dlogalpha is computed WITHOUT the prior
    let dlog_post_without_prior =
        dlog_posterior(log_alpha, y, mu, x, prior_mean, prior_var, false, use_cr);

    (ll_part + cr_term) * alpha * alpha + dlog_post_without_prior + prior_part
}

// ---------------------------------------------------------------------------
// Dispersion fitting via Armijo line search — matching DESeq2.cpp fitDisp
// ---------------------------------------------------------------------------

/// Result of fitting a single gene's dispersion parameter.
pub struct DispFitResult {
    pub log_alpha: f64,
    pub iterations: usize,
    pub converged: bool,
    pub final_log_posterior: f64,
}

/// Fit the dispersion parameter for a single gene using Armijo line search.
///
/// Matches DESeq2.cpp `fitDisp` (lines 164-277), operating on a single gene.
///
/// # Arguments
/// * `y` - observed counts for one gene
/// * `mu` - fitted means for one gene
/// * `x` - design matrix
/// * `log_alpha_init` - initial log(dispersion)
/// * `prior_mean` - prior mean for log(alpha)
/// * `prior_var` - prior variance for log(alpha)
/// * `use_prior` - whether to include the normal prior
/// * `use_cr` - whether to include Cox-Reid adjustment
pub fn fit_dispersion(
    y: &[f64],
    mu: &[f64],
    x: &Mat<f64>,
    log_alpha_init: f64,
    prior_mean: f64,
    prior_var: f64,
    use_prior: bool,
    use_cr: bool,
) -> DispFitResult {
    let kappa_0: f64 = 1.0;
    let epsilon: f64 = 1e-4;
    let tol: f64 = 1e-6;
    let maxit: usize = 100;
    let min_log_alpha: f64 = -30.0;

    let mut a = log_alpha_init;
    let mut lp = log_posterior(a, y, mu, x, prior_mean, prior_var, use_prior, use_cr);
    let mut dlp = dlog_posterior(a, y, mu, x, prior_mean, prior_var, use_prior, use_cr);
    let mut kappa = kappa_0;
    let mut converged = false;
    let mut iter_accept: usize = 0;
    let mut iterations: usize = 0;

    for _t in 0..maxit {
        iterations += 1;

        let mut a_propose = a + kappa * dlp;

        // Clamp proposals to [-30, 10]
        if a_propose < -30.0 {
            kappa = (-30.0 - a) / dlp;
        }
        if a_propose > 10.0 {
            kappa = (10.0 - a) / dlp;
        }

        // Recompute a_propose with possibly adjusted kappa
        a_propose = a + kappa * dlp;
        let _ = a_propose; // used implicitly via a + kappa*dlp in log_posterior

        let theta_kappa = -1.0
            * log_posterior(
                a + kappa * dlp,
                y,
                mu,
                x,
                prior_mean,
                prior_var,
                use_prior,
                use_cr,
            );
        let theta_hat_kappa = -1.0 * lp - kappa * epsilon * dlp * dlp;

        if theta_kappa <= theta_hat_kappa {
            // Accept step
            iter_accept += 1;
            a += kappa * dlp;
            let lpnew =
                log_posterior(a, y, mu, x, prior_mean, prior_var, use_prior, use_cr);
            let change = lpnew - lp;
            if change < tol {
                lp = lpnew;
                converged = true;
                break;
            }
            if a < min_log_alpha {
                lp = lpnew;
                break;
            }
            lp = lpnew;
            dlp = dlog_posterior(a, y, mu, x, prior_mean, prior_var, use_prior, use_cr);
            // Increase kappa, but cap at kappa_0
            kappa = (kappa * 1.1).min(kappa_0);
            // Every 5 accepts, halve kappa
            if iter_accept % 5 == 0 {
                kappa /= 2.0;
            }
        } else {
            // Reject: halve step size
            kappa /= 2.0;
        }
    }

    DispFitResult {
        log_alpha: a,
        iterations,
        converged,
        final_log_posterior: lp,
    }
}

/// Grid search for dispersion (fallback for non-converged genes).
/// Matches DESeq2.cpp fitDispGrid: evaluate log_posterior on coarse grid,
/// refine around maximum on finer grid.
pub fn fit_dispersion_grid(
    y: &[f64],
    mu: &[f64],
    x: &Mat<f64>,
    n_samples: usize,
    prior_mean: f64,
    prior_var: f64,
    use_prior: bool,
    use_cr: bool,
) -> f64 {
    // Match R's fitDispGridWrapper: 15 points from log(1e-9) to log(max(10, ncol))
    let n_grid = 15;
    let min_la = (1e-8 / 10.0_f64).ln(); // log(1e-9) ≈ -20.7
    let max_la = (10.0f64.max(n_samples as f64)).ln(); // log(10) ≈ 2.3
    let grid: Vec<f64> = (0..n_grid)
        .map(|i| min_la + (max_la - min_la) * i as f64 / (n_grid - 1) as f64)
        .collect();

    // Coarse grid
    let mut best_lp = f64::NEG_INFINITY;
    let mut best_idx = 0;
    for (i, &la) in grid.iter().enumerate() {
        let lp = log_posterior(la, y, mu, x, prior_mean, prior_var, use_prior, use_cr);
        if lp > best_lp {
            best_lp = lp;
            best_idx = i;
        }
    }

    // Fine grid around best
    let delta = if n_grid > 1 { grid[1] - grid[0] } else { 1.0 };
    let center = grid[best_idx];
    let fine_grid: Vec<f64> = (0..n_grid)
        .map(|i| center - delta + (2.0 * delta * i as f64 / (n_grid - 1) as f64))
        .collect();

    let mut best_la = center;
    best_lp = f64::NEG_INFINITY;
    for &la in &fine_grid {
        let lp = log_posterior(la, y, mu, x, prior_mean, prior_var, use_prior, use_cr);
        if lp > best_lp {
            best_lp = lp;
            best_la = la;
        }
    }

    best_la
}

// ---------------------------------------------------------------------------
// Rough dispersion estimate (method of moments)
// ---------------------------------------------------------------------------

/// Compute a rough per-gene dispersion estimate using the method of moments.
///
/// For each gene: `alpha = max(0, (sum_j ((y_j - mu_j)^2 - mu_j) / mu_j^2) / (m - p))`
///
/// # Arguments
/// * `normalized_counts` - genes x samples matrix (counts / size_factors)
/// * `mu` - genes x samples matrix of fitted means
/// * `n_params` - number of parameters in the model (p)
pub fn rough_disp_estimate(
    normalized_counts: &Mat<f64>,
    mu: &Mat<f64>,
    n_params: usize,
) -> Vec<f64> {
    let n_genes = normalized_counts.nrows();
    let n_samples = normalized_counts.ncols();
    let denom = (n_samples - n_params) as f64;

    (0..n_genes)
        .map(|i| {
            let sum: f64 = (0..n_samples)
                .map(|j| {
                    let y = normalized_counts[(i, j)];
                    let m = mu[(i, j)];
                    ((y - m).powi(2) - m) / (m * m)
                })
                .sum();
            (sum / denom).max(0.0)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Linear model mu via QR projection
// ---------------------------------------------------------------------------

/// Compute fitted values mu from a linear model using QR decomposition.
///
/// Computes mu = design * (design' * design)^{-1} * design' * Y' for each gene,
/// which is equivalent to the projection Y * Q * Q' where Q is the Q factor of QR(design).
/// Clamps mu values to be >= 0.5.
///
/// # Arguments
/// * `normalized_counts` - genes x samples matrix
/// * `design` - samples x p design matrix
///
/// # Returns
/// genes x samples matrix of fitted mean values
pub fn linear_model_mu(normalized_counts: &Mat<f64>, design: &Mat<f64>) -> Mat<f64> {
    let n_genes = normalized_counts.nrows();
    let n_samples = normalized_counts.ncols();
    let p = design.ncols();

    // Compute (X'X)^{-1} X' using the normal equations
    // X'X is p x p, X'Y_row is p x 1 for each gene
    // beta = (X'X)^{-1} X' y, then mu = X * beta

    // Build X'X as a dense p x p matrix
    let mut xtx = vec![vec![0.0; p]; p];
    for i in 0..p {
        for j in 0..p {
            let mut sum = 0.0;
            for k in 0..n_samples {
                sum += design[(k, i)] * design[(k, j)];
            }
            xtx[i][j] = sum;
        }
    }
    let xtx_inv = matrix_inverse(&xtx);

    // For each gene: beta = (X'X)^{-1} X' y, then mu = X * beta
    let mut mu = Mat::zeros(n_genes, n_samples);
    for gene in 0..n_genes {
        // X' y
        let mut xty = vec![0.0; p];
        for j in 0..p {
            let mut sum = 0.0;
            for k in 0..n_samples {
                sum += design[(k, j)] * normalized_counts[(gene, k)];
            }
            xty[j] = sum;
        }

        // beta = (X'X)^{-1} * X'y
        let mut beta = vec![0.0; p];
        for i in 0..p {
            let mut sum = 0.0;
            for j in 0..p {
                sum += xtx_inv[i][j] * xty[j];
            }
            beta[i] = sum;
        }

        // mu = X * beta, clamped to >= 0.5
        for k in 0..n_samples {
            let mut val = 0.0;
            for j in 0..p {
                val += design[(k, j)] * beta[j];
            }
            mu[(gene, k)] = val.max(0.5);
        }
    }

    mu
}

// ---------------------------------------------------------------------------
// Orchestrator: estimate_dispersions_gene
// ---------------------------------------------------------------------------

/// Moments-based dispersion estimate (R's momentsDispEstimate).
/// momentsDisp = (rowVar(norm) - mean(1/sf) * rowMean(norm)) / rowMean(norm)^2
pub fn moments_disp_estimate(
    normalized_counts: &Mat<f64>,
    size_factors: &[f64],
) -> Vec<f64> {
    let n_genes = normalized_counts.nrows();
    let n_samples = normalized_counts.ncols();
    let xim: f64 = size_factors.iter().map(|&s| 1.0 / s).sum::<f64>() / n_samples as f64;

    (0..n_genes)
        .map(|i| {
            let mean: f64 = (0..n_samples).map(|j| normalized_counts[(i, j)]).sum::<f64>()
                / n_samples as f64;
            if mean < 1e-8 {
                return 0.0;
            }
            let var: f64 = (0..n_samples)
                .map(|j| (normalized_counts[(i, j)] - mean).powi(2))
                .sum::<f64>()
                / (n_samples - 1) as f64;
            ((var - xim * mean) / (mean * mean)).max(0.0)
        })
        .collect()
}

/// Compute mu from NB GLM fit (unused in current pipeline, retained for future use).
#[allow(dead_code)]
fn compute_mu_from_glm(
    counts: &Mat<f64>,
    size_factors: &[f64],
    design_matrix: &Mat<f64>,
    dispersions: &[f64],
    minmu: f64,
) -> Mat<f64> {
    let n_genes = counts.nrows();
    let n_samples = counts.ncols();
    let n_params = design_matrix.ncols();
    let tol = 1e-6;
    let maxit = 100;
    let large = 30.0;
    let lambda = 1e-6;

    // Initialize betas via normal equations on log(normalized + 0.1)
    // Build (X'X)^{-1} X'
    let p = n_params;
    let mut xtx = vec![vec![0.0; p]; p];
    for i in 0..p {
        for j in 0..p {
            for k in 0..n_samples {
                xtx[i][j] += design_matrix[(k, i)] * design_matrix[(k, j)];
            }
        }
    }
    let xtx_inv = matrix_inverse(&xtx);

    let mut mu = Mat::zeros(n_genes, n_samples);

    for gene in 0..n_genes {
        let alpha = dispersions[gene];
        if alpha.is_nan() || alpha <= 0.0 {
            for j in 0..n_samples {
                mu[(gene, j)] = (counts[(gene, j)] / size_factors[j]).max(minmu) * size_factors[j];
            }
            continue;
        }

        // Init beta from log(normalized + 0.1)
        let mut xty = vec![0.0; p];
        for k in 0..p {
            for j in 0..n_samples {
                xty[k] += design_matrix[(j, k)] * (counts[(gene, j)] / size_factors[j] + 0.1).ln();
            }
        }
        let mut beta = vec![0.0; p];
        for i in 0..p {
            for j in 0..p {
                beta[i] += xtx_inv[i][j] * xty[j];
            }
        }

        // Compute initial mu
        let mut mu_g: Vec<f64> = (0..n_samples)
            .map(|j| {
                let eta: f64 = (0..p).map(|k| design_matrix[(j, k)] * beta[k]).sum();
                (size_factors[j] * eta.exp()).max(minmu)
            })
            .collect();

        // IRLS
        let mut dev_old = 0.0;
        for t in 0..maxit {
            let w: Vec<f64> = mu_g.iter().map(|&m| m / (1.0 + alpha * m)).collect();
            let z: Vec<f64> = (0..n_samples)
                .map(|j| (mu_g[j] / size_factors[j]).ln() + (counts[(gene, j)] - mu_g[j]) / mu_g[j])
                .collect();

            // Solve (X'WX + L) beta = X'Wz
            let mut xtwx = vec![vec![0.0; p]; p];
            let mut xtwz = vec![0.0; p];
            for a in 0..p {
                for b in 0..p {
                    for j in 0..n_samples {
                        xtwx[a][b] += design_matrix[(j, a)] * w[j] * design_matrix[(j, b)];
                    }
                }
                xtwx[a][a] += lambda;
                for j in 0..n_samples {
                    xtwz[a] += design_matrix[(j, a)] * w[j] * z[j];
                }
            }

            let xtwx_inv = matrix_inverse(&xtwx);
            let mut new_beta = vec![0.0; p];
            for i in 0..p {
                for j in 0..p {
                    new_beta[i] += xtwx_inv[i][j] * xtwz[j];
                }
            }

            if new_beta.iter().any(|&b| b.abs() > large) { break; }
            beta = new_beta;

            mu_g = (0..n_samples)
                .map(|j| {
                    let eta: f64 = (0..p).map(|k| design_matrix[(j, k)] * beta[k]).sum();
                    (size_factors[j] * eta.exp()).max(minmu)
                })
                .collect();

            let dev: f64 = -2.0 * (0..n_samples)
                .map(|j| {
                    let size = 1.0 / alpha;
                    let y = counts[(gene, j)];
                    ln_gamma(y + size) - ln_gamma(size) - ln_gamma(y + 1.0)
                        + y * (mu_g[j] / (mu_g[j] + size)).ln()
                        + size * (size / (mu_g[j] + size)).ln()
                })
                .sum::<f64>();

            let conv = (dev - dev_old).abs() / (dev.abs() + 0.1);
            if conv.is_nan() { break; }
            if t > 0 && conv < tol { break; }
            dev_old = dev;
        }

        for j in 0..n_samples {
            mu[(gene, j)] = mu_g[j];
        }
    }

    mu
}

/// Estimate gene-wise dispersions via maximum likelihood with Armijo line search.
///
/// Matches R's estimateDispersionsGeneEst:
/// 1. Filter all-zero genes (they get NaN dispersion)
/// 2. Normalize counts, compute mu via linear model, scale to raw
/// 3. Combine rough + moments dispersion estimates as starting values
/// 4. Per-gene line search optimization with convergence check (parallel via rayon)
pub fn estimate_dispersions_gene(
    counts: &Mat<f64>,
    size_factors: &[f64],
    design_matrix: &Mat<f64>,
) -> (Vec<f64>, Mat<f64>) {
    let n_genes = counts.nrows();
    let n_samples = counts.ncols();

    // Identify all-zero genes
    let all_zero: Vec<bool> = (0..n_genes)
        .map(|i| (0..n_samples).all(|j| counts[(i, j)] == 0.0))
        .collect();

    // Step 1: Normalize counts
    let normalized = Mat::from_fn(n_genes, n_samples, |i, j| {
        counts[(i, j)] / size_factors[j]
    });

    // Step 2: Compute mu on normalized scale, then scale to raw counts
    let mu_normalized = linear_model_mu(&normalized, design_matrix);
    let mu = Mat::from_fn(n_genes, n_samples, |i, j| {
        (mu_normalized[(i, j)] * size_factors[j]).max(0.5)
    });

    // Step 3: Starting dispersion = min(rough, moments), clamped
    let n_params = design_matrix.ncols();
    let rough_disps = rough_disp_estimate(&normalized, &mu_normalized, n_params);
    let moments_disps = moments_disp_estimate(&normalized, size_factors);

    let min_disp = 1e-8;
    let max_disp = 10.0f64.max(n_samples as f64);

    // Step 4: Per-gene optimization (parallel), with convergence check
    let dispersions: Vec<f64> = (0..n_genes)
        .into_par_iter()
        .map(|gene| {
            if all_zero[gene] {
                return f64::NAN;
            }

            let y: Vec<f64> = (0..n_samples).map(|j| counts[(gene, j)]).collect();
            let mu_row: Vec<f64> = (0..n_samples).map(|j| mu[(gene, j)]).collect();

            let alpha_init = rough_disps[gene]
                .min(moments_disps[gene])
                .clamp(min_disp, max_disp);
            let log_alpha_init = alpha_init.ln();

            // Compute initial log-posterior for convergence check
            let initial_lp = log_posterior(
                log_alpha_init, &y, &mu_row, design_matrix, log_alpha_init, 1.0, false, true,
            );

            let result = fit_dispersion(
                &y,
                &mu_row,
                design_matrix,
                log_alpha_init,
                log_alpha_init,
                1.0,
                false,
                true,
            );

            // R's convergence check: revert if log-posterior didn't increase
            if result.final_log_posterior < initial_lp + initial_lp.abs() * 1e-6 {
                return alpha_init;
            }

            result.log_alpha.exp().clamp(min_disp, max_disp)
        })
        .collect();

    (dispersions, mu)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn make_test_data() -> (Vec<f64>, Vec<f64>, Mat<f64>) {
        let y = vec![10.0, 20.0, 15.0, 25.0];
        let mu = vec![12.0, 18.0, 14.0, 22.0];
        let x = Mat::from_fn(4, 2, |i, j| {
            if j == 0 {
                1.0
            } else {
                [0.0, 1.0, 0.0, 1.0][i]
            }
        });
        (y, mu, x)
    }

    #[test]
    fn test_log_posterior_finite() {
        let (y, mu, x) = make_test_data();
        let lp = log_posterior(0.0, &y, &mu, &x, 0.0, 1.0, false, true);
        assert!(lp.is_finite());
    }

    #[test]
    fn test_derivative_numerical_check() {
        let (y, mu, x) = make_test_data();
        let log_alpha = 0.0;
        let eps = 1e-5;
        let lp_plus = log_posterior(log_alpha + eps, &y, &mu, &x, 0.0, 1.0, false, true);
        let lp_minus = log_posterior(log_alpha - eps, &y, &mu, &x, 0.0, 1.0, false, true);
        let numerical = (lp_plus - lp_minus) / (2.0 * eps);
        let analytical = dlog_posterior(log_alpha, &y, &mu, &x, 0.0, 1.0, false, true);
        assert_relative_eq!(analytical, numerical, epsilon = 1e-4);
    }

    #[test]
    fn test_second_derivative_numerical_check() {
        let (y, mu, x) = make_test_data();
        let log_alpha = 0.0;
        let eps = 1e-5;
        let d1_plus = dlog_posterior(log_alpha + eps, &y, &mu, &x, 0.0, 1.0, false, true);
        let d1_minus = dlog_posterior(log_alpha - eps, &y, &mu, &x, 0.0, 1.0, false, true);
        let numerical = (d1_plus - d1_minus) / (2.0 * eps);
        let analytical = d2log_posterior(log_alpha, &y, &mu, &x, 0.0, 1.0, false, true);
        assert_relative_eq!(analytical, numerical, epsilon = 1e-3);
    }

    #[test]
    fn test_fit_dispersion_converges() {
        // Use data with clear overdispersion so the MLE is finite
        let y = vec![5.0, 50.0, 3.0, 60.0];
        let mu = vec![15.0, 15.0, 15.0, 15.0];
        let x = Mat::from_fn(4, 2, |i, j| {
            if j == 0 {
                1.0
            } else {
                [0.0, 1.0, 0.0, 1.0][i]
            }
        });
        let result = fit_dispersion(&y, &mu, &x, 0.0, 0.0, 1.0, false, true);
        assert!(result.log_alpha.is_finite());
        assert!(result.log_alpha > -30.0 && result.log_alpha < 10.0);
        assert!(result.converged);
    }

    #[test]
    fn test_trigamma_values() {
        // trigamma(1) = pi^2/6
        assert_relative_eq!(
            trigamma(1.0),
            std::f64::consts::PI.powi(2) / 6.0,
            epsilon = 1e-10
        );
        // trigamma(2) = pi^2/6 - 1
        assert_relative_eq!(
            trigamma(2.0),
            std::f64::consts::PI.powi(2) / 6.0 - 1.0,
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_log_posterior_no_cr() {
        let (y, mu, x) = make_test_data();
        let lp = log_posterior(0.0, &y, &mu, &x, 0.0, 1.0, false, false);
        assert!(lp.is_finite());
    }

    #[test]
    fn test_log_posterior_with_prior() {
        let (y, mu, x) = make_test_data();
        let lp_no_prior = log_posterior(0.0, &y, &mu, &x, 0.0, 1.0, false, true);
        let lp_with_prior = log_posterior(0.0, &y, &mu, &x, 0.0, 1.0, true, true);
        // With prior at mean=0, log_alpha=0, the prior term is 0
        assert_relative_eq!(lp_no_prior, lp_with_prior, epsilon = 1e-10);

        // With a different prior mean, the posterior should differ
        let lp_shifted = log_posterior(0.0, &y, &mu, &x, 2.0, 1.0, true, true);
        assert!(lp_shifted < lp_no_prior);
    }

    #[test]
    fn test_rough_disp_estimate() {
        // Simple test: if y == mu, dispersion should be near 0
        let n = Mat::from_fn(2, 4, |i, j| {
            [[10.0, 20.0, 15.0, 25.0], [30.0, 60.0, 45.0, 75.0]][i][j]
        });
        let mu = n.clone();
        let disps = rough_disp_estimate(&n, &mu, 2);
        // When y == mu, (y-mu)^2 - mu = -mu, which is negative, so clamped to 0
        // Actually: (y-mu)^2 = 0, so sum = sum(-mu/mu^2) = sum(-1/mu) which is negative
        for d in &disps {
            assert!(*d >= 0.0);
        }
    }

    #[test]
    fn test_linear_model_mu() {
        let counts = Mat::from_fn(2, 4, |i, j| {
            [[10.0, 20.0, 15.0, 25.0], [5.0, 10.0, 7.0, 12.0]][i][j]
        });
        let design = Mat::from_fn(4, 2, |i, j| {
            if j == 0 {
                1.0
            } else {
                [0.0, 1.0, 0.0, 1.0][i]
            }
        });
        let mu = linear_model_mu(&counts, &design);
        assert_eq!(mu.nrows(), 2);
        assert_eq!(mu.ncols(), 4);
        // All fitted values should be positive
        for i in 0..2 {
            for j in 0..4 {
                assert!(mu[(i, j)] >= 0.5);
            }
        }
    }

    #[test]
    fn test_matrix_determinant() {
        let m = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        assert_relative_eq!(matrix_determinant(&m), -2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_matrix_inverse() {
        let m = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let inv = matrix_inverse(&m);
        assert_relative_eq!(inv[0][0], -2.0, epsilon = 1e-10);
        assert_relative_eq!(inv[0][1], 1.0, epsilon = 1e-10);
        assert_relative_eq!(inv[1][0], 1.5, epsilon = 1e-10);
        assert_relative_eq!(inv[1][1], -0.5, epsilon = 1e-10);
    }

    #[test]
    fn test_derivative_with_prior() {
        let (y, mu, x) = make_test_data();
        let log_alpha = 1.0;
        let eps = 1e-5;
        let lp_plus =
            log_posterior(log_alpha + eps, &y, &mu, &x, 0.0, 1.0, true, true);
        let lp_minus =
            log_posterior(log_alpha - eps, &y, &mu, &x, 0.0, 1.0, true, true);
        let numerical = (lp_plus - lp_minus) / (2.0 * eps);
        let analytical = dlog_posterior(log_alpha, &y, &mu, &x, 0.0, 1.0, true, true);
        assert_relative_eq!(analytical, numerical, epsilon = 1e-4);
    }

    #[test]
    fn test_estimate_dispersions_gene() {
        // Small synthetic dataset
        let counts = Mat::from_fn(3, 4, |i, j| {
            [
                [100.0, 120.0, 80.0, 90.0],
                [200.0, 250.0, 180.0, 220.0],
                [50.0, 60.0, 40.0, 55.0],
            ][i][j]
        });
        let size_factors = vec![1.0, 1.1, 0.9, 1.0];
        let design = Mat::from_fn(4, 2, |i, j| {
            if j == 0 {
                1.0
            } else {
                [0.0, 0.0, 1.0, 1.0][i]
            }
        });
        let (disps, mu) = estimate_dispersions_gene(&counts, &size_factors, &design);
        assert_eq!(disps.len(), 3);
        assert_eq!(mu.nrows(), 3);
        assert_eq!(mu.ncols(), 4);
        for d in &disps {
            assert!(d.is_finite());
            assert!(*d > 0.0);
        }
    }

    #[test]
    fn test_trigamma_large_x() {
        // For large x, trigamma(x) ~ 1/x
        let x = 100.0;
        let tg = trigamma(x);
        assert_relative_eq!(tg, 1.0 / x, epsilon = 1e-3);
    }

    #[test]
    fn test_trigamma_recurrence() {
        // trigamma(x) = trigamma(x+1) + 1/x^2
        let x = 3.5;
        assert_relative_eq!(
            trigamma(x),
            trigamma(x + 1.0) + 1.0 / (x * x),
            epsilon = 1e-12
        );
    }
}
