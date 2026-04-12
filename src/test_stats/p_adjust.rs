/// Benjamini-Hochberg p-value adjustment.
/// NaN p-values remain NaN and are excluded from the count.
pub fn p_adjust_bh(p_values: &[f64]) -> Vec<f64> {
    let n = p_values.len();
    let mut result = vec![f64::NAN; n];

    // Collect non-NaN (index, p-value) pairs
    let mut indexed: Vec<(usize, f64)> = p_values
        .iter()
        .enumerate()
        .filter(|(_, p)| !p.is_nan())
        .map(|(i, &p)| (i, p))
        .collect();

    if indexed.is_empty() {
        return result;
    }
    let m = indexed.len();

    // Sort by p-value descending
    indexed.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    // Adjust with cumulative minimum (ensures monotonicity)
    let mut cummin = f64::INFINITY;
    for (rank_from_top, &(orig_idx, p)) in indexed.iter().enumerate() {
        let rank_from_bottom = m - rank_from_top;
        let adjusted = (p * m as f64 / rank_from_bottom as f64).min(1.0);
        cummin = cummin.min(adjusted);
        result[orig_idx] = cummin;
    }
    result
}

/// Independent filtering: find optimal baseMean threshold, then apply BH.
/// Returns adjusted p-values (NaN for filtered-out genes).
pub fn independent_filtering(base_means: &[f64], p_values: &[f64], alpha: f64) -> Vec<f64> {
    let n = base_means.len();

    // Quantile grid: 50 points from fraction_zero to 0.95
    let fraction_zero =
        base_means.iter().filter(|&&m| m == 0.0).count() as f64 / n as f64;
    let upper = if fraction_zero < 0.95 { 0.95 } else { 1.0 };
    let n_theta = 50;
    let thetas: Vec<f64> = (0..n_theta)
        .map(|i| {
            fraction_zero
                + (upper - fraction_zero) * i as f64 / (n_theta - 1).max(1) as f64
        })
        .collect();

    // Sorted means for quantile computation
    let mut sorted_means: Vec<f64> = base_means.to_vec();
    sorted_means.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let quantile = |q: f64| -> f64 {
        if sorted_means.is_empty() {
            return 0.0;
        }
        if q <= 0.0 {
            return sorted_means[0];
        }
        if q >= 1.0 {
            return sorted_means[n - 1];
        }
        let idx = q * (n - 1) as f64;
        let lo = idx.floor() as usize;
        let hi = (lo + 1).min(n - 1);
        let frac = idx - lo as f64;
        sorted_means[lo] * (1.0 - frac) + sorted_means[hi] * frac
    };

    // For each theta, count rejections after filtering
    let rejections: Vec<usize> = thetas
        .iter()
        .map(|&theta| {
            let threshold = quantile(theta);
            let filtered_p: Vec<f64> = (0..n)
                .map(|i| {
                    if base_means[i] >= threshold {
                        p_values[i]
                    } else {
                        f64::NAN
                    }
                })
                .collect();
            let padj = p_adjust_bh(&filtered_p);
            padj.iter().filter(|&&p| !p.is_nan() && p < alpha).count()
        })
        .collect();

    // Find best threshold (maximizes rejections)
    let mut best_idx = 0;
    if let Some(&max_rej) = rejections.iter().max() {
        if max_rej > 10 {
            best_idx = rejections.iter().position(|&r| r == max_rej).unwrap_or(0);
        }
    }

    let threshold = quantile(thetas[best_idx]);
    let filtered_p: Vec<f64> = (0..n)
        .map(|i| {
            if base_means[i] >= threshold {
                p_values[i]
            } else {
                f64::NAN
            }
        })
        .collect();
    p_adjust_bh(&filtered_p)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_bh_adjustment() {
        let pvals = vec![0.01, 0.04, 0.03, 0.005];
        let padj = p_adjust_bh(&pvals);
        // Sorted: 0.005(3), 0.01(0), 0.03(2), 0.04(1)
        // Adjusted: 0.02, 0.02, 0.04, 0.04
        assert_relative_eq!(padj[3], 0.02, epsilon = 1e-10);
        assert_relative_eq!(padj[0], 0.02, epsilon = 1e-10);
        assert_relative_eq!(padj[2], 0.04, epsilon = 1e-10);
        assert_relative_eq!(padj[1], 0.04, epsilon = 1e-10);
    }

    #[test]
    fn test_bh_with_nan() {
        let pvals = vec![0.01, f64::NAN, 0.03];
        let padj = p_adjust_bh(&pvals);
        assert!(!padj[0].is_nan());
        assert!(padj[1].is_nan());
        assert!(!padj[2].is_nan());
        assert_relative_eq!(padj[0], 0.02, epsilon = 1e-10);
        assert_relative_eq!(padj[2], 0.03, epsilon = 1e-10);
    }

    #[test]
    fn test_independent_filtering_no_filter() {
        // All genes have same baseMean — filtering should not help
        let means = vec![100.0; 10];
        let pvals: Vec<f64> = (0..10).map(|i| (i as f64 + 1.0) / 100.0).collect();
        let padj = independent_filtering(&means, &pvals, 0.1);
        // Should be same as plain BH
        let bh = p_adjust_bh(&pvals);
        for i in 0..10 {
            assert_relative_eq!(padj[i], bh[i], epsilon = 1e-10);
        }
    }
}
