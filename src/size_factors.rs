use faer::Mat;

/// Estimate size factors using the median-of-ratios method.
pub fn estimate_size_factors(counts: &Mat<f64>) -> Vec<f64> {
    let n_genes = counts.nrows();
    let n_samples = counts.ncols();

    // Step 1: log geometric mean per gene = rowMeans(log(counts))
    let log_geo_means: Vec<f64> = (0..n_genes)
        .map(|i| {
            let sum: f64 = (0..n_samples).map(|j| counts[(i, j)].ln()).sum();
            sum / n_samples as f64
        })
        .collect();

    // Step 2: for each sample, collect log-ratios and take median
    (0..n_samples)
        .map(|j| {
            let mut ratios: Vec<f64> = (0..n_genes)
                .filter(|&i| log_geo_means[i].is_finite() && counts[(i, j)] > 0.0)
                .map(|i| counts[(i, j)].ln() - log_geo_means[i])
                .collect();
            ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let median = if ratios.is_empty() {
                1.0
            } else {
                let mid = ratios.len() / 2;
                if ratios.len() % 2 == 0 {
                    (ratios[mid - 1] + ratios[mid]) / 2.0
                } else {
                    ratios[mid]
                }
            };
            median.exp()
        })
        .collect()
}

/// Compute base mean per gene after size factor normalization.
pub fn base_means(counts: &Mat<f64>, size_factors: &[f64]) -> Vec<f64> {
    let n_genes = counts.nrows();
    let n_samples = counts.ncols();
    (0..n_genes)
        .map(|i| {
            let sum: f64 = (0..n_samples)
                .map(|j| counts[(i, j)] / size_factors[j])
                .sum();
            sum / n_samples as f64
        })
        .collect()
}

/// Compute size-factor normalized counts (`counts(dds, normalized=TRUE)` in DESeq2).
pub fn normalized_counts(counts: &Mat<f64>, size_factors: &[f64]) -> Mat<f64> {
    let n_genes = counts.nrows();
    let n_samples = counts.ncols();
    assert_eq!(
        size_factors.len(),
        n_samples,
        "size_factors length must match count matrix columns"
    );

    Mat::from_fn(n_genes, n_samples, |i, j| counts[(i, j)] / size_factors[j])
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use faer::Mat;

    #[test]
    fn test_size_factors_simple() {
        // 3 genes, 2 samples — all ratios are 0.5 for sample 0, 2.0 for sample 1
        let counts = Mat::from_fn(3, 2, |i, j| {
            [[10.0, 20.0], [30.0, 60.0], [50.0, 100.0]][i][j]
        });
        let sf = estimate_size_factors(&counts);
        assert_eq!(sf.len(), 2);
        assert_relative_eq!(sf[0] / sf[1], 0.5, epsilon = 1e-10);
    }

    #[test]
    fn test_size_factors_with_zeros() {
        // Gene with zero excluded from geo mean
        let counts = Mat::from_fn(3, 2, |i, j| {
            [[0.0, 20.0], [30.0, 60.0], [50.0, 100.0]][i][j]
        });
        let sf = estimate_size_factors(&counts);
        assert_eq!(sf.len(), 2);
        assert_relative_eq!(sf[0] / sf[1], 0.5, epsilon = 1e-10);
    }

    #[test]
    fn test_base_means() {
        let counts = Mat::from_fn(2, 2, |i, j| [[10.0, 20.0], [30.0, 60.0]][i][j]);
        let sf = vec![1.0, 2.0];
        let bm = base_means(&counts, &sf);
        assert_relative_eq!(bm[0], 10.0, epsilon = 1e-10); // (10/1 + 20/2) / 2 = 10
        assert_relative_eq!(bm[1], 30.0, epsilon = 1e-10); // (30/1 + 60/2) / 2 = 30
    }

    #[test]
    fn test_normalized_counts() {
        let counts = Mat::from_fn(2, 2, |i, j| [[10.0, 20.0], [30.0, 60.0]][i][j]);
        let sf = vec![1.0, 2.0];
        let normalized = normalized_counts(&counts, &sf);
        assert_relative_eq!(normalized[(0, 0)], 10.0, epsilon = 1e-10);
        assert_relative_eq!(normalized[(0, 1)], 10.0, epsilon = 1e-10);
        assert_relative_eq!(normalized[(1, 0)], 30.0, epsilon = 1e-10);
        assert_relative_eq!(normalized[(1, 1)], 30.0, epsilon = 1e-10);
    }
}
