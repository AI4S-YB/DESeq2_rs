use crate::data::DataFrame;
use faer::Mat;

/// Build a design matrix from a DataFrame column.
/// Creates intercept + indicator columns for non-reference levels (sorted alphabetically).
pub fn build_design_matrix(
    col_data: &DataFrame,
    column: &str,
    reference: &str,
) -> Result<Mat<f64>, String> {
    let values = col_data
        .get_column(column)
        .ok_or_else(|| format!("column '{}' not found in coldata", column))?;
    let levels = col_data
        .levels(column)
        .ok_or_else(|| format!("column '{}' not found", column))?;

    if !levels.contains(&reference.to_string()) {
        return Err(format!("reference level '{}' not found in column '{}'", reference, column));
    }

    let non_ref: Vec<&String> = levels.iter().filter(|l| l.as_str() != reference).collect();
    let n_samples = values.len();
    let n_cols = 1 + non_ref.len();

    let mat = Mat::from_fn(n_samples, n_cols, |i, j| {
        if j == 0 {
            1.0
        } else {
            if values[i] == *non_ref[j - 1] { 1.0 } else { 0.0 }
        }
    });

    Ok(mat)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::DataFrame;

    #[test]
    fn test_design_matrix_two_levels() {
        let mut df = DataFrame::new();
        df.add_column("dex", vec![
            "untrt".into(), "trt".into(), "untrt".into(), "trt".into()
        ]);
        let mat = build_design_matrix(&df, "dex", "untrt").unwrap();
        assert_eq!(mat.nrows(), 4);
        assert_eq!(mat.ncols(), 2);
        assert_eq!(mat[(0, 0)], 1.0); assert_eq!(mat[(0, 1)], 0.0);
        assert_eq!(mat[(1, 0)], 1.0); assert_eq!(mat[(1, 1)], 1.0);
        assert_eq!(mat[(2, 0)], 1.0); assert_eq!(mat[(2, 1)], 0.0);
        assert_eq!(mat[(3, 0)], 1.0); assert_eq!(mat[(3, 1)], 1.0);
    }

    #[test]
    fn test_design_matrix_bad_reference() {
        let mut df = DataFrame::new();
        df.add_column("dex", vec!["untrt".into(), "trt".into()]);
        assert!(build_design_matrix(&df, "dex", "missing").is_err());
    }
}
