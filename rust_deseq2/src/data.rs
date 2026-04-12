use faer::Mat;
use std::collections::HashMap;

/// Lightweight sample metadata container.
#[derive(Debug, Clone)]
pub struct DataFrame {
    columns: HashMap<String, Vec<String>>,
    pub row_names: Vec<String>,
}

impl DataFrame {
    pub fn new() -> Self {
        Self {
            columns: HashMap::new(),
            row_names: Vec::new(),
        }
    }

    pub fn add_column(&mut self, name: &str, values: Vec<String>) {
        self.columns.insert(name.to_string(), values);
    }

    pub fn get_column(&self, name: &str) -> Option<&[String]> {
        self.columns.get(name).map(|v| v.as_slice())
    }

    /// Return sorted unique values for a column.
    pub fn levels(&self, name: &str) -> Option<Vec<String>> {
        self.columns.get(name).map(|v| {
            let mut lvls: Vec<String> = v.iter().cloned().collect::<std::collections::BTreeSet<_>>().into_iter().collect();
            lvls.sort();
            lvls
        })
    }

    pub fn nrow(&self) -> usize {
        self.columns.values().next().map_or(0, |v| v.len())
    }
}

/// Dispersion estimation results.
#[derive(Debug, Clone)]
pub struct DispersionResult {
    pub gene_estimates: Vec<f64>,
    pub trend_values: Vec<f64>,
    pub map_estimates: Vec<f64>,
    pub prior_var: f64,
}

/// Parametric trend coefficients: dispersion = asympt_disp + extra_pois / mean
#[derive(Debug, Clone, Copy)]
pub struct DispTrendCoeffs {
    pub asympt_disp: f64,
    pub extra_pois: f64,
}

impl DispTrendCoeffs {
    pub fn eval(&self, mean: f64) -> f64 {
        self.asympt_disp + self.extra_pois / mean
    }
}

/// GLM fitting result for all genes.
#[derive(Debug, Clone)]
pub struct BetaFitResult {
    pub coefficients: Mat<f64>,    // genes x p, natural log scale
    pub standard_errors: Mat<f64>, // genes x p
    pub deviances: Vec<f64>,
    pub converged: Vec<bool>,
}

/// Single-gene Wald test output.
#[derive(Debug, Clone)]
pub struct WaldTestResult {
    pub log2_fold_change: f64,
    pub lfc_se: f64,
    pub stat: f64,
    pub p_value: f64,
}

/// Final per-gene result row.
#[derive(Debug, Clone)]
pub struct DESeqResult {
    pub gene: String,
    pub base_mean: f64,
    pub log2_fold_change: f64,
    pub lfc_se: f64,
    pub stat: f64,
    pub p_value: f64,
    pub p_adjusted: f64,
}

/// Contrast specification for hypothesis testing.
#[derive(Debug, Clone)]
pub enum Contrast {
    LastCoefficient,
    ColumnLevels(String, String, String),
    Vector(Vec<f64>),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dataframe_get_column() {
        let mut df = DataFrame::new();
        df.add_column("condition", vec!["untrt".into(), "trt".into(), "untrt".into()]);
        assert_eq!(df.get_column("condition").unwrap(), &["untrt", "trt", "untrt"]);
        assert!(df.get_column("missing").is_none());
    }

    #[test]
    fn test_dataframe_levels() {
        let mut df = DataFrame::new();
        df.add_column("dex", vec!["trt".into(), "untrt".into(), "trt".into(), "untrt".into()]);
        let levels = df.levels("dex").unwrap();
        assert!(levels.contains(&"trt".to_string()));
        assert!(levels.contains(&"untrt".to_string()));
        assert_eq!(levels.len(), 2);
    }
}
