use crate::data::{DataFrame, DESeqResult};
use faer::Mat;
use std::error::Error;
use std::io::BufRead;
use std::path::Path;

/// Read a tab-separated count matrix. First column is gene names, first row is sample names.
pub fn read_count_matrix(path: &Path) -> Result<(Mat<f64>, Vec<String>, Vec<String>), Box<dyn Error>> {
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    let mut lines = reader.lines();

    let header = lines.next().ok_or("empty file")??;
    let sample_names: Vec<String> = header.split('\t').skip(1).map(|s| s.trim().to_string()).collect();
    let n_samples = sample_names.len();

    let mut gene_names = Vec::new();
    let mut data: Vec<Vec<f64>> = Vec::new();

    for line in lines {
        let line = line?;
        if line.is_empty() { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        gene_names.push(fields[0].trim().to_string());
        let row: Vec<f64> = fields[1..]
            .iter()
            .map(|s| s.trim().parse::<f64>())
            .collect::<Result<_, _>>()?;
        if row.len() != n_samples {
            return Err(format!("row {} has {} values, expected {}", gene_names.last().unwrap(), row.len(), n_samples).into());
        }
        data.push(row);
    }

    let n_genes = gene_names.len();
    let counts = Mat::from_fn(n_genes, n_samples, |i, j| data[i][j]);
    Ok((counts, gene_names, sample_names))
}

/// Read tab-separated coldata. First column is sample names (row names).
pub fn read_coldata(path: &Path) -> Result<(DataFrame, Vec<String>), Box<dyn Error>> {
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    let mut lines = reader.lines();

    let header = lines.next().ok_or("empty file")??;
    let col_names: Vec<String> = header.split('\t').skip(1).map(|s| s.trim().to_string()).collect();

    let mut sample_names = Vec::new();
    let mut columns: Vec<Vec<String>> = vec![Vec::new(); col_names.len()];

    for line in lines {
        let line = line?;
        if line.is_empty() { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        sample_names.push(fields[0].trim().to_string());
        for (i, val) in fields[1..].iter().enumerate() {
            if i < columns.len() {
                columns[i].push(val.trim().to_string());
            }
        }
    }

    let mut df = DataFrame::new();
    df.row_names = sample_names.clone();
    for (name, vals) in col_names.iter().zip(columns) {
        df.add_column(name, vals);
    }
    Ok((df, sample_names))
}

/// Write results to a TSV file.
pub fn write_results(path: &Path, results: &[DESeqResult]) -> Result<(), Box<dyn Error>> {
    let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    wtr.write_record(&["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"])?;
    for r in results {
        wtr.write_record(&[
            &r.gene,
            &format_float(r.base_mean),
            &format_float(r.log2_fold_change),
            &format_float(r.lfc_se),
            &format_float(r.stat),
            &format_float(r.p_value),
            &format_float(r.p_adjusted),
        ])?;
    }
    Ok(())
}

fn format_float(v: f64) -> String {
    if v.is_nan() { "NA".to_string() } else { format!("{:.10e}", v) }
}
