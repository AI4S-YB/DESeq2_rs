use deseq2_rs::io::read_count_matrix;
use faer::Mat;
use std::path::Path;
use std::process::Command;

fn read_r_size_factors(path: &Path) -> (Vec<String>, Vec<f64>) {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .unwrap();
    let mut sample_names = Vec::new();
    let mut size_factors = Vec::new();

    for record in rdr.records() {
        let record = record.unwrap();
        sample_names.push(record[0].to_string());
        size_factors.push(record[1].parse::<f64>().unwrap());
    }

    (sample_names, size_factors)
}

fn assert_matrix_close(
    actual: &Mat<f64>,
    expected: &Mat<f64>,
    gene_names: &[String],
    sample_names: &[String],
) {
    assert_eq!(actual.nrows(), expected.nrows());
    assert_eq!(actual.ncols(), expected.ncols());

    for i in 0..expected.nrows() {
        for j in 0..expected.ncols() {
            let expected_value = expected[(i, j)];
            let actual_value = actual[(i, j)];
            let abs_diff = (actual_value - expected_value).abs();
            let rel_diff = if expected_value.abs() < 1.0 {
                abs_diff
            } else {
                abs_diff / expected_value.abs()
            };
            assert!(
                abs_diff <= 1e-8 || rel_diff <= 1e-8,
                "normalized count mismatch for gene {} sample {}: expected {}, got {}",
                gene_names[i],
                sample_names[j],
                expected_value,
                actual_value
            );
        }
    }
}

#[test]
fn test_cli_normalized_counts_match_deseq2_airway_reference() {
    let tmp = tempfile::tempdir().unwrap();
    let results_path = tmp.path().join("results.tsv");
    let normalized_path = tmp.path().join("normalized_counts.tsv");

    let output = Command::new(env!("CARGO_BIN_EXE_deseq2-rs"))
        .args([
            "--counts",
            "tests/reference_data/airway_counts.tsv",
            "--coldata",
            "tests/reference_data/airway_coldata.tsv",
            "--design",
            "dex",
            "--reference",
            "untrt",
            "--output",
            results_path.to_str().unwrap(),
            "--normalized-counts-output",
            normalized_path.to_str().unwrap(),
            "--threads",
            "1",
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "deseq2-rs failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let (raw_counts, raw_gene_names, raw_sample_names) =
        read_count_matrix(Path::new("tests/reference_data/airway_counts.tsv")).unwrap();
    let (normalized_counts, normalized_gene_names, normalized_sample_names) =
        read_count_matrix(&normalized_path).unwrap();

    assert_eq!(normalized_gene_names, raw_gene_names);
    assert_eq!(normalized_sample_names, raw_sample_names);

    let r_normalized_path = Path::new("tests/reference_data/r_normalized_counts.tsv");
    let expected_counts = if r_normalized_path.exists() {
        let (r_counts, r_gene_names, r_sample_names) =
            read_count_matrix(r_normalized_path).unwrap();
        assert_eq!(r_gene_names, raw_gene_names);
        assert_eq!(r_sample_names, raw_sample_names);
        r_counts
    } else {
        let (r_sample_names, r_size_factors) =
            read_r_size_factors(Path::new("tests/reference_data/r_size_factors.tsv"));
        assert_eq!(r_sample_names, raw_sample_names);
        Mat::from_fn(raw_counts.nrows(), raw_counts.ncols(), |i, j| {
            raw_counts[(i, j)] / r_size_factors[j]
        })
    };

    assert_matrix_close(
        &normalized_counts,
        &expected_counts,
        &raw_gene_names,
        &raw_sample_names,
    );
}
