use deseq2_rs::io::{read_count_matrix, read_coldata};

#[test]
fn test_read_count_matrix_small() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("counts.tsv");
    std::fs::write(&path, "gene\tS1\tS2\nGENE1\t10\t20\nGENE2\t30\t40\n").unwrap();

    let (counts, gene_names, sample_names) = read_count_matrix(&path).unwrap();
    assert_eq!(gene_names, vec!["GENE1", "GENE2"]);
    assert_eq!(sample_names, vec!["S1", "S2"]);
    assert_eq!(counts.nrows(), 2);
    assert_eq!(counts.ncols(), 2);
    assert_eq!(counts[(0, 0)], 10.0);
    assert_eq!(counts[(1, 1)], 40.0);
}

#[test]
fn test_read_coldata_small() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("coldata.tsv");
    std::fs::write(&path, "\tdex\nSRR1039508\tuntrt\nSRR1039509\ttrt\n").unwrap();

    let (df, sample_names) = read_coldata(&path).unwrap();
    assert_eq!(sample_names, vec!["SRR1039508", "SRR1039509"]);
    assert_eq!(df.get_column("dex").unwrap(), &["untrt", "trt"]);
}
