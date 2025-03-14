use assert_cmd::Command;
use std::fs;
use std::path::Path;
use std::process::Command as StdCommand;
use tempfile::tempdir;

// Helper to create a small FASTA file for testing
fn create_test_fasta(path: &Path) {
    let fasta_content = ">seq1\nACGTACGTACGT\n>seq2\nGTACGTACGTAC\n";
    fs::write(path, fasta_content).unwrap();
}

// Helper to create a small FASTQ file for testing
fn create_test_fastq(path: &Path) {
    let fastq_content =
        "@seq1\nACGTACGTACGT\n+\n~~~~~~~~~~~~\n@seq2\nGTACGTACGTAC\n+\n~~~~~~~~~~~~\n";
    fs::write(path, fastq_content).unwrap();
}

// Helper to run index build and capture output to a file
fn build_index(fasta_path: &Path, bin_path: &Path) {
    let output = StdCommand::new(assert_cmd::cargo::cargo_bin("deacon"))
        .arg("index")
        .arg("build")
        .arg(fasta_path)
        .output()
        .expect("Failed to execute command");

    fs::write(bin_path, output.stdout).expect("Failed to write index file");
    assert!(output.status.success(), "Index build command failed");
}

#[test]
fn test_filter_to_file() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq");
    let report_path = temp_dir.path().join("report.json");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);

    build_index(&fasta_path, &bin_path);
    assert!(bin_path.exists(), "Index file wasn't created");

    // Run filtering command with writing raw (dog) fastq file
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .arg("--report")
        .arg(&report_path)
        .assert()
        .success();

    // Check output and report creation
    assert!(output_path.exists(), "Output file wasn't created");
    assert!(report_path.exists(), "Report file wasn't created");

    // Validate output content
    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(!output_content.is_empty(), "Output file is empty");
}

#[test]
fn test_filter_to_file_gzip() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq.gz");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    // Check gzipped output file creation
    assert!(output_path.exists(), "Gzipped output file wasn't created");
    assert!(
        fs::metadata(&output_path).unwrap().len() > 0,
        "Gzipped output file is empty"
    );
}

#[test]
fn test_filter_to_file_zstd() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered.fastq.zst");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    // Check that zstd output file was created
    assert!(output_path.exists(), "Zstd output file wasn't created");
    assert!(
        fs::metadata(&output_path).unwrap().len() > 0,
        "Zstd output file is empty"
    );
}

#[test]
fn test_filter_invert() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_inverted.fastq");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--invert")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    assert!(
        output_path.exists(),
        "Output file with inverted flag wasn't created"
    );
}

#[test]
fn test_filter_rename() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_renamed.fastq");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--rename")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    assert!(
        output_path.exists(),
        "Output file with rename flag wasn't created"
    );

    let output_content = fs::read_to_string(&output_path).unwrap();
    assert!(
        output_content.contains("@1\n") || output_content.contains("@2\n"),
        "Output does not contain renamed sequences"
    );
}

#[test]
fn test_filter_min_matches() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_min_matches.fastq");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--matches")
        .arg("2")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    assert!(
        output_path.exists(),
        "Output file with min_matches parameter wasn't created"
    );
}

#[test]
fn test_filter_prefix_length() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("ref.fasta");
    let fastq_path = temp_dir.path().join("reads.fastq");
    let bin_path = temp_dir.path().join("ref.bin");
    let output_path = temp_dir.path().join("filtered_prefix.fastq");

    create_test_fasta(&fasta_path);
    create_test_fastq(&fastq_path);
    build_index(&fasta_path, &bin_path);

    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("filter")
        .arg("--nucleotides")
        .arg("6")
        .arg(&bin_path)
        .arg(&fastq_path)
        .arg("--output")
        .arg(&output_path)
        .assert()
        .success();

    assert!(
        output_path.exists(),
        "Output file with prefix_length parameter wasn't created"
    );
}
