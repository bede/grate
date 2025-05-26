use assert_cmd::Command;
use std::fs;
use std::path::Path;
use tempfile::tempdir;

// Create one of several test fastas
fn create_test_fasta(path: &Path, variant: usize) {
    let fasta_content = match variant {
        1 => ">seq1\nACGTACGTACGT\n>seq2\nCGTACGTACGTA\n",
        2 => ">seq1\nTGCATGCATGCA\n>seq2\nGCATGCATGCAT\n",
        _ => ">seq1\nACGTACGTACGT\n>seq2\nGTACGTACGTAC\n",
    };
    fs::write(path, fasta_content).unwrap();
}

// Index builder helper
fn build_index(fasta_path: &Path, bin_path: &Path) {
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("index")
        .arg("build")
        .arg(fasta_path)
        .arg("-o")
        .arg(bin_path)
        .assert()
        .success();

    // Check file exists and isn't empty
    assert!(
        bin_path.exists(),
        "Index file wasn't created at {:?}",
        bin_path
    );
    assert!(
        fs::metadata(bin_path).unwrap().len() > 0,
        "Index file is empty"
    );
}

#[test]
fn test_index_build() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("test.fasta");
    let bin_path = temp_dir.path().join("test.bin");

    create_test_fasta(&fasta_path, 1);

    // Build index and save to file using -o
    build_index(&fasta_path, &bin_path);
}

#[test]
fn test_index_build_with_custom_kmer_window() {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join("test.fasta");
    let bin_path = temp_dir.path().join("test.bin");

    create_test_fasta(&fasta_path, 1);

    // Build index with custom k-mer length and window size using -o
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("index")
        .arg("build")
        .arg(fasta_path)
        .arg("-k")
        .arg("15")
        .arg("-w")
        .arg("10")
        .arg("-o")
        .arg(&bin_path)
        .assert()
        .success();

    // Check file exists and isn't empty
    assert!(bin_path.exists());
    assert!(fs::metadata(&bin_path).unwrap().len() > 0);
}

#[test]
fn test_index_union() {
    let temp_dir = tempdir().unwrap();
    let fasta1_path = temp_dir.path().join("test1.fasta");
    let fasta2_path = temp_dir.path().join("test2.fasta");
    let bin1_path = temp_dir.path().join("test1.bin");
    let bin2_path = temp_dir.path().join("test2.bin");
    let combined_path = temp_dir.path().join("combined.bin");

    // Create different test FASTA files
    create_test_fasta(&fasta1_path, 1);
    create_test_fasta(&fasta2_path, 2);

    // Build indexes
    build_index(&fasta1_path, &bin1_path);
    build_index(&fasta2_path, &bin2_path);

    // Combine indexes
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("index")
        .arg("union")
        .arg("-o")
        .arg(&combined_path)
        .arg(&bin1_path)
        .arg(&bin2_path)
        .assert()
        .success();

    // Check combined file exists
    assert!(combined_path.exists());

    // The combined size, should be larger than either index
    let combined_size = fs::metadata(&combined_path).unwrap().len();
    let bin1_size = fs::metadata(&bin1_path).unwrap().len();
    let bin2_size = fs::metadata(&bin2_path).unwrap().len();

    let max_individual_size = std::cmp::max(bin1_size, bin2_size);
    assert!(
        combined_size >= max_individual_size,
        "Combined index size {} should be at least as large as the largest individual index size {}",
        combined_size,
        max_individual_size
    );
}

#[test]
fn test_index_diff() {
    let temp_dir = tempdir().unwrap();
    let fasta1_path = temp_dir.path().join("test1.fasta");
    let fasta2_path = temp_dir.path().join("test2.fasta");
    let bin1_path = temp_dir.path().join("test1.bin");
    let bin2_path = temp_dir.path().join("test2.bin");
    let result_path = temp_dir.path().join("result.bin");

    // Create test FASTAs
    create_test_fasta(&fasta1_path, 1);
    create_test_fasta(&fasta2_path, 2);

    // Build indexes
    build_index(&fasta1_path, &bin1_path);
    build_index(&fasta2_path, &bin2_path);

    // Diff second index from first
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("index")
        .arg("diff")
        .arg("-o")
        .arg(&result_path)
        .arg(&bin1_path)
        .arg(&bin2_path)
        .assert()
        .success();

    // Check diffed file exists
    assert!(result_path.exists());

    // The diffed should be smaller than or equal to the first index
    let result_size = fs::metadata(&result_path).unwrap().len();
    let bin1_size = fs::metadata(&bin1_path).unwrap().len();

    assert!(
        result_size <= bin1_size,
        "Result index size {} should be less than or equal to the first index size {}",
        result_size,
        bin1_size
    );
}
