use assert_cmd::Command;
use predicates::prelude::*;

#[test]
fn test_version() {
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("--version")
        .assert()
        .success()
        .stdout(predicate::str::contains(env!("CARGO_PKG_VERSION")));
}

#[test]
fn test_no_args() {
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("Usage"));
}
