use assert_cmd::Command;
use predicates::str;

#[test]
fn test_version() {
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.arg("--version")
        .assert()
        .success()
        .stdout(str::contains(env!("CARGO_PKG_VERSION")));
}

#[test]
fn test_no_args() {
    let mut cmd = Command::cargo_bin("deacon").unwrap();
    cmd.assert().failure().stderr(str::contains("Usage"));
}
