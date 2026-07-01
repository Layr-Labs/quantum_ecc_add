use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

const GOOD_OPS_FIXTURE: &str = "tests/testdata/good_current_ops.bin";
const BAD_OPS_FIXTURE: &str = "tests/testdata/bad_phase_classical_ancilla_ops.bin";
const BAD_ANCILLA_ONLY_OPS_FIXTURE: &str = "tests/testdata/bad_ancilla_only_ops.bin";

fn temp_eval_dir(name: &str) -> PathBuf {
    let mut p = std::env::temp_dir();
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    p.push(format!("quantum_ecc_{name}_{}_{}", std::process::id(), now));
    std::fs::create_dir_all(&p).unwrap();
    p
}

fn install_ops_fixture(dir: &Path, fixture: &str) {
    let src = Path::new(env!("CARGO_MANIFEST_DIR")).join(fixture);
    let dst = dir.join("ops.bin");
    if std::fs::hard_link(&src, &dst).is_err() {
        std::fs::copy(&src, &dst).unwrap();
    }
}

fn run_eval_circuit(dir: &Path, shots: usize) -> Output {
    Command::new(env!("CARGO_BIN_EXE_eval_circuit"))
        .current_dir(dir)
        .env("QECC_EVAL_SUPPRESS_RESULTS", "1")
        .env("QECC_EVAL_NUM_TESTS", shots.to_string())
        .arg("--sample-seed=op-bin-test-seed")
        .output()
        .unwrap()
}

#[test]
fn eval_circuit_accepts_good_prebuilt_circuit() {
    let dir = temp_eval_dir("good_ops_bin");
    install_ops_fixture(&dir, GOOD_OPS_FIXTURE);

    let output = run_eval_circuit(&dir, 8);
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(
        output.status.success(),
        "stdout:\n{stdout}\nstderr:\n{stderr}"
    );
    assert!(
        stdout.contains("loaded ops  : 28687418"),
        "stdout:\n{stdout}"
    );
    assert!(
        stdout.contains("statistical correctness OK"),
        "stdout:\n{stdout}"
    );
    assert!(
        stdout.contains("ancilla-garbage batches : 0"),
        "stdout:\n{stdout}"
    );

    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn eval_circuit_rejects_bad_prebuilt_circuit() {
    let dir = temp_eval_dir("bad_ops_bin");
    install_ops_fixture(&dir, BAD_OPS_FIXTURE);

    let output = run_eval_circuit(&dir, 1000);
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(!output.status.success());
    assert!(stdout.contains("loaded ops  : 1030"), "stdout:\n{stdout}");
    assert!(
        stdout.contains("classical mismatches    : 1000"),
        "stdout:\n{stdout}"
    );
    assert!(
        stdout.contains("phase-failed shots      : 1000"),
        "stdout:\n{stdout}"
    );
    assert!(
        stdout.contains("ancilla-garbage batches : 16"),
        "stdout:\n{stdout}"
    );
    assert!(
        stdout.contains("ANCILLA GARBAGE"),
        "stdout:\n{stdout}\nstderr:\n{stderr}"
    );

    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn eval_circuit_rejects_ancilla_only_failure() {
    let dir = temp_eval_dir("bad_ancilla_only_ops_bin");
    install_ops_fixture(&dir, BAD_ANCILLA_ONLY_OPS_FIXTURE);

    let output = run_eval_circuit(&dir, 8);
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(!output.status.success());
    assert!(
        stdout.contains("loaded ops  : 28687419"),
        "stdout:\n{stdout}"
    );
    assert!(
        stdout.contains("sampled failures        : 0"),
        "stdout:\n{stdout}"
    );
    assert!(
        stdout.contains("classical mismatches    : 0"),
        "stdout:\n{stdout}"
    );
    assert!(
        stdout.contains("phase-failed shots      : 0"),
        "stdout:\n{stdout}"
    );
    assert!(
        stdout.contains("phase-garbage batches   : 0"),
        "stdout:\n{stdout}"
    );
    assert!(
        stdout.contains("ancilla-garbage batches : 1"),
        "stdout:\n{stdout}"
    );
    assert!(
        stdout.contains("ANCILLA GARBAGE"),
        "stdout:\n{stdout}\nstderr:\n{stderr}"
    );

    let _ = std::fs::remove_dir_all(dir);
}
