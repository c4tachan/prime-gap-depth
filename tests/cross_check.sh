#!/usr/bin/env bash
# Cross-check the Rust and PARI implementations of gap-depth.
#
# Runs:
#   1. Rust unit tests (cargo test).
#   2. PARI sanity tests (n=100 + first 50 of A395913).
#   3. Diff of m-values for the first 1000 primes between Rust and PARI.
#   4. Self-contained OEIS PROG snippet (regression check that what we paste
#      into OEIS still produces the correct first 50 terms of A395913).
#
# Exit non-zero on any failure.

set -euo pipefail

cd "$(dirname "$0")/.."

if ! command -v gp >/dev/null 2>&1; then
    echo "FAIL: gp (PARI/GP) not found on PATH; install with 'sudo apt install pari-gp'" >&2
    exit 1
fi

echo "=== 1. Rust unit tests ==="
cargo test --release --quiet

echo
echo "=== 2. PARI sanity tests ==="
gp -q pari/run_tests.gp

echo
echo "=== 3. Cross-check Rust vs PARI at N=1000 ==="
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

cargo run --release --quiet -- -n 1000 --outdir "$TMPDIR/rust" >/dev/null
awk -F, 'NR>1 {print $1","$2","$3}' "$TMPDIR/rust/results.csv" > "$TMPDIR/rust.csv"
gp -q pari/cross_check.gp > "$TMPDIR/pari.csv"

if diff -q "$TMPDIR/rust.csv" "$TMPDIR/pari.csv" >/dev/null; then
    echo "PASS: Rust and PARI m-values agree on first 1000 primes"
else
    echo "FAIL: m-values diverge"
    diff "$TMPDIR/rust.csv" "$TMPDIR/pari.csv" | head -20
    exit 1
fi

echo
echo "=== 4. OEIS PROG snippet regression ==="
EXPECTED="[19, 23, 43, 47, 67, 73, 101, 107, 109, 131, 139, 151, 163, 173, 179, 181, 193, 227, 229, 233, 241, 263, 269, 277, 293, 311, 317, 337, 353, 379, 383, 389, 401, 433, 449, 463, 467, 487, 491, 499, 503, 509, 563, 569, 577, 593, 599, 601, 613, 619]"
ACTUAL=$(gp -q pari/oeis_prog.gp)
if [[ "$ACTUAL" == "$EXPECTED" ]]; then
    echo "PASS: pari/oeis_prog.gp produces the OEIS A395913 draft's first 50 terms"
else
    echo "FAIL: oeis_prog.gp output does not match the OEIS draft"
    echo "  expected: $EXPECTED"
    echo "  actual:   $ACTUAL"
    exit 1
fi

echo
echo "All cross-checks passed."
