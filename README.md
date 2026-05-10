# prime-gap-depth

A Rust CLI (`pgd`) that computes a depth function `m(p)` on the primes via
iterated regrouping by gap value, and provides several analyses of the
resulting classes.

For motivation and empirical observations, see [docs/exploration.md](docs/exploration.md).
For the formal construction, theorem, and complexity analysis, see
[docs/algorithm.md](docs/algorithm.md). This README is the operational
overview.

## Algorithm

Let `S = (s_1, s_2, …)` be a strictly increasing sequence of positive
integers (the primes by default). Build a sequence of multisets of rows
`L_0, L_1, …`:

- **Level 0:** a single row `(s_1, s_2, …, s_N)` containing the entire
  prefix in order.
- **Recursive step:** for each row `R` at the current level, compute the
  internal gaps `g_j = R[j+1] − R[j]` and bucket the non-leader entries by
  the value of `g_j`. Each bucket (sorted ascending) becomes a child row at
  the next level. The first element of `R` is its **leader**.

Every element of `S` is the leader of exactly one row, at exactly one
level. Define

```
m(s_i) := the level at which s_i is a leader
```

i.e. the number of gap-regroupings needed before `s_i` is the smallest
element of its group. Empirically `m` grows extremely slowly on the
primes — depth 6 is the maximum observed at `N = 10⁸`.

A theorem (proved in [docs/algorithm.md §3](docs/algorithm.md)) shows
`m(s_i)` is independent of the prefix length `N` used to compute it, so
`m : S → ℕ` is well-defined as a function of the full sequence. The
`stability` subcommand checks this empirically.

The reference implementation is `compute_m` in [src/depth.rs](src/depth.rs);
it runs the recursion as a worklist over `Vec<usize>` index-rows with
`BTreeMap<u64, Vec<usize>>` for bucketing, giving `O(N · D · log N)`
runtime where `D = max m + 1`. A second implementation in PARI/GP lives in
[pari/gap_depth.gp](pari/gap_depth.gp) — used as an independent
cross-check and as the source of OEIS `PROG` entries.

## Building

Requires a stable Rust toolchain.

```
cargo build --release
```

The binary is produced at `target/release/pgd`.

## Running

The default invocation computes `m` (and the alternative pi-chain depth)
on the first 1,000,000 primes, prints histograms, and writes a CSV:

```
./target/release/pgd
```

Common flags (apply to most subcommands):

| flag | description |
| --- | --- |
| `-n, --count N` | number of elements to use (default `1_000_000`) |
| `--seed-set FILE` | use the ascending integers in `FILE` (one per line) instead of primes |
| `-o, --outdir DIR` | output directory for CSV/TSV (default `out/`) |

### Subcommands

| command | what it does |
| --- | --- |
| *(none)* | compute `m` and pi-chain depth, print histograms, write `out/results.csv` |
| `stability` | recompute `m` at `N ∈ {10³, 10⁴, 10⁵, 10⁶}` and verify agreement on the common prefix |
| `mod-residue [MOD]` | distribution of primes mod `MOD` per `m`-class with chi-squared p-values (default mod 30) |
| `growth` | how `m`-class counts grow with `N` from 100 to 100M |
| `oeis-export` | export each small `m`-class as an OEIS b-file (`out/oeis_m*.txt`) |
| `first-at [MAX_M]` | find the first prime to achieve each `m`-value up to `MAX_M` (default 6) |
| `class-quantiles` | within each `m`-class, the 1st, 10th, 100th, … prime that hits that level |
| `overlay` | log-log CDFs of `m`-classes; estimate horizontal shift between consecutive classes |
| `predict [--m-min M] [--m-max M]` | fit intercept(`m`) for converged classes and project forward |
| `pi-chain` | pi-chain depth: family counts, first appearances, `C(m,k)`, ratios |

### Examples

Compute on 10M primes, write outputs to `./results/`:

```
./target/release/pgd -n 10000000 -o results
```

Distribution of `m`-classes mod 4:

```
./target/release/pgd mod-residue 4
```

Run the construction on a custom ascending integer sequence:

```
./target/release/pgd --seed-set my_sequence.txt
```

Verify `m` is prefix-independent:

```
./target/release/pgd stability
```

## Output

The default run writes `out/results.csv` with one row per input number,
containing the number, its `m` value, and (when running on primes) its
pi-chain depth. Subcommands write additional TSV/CSV files into the same
directory; see the `out/` folder for examples produced by recent runs.

## Testing

Rust unit tests pin the n=100 histogram and the m-class membership for
m ∈ {0, 1, 4} against hard-coded values:

```
cargo test --release
```

For end-to-end cross-validation against the independent PARI/GP
implementation, run:

```
bash tests/cross_check.sh
```

This runs four layers:

1. `cargo test` (Rust correctness against pinned values).
2. PARI sanity tests at n=100 (same pinned values, plus first 50 of A395913).
3. Diff of `(pi, p, m)` for the first 1000 primes between the Rust binary
   and the PARI implementation — must be byte-identical.
4. Regression check that [pari/oeis_prog.gp](pari/oeis_prog.gp) — the
   self-contained snippet pasted into OEIS submissions — still produces
   the exact first 50 terms of A395913.

Requires `pari-gp` on `PATH` (`sudo apt install pari-gp` on Debian/Ubuntu).

## OEIS submissions

The b-files in `out/oeis_m*.txt` are produced by `pgd oeis-export`. The
PARI `PROG` snippet for any gap-depth sequence is in
[pari/oeis_prog.gp](pari/oeis_prog.gp); the body of that file (between
the function definition and the trailing `print(...)` line) is what gets
pasted into the OEIS PROG entry, with the final `print` line adjusted for
the target depth and term count.

## Performance

The construction is effectively linear up to log factors. As a reference
point, `N = 10⁸` runs in ~65 s wall-clock with peak RSS ~2.7 GB on
commodity hardware (the input prime array dominates memory).
