# Prime Gap Depth — working notes

Last updated: 2026-05-10

## Where this came from

I started out wanting a "periodic table" of the primes: some discrete
classification that groups primes into a small number of structurally meaningful
buckets — analogous to periods/groups for the elements — where membership in a
bucket would predict shared behaviour (residues, gap statistics, distribution,
something).

The natural starting point is the gap sequence. Twin primes, cousin primes,
sexy primes, etc. already sort primes by predecessor-gap. That gives countably
many classes but the classes are unbounded and the structure is flat — not
really a periodic table, just a partition by a single feature.

So I tried iterating: partition by gap, then partition each part by its own
internal gap, recursively. The thing that fell out is not a periodic table,
but a depth function `m(p)` that turns out to be interesting in its own right
and is what this repo computes. These notes document the construction, what I
think I'm seeing, and what I don't yet understand.

## The construction

Let `P = (p_1, p_2, p_3, …)` be the primes in ascending order. Define rows
recursively, starting with a single level-0 row containing all of `P`.

For each row `R` at the current level:

1. The first element of `R` is the **leader** and is assigned `m(leader) = current level`.
2. Compute consecutive differences inside `R`: `d_i = R[i] - R[i-1]` for `i ≥ 1`.
3. Bucket the non-leader entries `R[1], R[2], …` by the value of `d_i`.
   Each bucket, sorted ascending, becomes a new row at level `current + 1`.

Repeat until no rows remain. Each prime is the leader of exactly one row, at
exactly one level, so `m : P → ℕ` is total.

In words: `m(p)` is the number of gap-regroupings it takes before `p` is the
smallest element of its group.

### Worked example, first few levels

Level 0: `[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, …]`. Leader is 2, so `m(2) = 0`.
Gaps: `1, 2, 2, 4, 2, 4, 2, 4, 6, …`.

Level 1 rows (one per gap value seen):

- gap 1: `[3]` — leader 3, `m(3) = 1`.
- gap 2: `[5, 7, 13, 19, 31, 43, 61, 73, 103, 109, …]` — leader 5, `m(5) = 1`.
- gap 4: `[11, 17, 23, 41, 47, 53, …]` — leader 11, `m(11) = 1`.
- gap 6: `[29, 37, 53, …]` — leader 29, `m(29) = 1`.
- gap 8: `[97, …]` — leader 97, `m(97) = 1`.
- …

The level-1 leaders are exactly the **first prime to be preceded by gap g**,
for each gap value `g` that occurs at least once. For `N = 100` these are
`[3, 5, 11, 29, 97, 127, 149, 211, 541]`, which matches OEIS A000101 (the
first prime *after* a gap of `2n`) for the gap values that have appeared by
the 100th prime.

Inside a level-1 row, take the gap-2 row `[5, 7, 13, 19, 31, …]`. Its internal
gaps are `2, 6, 6, 12, …`. So at level 2 it spawns:

- "preceded-by-2 then internal-gap-2": `[7]` — `m(7) = 2`.
- "preceded-by-2 then internal-gap-6": `[13, 19, 109, …]` — `m(13) = 2`.
- "preceded-by-2 then internal-gap-12": `[31, 43, 73, …]` — `m(31) = 2`.

And so on.

## Is `m(p)` well-defined independently of `N`?

Empirically, yes — that's what the `stability` subcommand checks. For the
first 1000 primes, computing `m` with cutoffs `N = 1k, 10k, 100k, 1M` gives
the same `m`-value for every prime.

The reason is that the level-0 row is just an ascending prefix of `P`, so the
level-1 buckets are determined by the gaps appearing in that prefix. Adding
more primes can only add new entries to the *tail* of existing buckets or
create new buckets for previously-unseen gap values; it cannot change which
prime is the leader of an existing bucket, nor can it change the order of
primes already in a bucket. By induction this propagates through all levels.

So `m(p)` is a property of `p` alone (assuming the construction is run on a
prefix that includes `p`). This is what makes it worth studying — it's a
genuine arithmetic-flavoured invariant, not an artefact of where we truncated.

## Empirical observations

All from `out/`, computed against the first 1M primes (and 10M for the growth
table).

### Histogram is sharply peaked

| N        | m=0 | m=1 | m=2    | m=3      | m=4      | m=5  |
| -------- | --- | --- | ------ | -------- | -------- | ---- |
| 100      | 1   | 9   | 42     | 42       | 6        | 0    |
| 1,000    | 1   | 18  | 333    | 595      | 53       | 0    |
| 10,000   | 1   | 34  | 1,887  | 7,401    | 674      | 3    |
| 100,000  | 1   | 54  | 9,214  | 80,388   | 10,299   | 44   |
| 1,000,000 | 1  | 78  | 40,935 | 786,778  | 171,435  | 773  |
| 10,000,000 | 1 | 105 | 167,799 | 6,976,019 | 2,846,377 | 9,699 |

The bulk of primes sit at `m = 3`. `m = 0` is always a singleton (just 2),
and `m = 1` grows extremely slowly — clearly sub-logarithmic in `N`.

### `m = 1` is asymptotically the prime-gap-record set

The `m = 1` primes are exactly the first-occurrence primes for each
predecessor-gap that has appeared so far. The number of such primes up to `N`
is the number of distinct gaps observed in `(p_1, …, p_N)`, which is
controlled by the **maximum prime gap** below `p_N`. Cramér's conjecture puts
this at `O((log p_N)^2)`, which would make `|{p ≤ x : m(p) = 1}| = O((log x)^2)`.
The observed counts (9, 18, 34, 54, 78, 105) are consistent with that order of
magnitude though I haven't fit a curve.

### Growth exponents (rough, eyeballed from the table above)

- `m = 2` looks sublinear: ratio per decade of `N` drops from ~5.7 to ~4.4 to
  ~4.1. Possibly `N^α` with `α ≈ 0.6`.
- `m = 3` looks roughly linear. It's the dominant class throughout, and the
  ratio per decade has settled near 9 — i.e. growing slightly slower than `N`.
- `m = 4` is currently *superlinear* in this range (171k → 2.85M is a 16.6×
  jump for 10×N). Either there's a regime change coming, or the dominant class
  shifts from `m = 3` to `m = 4` somewhere ahead.
- `m = 5` first appears at `N = 10⁴` and is exploding: 3 → 44 → 773 → 9,699.

The natural conjecture is that for each fixed `m ≥ 2` the class has positive
density asymptotically, but the densities are reached at very different scales
of `N`, and finite `N` always shows a peak at one particular `m`. The depth
itself is presumably unbounded (`max m` should keep climbing as `N` grows),
but I don't have a proof or even a strong heuristic.

### `m = 5` has a striking residue skew mod 30

From `out/mod_30.tsv`, residues are percentages within each m-class:

| m | r=1   | r=7   | r=11  | r=13  | r=17  | r=19  | r=23  | r=29  |
| - | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- |
| 2 | 12.27 | 13.75 | 12.88 | 12.15 | 13.08 | 11.78 | 11.46 | 12.63 |
| 3 | 12.58 | 13.05 | 12.62 | 12.11 | 12.39 | 12.34 | 12.08 | 12.82 |
| 4 | 12.14 |  9.75 | 11.95 | 14.42 | 12.89 | 13.20 | 14.75 | 10.90 |
| 5 |  0.39 |  0.13 |  0.65 |  3.88 |  0.39 | 47.22 |  2.98 | 44.37 |

(Other residues — 0, 2, 3, 5 — are essentially zero for large primes, as
expected.)

For `m = 2, 3` the distribution is close to uniform across the eight residues
coprime to 30 (12.5% each). For `m = 4` there's a mild skew toward 13, 19, 23.
For **`m = 5` the distribution is wildly non-uniform**: roughly 92% of the
773 m=5 primes are at residue 19 or 29 mod 30 — i.e. `≡ −11` or `≡ −1` mod 30.

I don't have an explanation. Possible directions:

- This could be a finite-`N` artefact: with only 773 primes the distribution
  hasn't equidistributed yet, and a small constellation of "early" m=5 primes
  with a particular shape is dragging the distribution. The same effect should
  appear at lower m at smaller `N`. A scan of `m = 4` at `N = 10⁴` (only 674
  primes) would be a good control.
- It could be real: the regrouping process at depth 5 might preserve a residue
  bias inherited from the gap structure. Gaps of 2, 6, 12, 30 are all `≡ 0 mod
  6`, and chained gap constraints can force residue classes — see Hardy–
  Littlewood prime constellation densities.
- Worth re-running mod-30 at `N = 10M` once `m = 5` has ~10× more members and
  seeing whether the skew softens or hardens.

The chi-squared p-values are all reported as 0 — that's the large-N statistical
significance pathology, not effect size. What matters is the absolute deviation
from 12.5%.

## Connections to classical prime constellations

An OEIS reviewer asked whether the iteration-2 construction is "just the prime
quadruples." It isn't — but the question turns out to be sharper than it looks,
and the answer reveals a clean correspondence worth recording.

### One specific level-2 row equals A007530 + 8

At level 2, the row indexed by `(predecessor-gap = 2, within-row gap = 6)` is
exactly the set of largest elements of prime quadruples. That is,

```
{ q ∈ row (2, 6) at level 2 } = { p + 8 : (p, p+2, p+6, p+8) is a prime quadruple } = A007530 + 8.
```

The row begins `[13, 19, 109, 199, 829, 1489, 1879, 2089, …]`, which is OEIS
A007530 = `[5, 11, 101, 191, 821, 1481, 1871, 2081, …]` shifted by 8.

*Proof sketch.* In any prime quadruple `(p, p+2, p+6, p+8)` with `p > 3`, `p+4`
is divisible by 3 (since `p` and `p+2` cover the two nonzero residues mod 3,
forcing `p+4 ≡ 0 mod 3`). So `p+4` is composite, which means `p+2` and `p+8`
are *consecutive* in the level-1 gap-2 row (no other prime between them has
predecessor-gap 2). Their within-row gap is `(p+8) − (p+2) = 6`. The converse
runs the same argument backwards: any element `q` of the level-2 `(2, 6)` row
forces `(q−8, q−6, q−2, q)` to be a prime quadruple.

### Iteration 2 is broader than this single row

The reviewer's question generalises to "is *iteration 2* the prime quadruples?"
— and the answer is no, because iteration 2 produces many rows, of which
`(2, 6)` is one. Other early rows at level 2:

| `(g, h)` | First few elements | Connection to a constellation? |
| --- | --- | --- |
| `(2, 2)` | `[7]` | trivial singleton |
| `(2, 6)` | `[13, 19, 109, 199, 829, …]` | **A007530 + 8 (prime quadruples)** |
| `(2, 12)` | `[31, 43, 73, 313, …]` | not a known constellation I've identified |
| `(4, 6)` | `[17, 23, 47, 53, …]` | not a known constellation I've identified |
| `(4, 18)` | `[41, 191, …]` | not identified |

So 31 is at iteration 2 but is not in any prime quadruple, which kills the
"iteration 2 = quadruples" reading.

### What this opens up

The interesting follow-up: **which level-2 (and level-3+) rows correspond to
known prime k-tuple constellations, and which are genuinely new?**

A prime k-tuple constellation is a fixed admissible pattern
`(0, b_1, b_2, …, b_{k−1})` that primes can simultaneously fit. Each
constellation imposes a specific gap-pattern on consecutive primes that fit it.
Conjecturally, the iterated-gap construction's rows at sufficient depth should
encode all admissible constellations, with deeper rows encoding longer
patterns. Some plausible identifications to check:

- `(2, 6)` ≅ A007530 + 8 (prime quadruples) ✓ proved above.
- `(4, 6)` and `(4, 18)` involve sequences of primes preceded by gap 4 with
  specific within-row spacings — likely correspond to admissible 4- or 5-tuple
  patterns starting with a gap of 4. I haven't identified them.
- Level-3 rows would correspond to longer constellations (5- or 6-tuples).
  The first few `m = 3` primes (19, 23, 43, 47, 67, 73, …) are the leaders of
  level-3 rows; checking which OEIS k-tuple sequences contain these would be a
  cheap next step.

This is the cleanest "structure" result the construction has produced so far:
the row-coordinate `(g_1, g_2, …, g_m)` is **not arbitrary** — it encodes
constellation-membership in at least some cases. If this generalises, the
gap-path address discussed in [§"Where the periodic-table dream stands"](#where-the-periodic-table-dream-stands)
might literally be a coordinate system on the prime constellations.

## Where the periodic-table dream stands

Still open — but the shape it would take is different from what I first
imagined.

What I originally wanted was a small fixed number of classes with shared
chemical-property-like behaviour. The naive reading of `m(p)` doesn't deliver
that: the number of non-empty m-classes grows with `N`, and the bulk class
shifts as `N` grows. So `m` alone is not "the period."

But there are several places periodicity could still be hiding, and I haven't
ruled them out:

1. **The full address, not just the depth.** Each prime is the leader of
   exactly one row, and that row has a canonical address: the gap-path
   `(g_1, g_2, …, g_m)` you take through the level-1, level-2, … buckets to
   reach it. `m` is just the length of that address. The full address is
   much closer to a periodic-table coordinate — depth could play the role of
   "period" and the gap-path the role of "group". Two primes that share the
   same gap-path tail might be analogous to elements in the same group: same
   outer structure, different scale. I haven't computed the full address yet,
   only `m`. This is the next thing to try.
2. **Residue periodicity within m-classes.** The mod-30 table shows that
   higher m-classes have non-uniform residues mod 30, with m=5 concentrating
   sharply on `{19, 29}`. Residues mod a primorial are inherently periodic
   structure. If the residue signature of class `m = k` stabilises as `N→∞`
   to a specific subset of `(ℤ/30ℤ)*` (or `(ℤ/210ℤ)*`), that *is* a periodic
   classification — primes are partitioned by depth into residue-flavoured
   groups, with the groups being literal subgroups/cosets of the unit group
   mod a primorial.
3. **Periodic structure in `m(p_n)` as a function of `n`.** Treating `m` as
   a sequence indexed by prime index, are there autocorrelations, beat
   patterns, or structure visible in a Fourier transform? I haven't looked.
   A quick `np.fft` on `m(p_1), m(p_2), …` would either show something or
   rule it out cheaply.
4. **The construction on a primorial-aligned starting set.** Starting from
   the primes coprime to a primorial (say everything coprime to 30 below
   `x`) instead of from all primes might expose a period-30 structure that
   the all-primes construction smears out.

So the working hypothesis isn't "no periodic table" — it's "the periodic
table, if it exists, is in `(m, gap-path)` or in `(m, residue mod
primorial)`, not in `m` alone." The depth function may turn out to be one
axis of the table rather than the whole thing.

What's already true and worth noting: the m-classes are stratified — a thin
top layer (m=0, m=1) of canonical "early" primes, a small intermediate layer
(m=2), a fat middle (m=3, m=4), and a long tail — and the layers have
different residue signatures and very different growth rates, both unexpected
for a construction defined purely in terms of gap differences. Whatever this
turns into, the layers are real.

## Open questions

1. **Is `m` unbounded?** Conjecturally yes; the max-m for `N = 10⁷` is 5,
   so the growth is slow. A sieve up to `N = 10⁹` or `10¹⁰` would be the test.
2. **Does each m-class have positive natural density?** Or do all classes
   above some threshold have density zero, with a single "winning" class
   asymptotically?
3. **Why the residue skew at m=5?** Finite-N artefact or persistent feature?
   See controls above.
4. **Connection to gap-distribution heuristics.** Can the asymptotics of
   `|{p ≤ x : m(p) = k}|` be derived from Cramér / Hardy–Littlewood
   conjectures about the gap distribution? The level-1 case reduces cleanly
   (it's the prime-gap-record set), so maybe level-2 does too.
5. **Construction on non-prime sequences.** The `--seed-set FILE` flag exists
   precisely to test whether the depth-distribution shape is intrinsic to the
   primes or generic to "sufficiently irregular ascending integer sequences."
   Random sequences with matching gap distribution would be the cleanest
   control.
6. **Is there a periodic table hiding in the `(m, gap-path)` coordinate?**
   Compute the full gap-path address for each prime, not just its length.
   Group primes by tail-of-path or by path-shape (multiset of gaps along the
   path) and look for shared residue / density / distribution properties.
   This is the most direct attack on the original goal.
7. **Does residue-class membership stabilise per m-class?** As `N → ∞`, does
   `{r ∈ (ℤ/30ℤ)* : Pr[p ≡ r | m(p) = k]}` converge to a *subset* of the
   coprime residues (rather than the full set with uniform measure)? If so,
   that subset is the "group" for that period.
8. **Spectral structure in `m` as a sequence.** FFT of
   `(m(p_1), m(p_2), …, m(p_N))` — periodic peaks would be a smoking gun.

## Software map

The CLI subcommands map directly to the questions above:

- `pgd` (no subcommand) — compute `m(p)` for the first `N` primes and dump
  histogram + CSV.
- `pgd stability` — sanity check that `m` is dataset-independent.
  Question 0 (does the construction even define a function on primes?).
- `pgd mod-residue 30 | 210` — residue distributions per m-class.
  Question 3.
- `pgd growth` — class counts at `N = 10², …, 10⁸`.
  Questions 1 and 2.
- `pgd oeis-export` — small classes as OEIS b-files.
  Submission / external corroboration.
- `pgd first-at MAX_M` — smallest prime achieving each depth.
  Question 1 (the depth-record set).

## Things to try next

- **Compute and dump the full gap-path address per prime**, not just `m`.
  This unlocks the periodic-table-as-`(m, path)` line of attack and is
  probably the highest-information next step relative to its cost.
- Run `growth` out to `N = 10⁸` and re-fit the per-class growth exponents.
- Re-run `mod-residue 30` at `N = 10⁷` and `10⁸` and watch whether the m=5
  skew persists or softens. If it persists, check whether each m-class has a
  *stable* coprime-residue support set.
- `mod-residue 210` for finer residue structure (each gap-2 chain forces a
  mod-2,3,5,7 constraint).
- FFT of `m(p_n)` indexed by `n` — a cheap test for spectral periodicity.
- Run on a `--seed-set` of "random ascending integers with matching gap
  distribution" as a null model. If the depth structure replicates, it's
  generic to gap statistics, not specific to primes.
- Run on a `--seed-set` of "primes coprime to a primorial" or "primes in a
  fixed residue class mod 30" to expose any primorial-aligned periodicity.
- Identify A-numbers for the m-class b-files in `out/` and check whether any
  are already in OEIS. The m=1 class is A000101; the m=2 class onward I
  haven't checked.
- **Identify which level-2 and level-3 rows correspond to known OEIS
  k-tuple sequences.** The `(g=2, h=6)` row equals A007530 + 8 (see
  "Connections to classical prime constellations" above); systematically
  dump each row and search OEIS for the first few elements. Each match is
  a confirmed constellation/row correspondence; non-matches may be new
  sequences or admissible patterns whose primes have not been catalogued.
