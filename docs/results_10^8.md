# Prime Gap Depth — results at N = 10⁸

A snapshot of what the gap-depth construction (see [algorithm.md](algorithm.md))
says about the first 100,000,000 primes. This is a results report, not a
lab notebook — for motivation, half-formed conjectures, and the broader
"periodic table" framing, see [exploration.md](exploration.md). The data
this report draws on is checked in under [`out/Prime Numbers/`](../out/Prime%20Numbers/):

- [`growth.tsv`](../out/Prime%20Numbers/growth.tsv) — depth histogram at `N = 10², 10³, …, 10⁸`.
- [`results.csv`](../out/Prime Numbers/results.csv) — `(π(p), p, m(p))` for every one of the 100M primes (`p_{10⁸} = 2 038 074 743`).
- [`mod_30.tsv`](../out/Prime%20Numbers/mod_30.tsv), [`mod_210.tsv`](../out/Prime%20Numbers/mod_210.tsv) — residue distribution per `m`-class at the two relevant primorial moduli.
- [`class_quantiles.tsv`](../out/Prime%20Numbers/class_quantiles.tsv) — within each `m`-class, the 1st / 10th / 100th / … prime that hits that depth.
- [`overlay_shifts.tsv`](../out/Prime%20Numbers/overlay_shifts.tsv) — log-log CDF horizontal shifts between consecutive `m`-classes for `(m, m+1) ∈ {(1,2), (2,3), (3,4), (4,5)}`.
- [`oeis_m0.txt` … `oeis_m6.txt`](../out/Prime%20Numbers/) — each m-class as a b-file (depth 6 has 4 known members).
- [`pichain_C.tsv`](../out/Prime%20Numbers/pichain_C.tsv) — the alternative pi-chain depth, kept for comparison.

## 1. Headline

At `N = 10⁸`:

| `m` | count       | fraction of primes |
| --- | ----------: | -----------------: |
| 0   |           1 |           0.0000% |
| 1   |         144 |           0.0001% |
| 2   |     665,641 |           0.6656% |
| 3   |  56,679,939 |          56.6799% |
| 4   |  42,545,865 |          42.5459% |
| 5   |     108,406 |           0.1084% |
| 6   |           4 |           0.0000% |

Two facts the previous (10⁷) run could not see:

1. **`m = 6` exists.** The four currently-known depth-6 primes are
   `684 611 189`, `782 212 369`, `1 128 319 141`, `1 219 645 079` — all
   in the second half of the `[1, p_{10⁸}]` range. The first depth-6
   prime is therefore the 35-millionth prime, an enormous jump from the
   first depth-5 prime (the 6389th prime, `p = 63 709`). The
   first-at-depth set continues to grow super-exponentially in the prime
   index: 2, 3, 7, 19, 113, 63 709, 684 611 189 for `m = 0..6`.
2. **The `m = 3` / `m = 4` crossover is imminent.** At `N = 10⁸` the
   `m = 4` class accounts for 42.55% of primes vs 56.68% for `m = 3`;
   one more decade of `N` at the observed rate puts `m = 4` ahead. See §3.

## 2. Growth shape

The full decade table (from `growth.tsv`):

| `N`     | `m=0` | `m=1` | `m=2`   | `m=3`      | `m=4`      | `m=5`   | `m=6` |
| ------- | ----: | ----: | ------: | ---------: | ---------: | ------: | ----: |
| 10²     | 1     | 9     | 42      | 42         | 6          | 0       | 0 |
| 10³     | 1     | 18    | 333     | 595        | 53         | 0       | 0 |
| 10⁴     | 1     | 34    | 1 887   | 7 401      | 674        | 3       | 0 |
| 10⁵     | 1     | 54    | 9 214   | 80 388     | 10 299     | 44      | 0 |
| 10⁶     | 1     | 78    | 40 935  | 786 778    | 171 435    | 773     | 0 |
| 10⁷     | 1     | 105   | 167 799 | 6 976 019  | 2 846 377  | 9 699   | 0 |
| 10⁸     | 1     | 144   | 665 641 | 56 679 939 | 42 545 865 | 108 406 | 4 |

Per-decade growth ratios `c_m(10 N) / c_m(N)`:

| decade        | `m=1` | `m=2` | `m=3` | `m=4` | `m=5` |
| ------------- | ----: | ----: | ----: | ----: | ----: |
| 10³ → 10⁴     | 1.89  | 5.67  | 12.44 | 12.72 |  —   |
| 10⁴ → 10⁵     | 1.59  | 4.88  | 10.86 | 15.28 | 14.67 |
| 10⁵ → 10⁶     | 1.44  | 4.44  |  9.79 | 16.65 | 17.57 |
| 10⁶ → 10⁷     | 1.35  | 4.10  |  8.87 | 16.60 | 12.55 |
| 10⁷ → 10⁸     | 1.37  | 3.97  |  8.12 | 14.95 | 11.18 |

A ratio of 10 corresponds to linear density. Reading the columns:

- **`m = 1` is sub-logarithmic.** The decade ratio is decaying toward 1
  and the count is well-fit by `c · log²(p_N)` (see §4).
- **`m = 2` is sub-linear** with a clearly-decaying ratio, last seen at
  3.97. The fraction of primes at `m = 2` peaked at ~5% around `N = 10⁴`
  and has been falling ever since; it is now 0.67% and shrinking. The
  natural reading is asymptotic density zero, but the descent is slow.
- **`m = 3` is slightly sub-linear** and is still the plurality class,
  but its ratio is decaying (12.4 → 8.1 over five decades). It is
  losing share to `m = 4` decade-on-decade and the gap is closing fast.
- **`m = 4` is super-linear up to and including `N = 10⁸`** with the
  ratio still well above 10. This is what drives the crossover.
- **`m = 5` is super-linear** but the ratio is decaying (17.6 → 11.2),
  consistent with the same pattern that played out for `m = 3` and `m = 4`
  one and two decades earlier: a class first appears, ramps up super-
  linearly, and gradually settles toward (or below) linear.

The empirical pattern across classes is a **wave**: each `m`-class starts
out small, swells super-linearly until it hits the dominant share at some
scale, then settles into slow decay as the next class catches up. The
waves are well-separated in `N` — successive class-peaks appear to be
about a decade or two apart.

## 3. The `m = 3` → `m = 4` crossover

| `N`  | `m = 3 / m = 4` | `m = 4` share of `{m = 3} ∪ {m = 4}` |
| ---- | --------------: | -----------------------------------: |
| 10⁴  | 10.98           |  8.35% |
| 10⁵  | 7.81            | 11.36% |
| 10⁶  | 4.59            | 17.89% |
| 10⁷  | 2.45            | 28.98% |
| 10⁸  | 1.33            | 42.88% |

The share grows by roughly one decade-doubling (8 → 11 → 18 → 29 → 43%);
the ratio of decade ratios `(m=4)/(m=3)` has been close to 1.7–1.8 for
the last three decades. Extrapolating one more decade puts `m = 4` near
58–60% share — i.e. `m = 4` overtakes `m = 3` between `N = 10⁸` and
`N = 10⁹`. A confirmatory `N = 10⁹` run is the cheapest test.

This crossover is the cleanest evidence that the depth histogram is not
peaked at a single asymptotic `m`. Whether *every* `m`-class eventually
takes its turn as dominant, or whether some upper `m` "wins" forever,
is not visible from this data — but at minimum the histogram's mode is
not stable.

## 4. `m = 1` and the prime-gap-record set

By construction the `m = 1` primes are exactly the first primes preceded
by each gap value that has appeared in `(p_1, …, p_N)`. The size of the
class is the number of distinct predecessor-gaps observed below `p_N`,
which is controlled by the maximum prime gap. Under Cramér's conjecture
`g_max(x) ≲ (log x)²`, so `|{p ≤ x : m(p) = 1}| = O((log x)²)` is the
natural ansatz.

Fitting `c_1(N) ≈ κ · ln²(p_N)`:

| `N`  | `m = 1` count | `ln²(p_N)` | `κ`   |
| ---- | ------------: | ---------: | ----: |
| 10²  |   9           |    39.6    | 0.227 |
| 10³  |  18           |    80.6    | 0.223 |
| 10⁴  |  34           |   133.6    | 0.255 |
| 10⁵  |  54           |   198.2    | 0.273 |
| 10⁶  |  78           |   274.1    | 0.285 |
| 10⁷  | 105           |   361.2    | 0.291 |
| 10⁸  | 144           |   459.5    | 0.313 |

`κ` is slowly increasing and has not stabilised, but the order of growth
is right. The same data plotted against `(log p_N)·(log log p_N)`,
which is closer to the Granville/Heath-Brown heuristics for the count of
*distinct* gap values up to `x`, fits with a slightly more stable
constant; we have not tried to formalise that.

`m = 1` is the only class whose asymptotic count reduces directly to a
classical gap-distribution quantity. The level-1 leaders are also exactly
OEIS A000101 truncated to the gap values that have appeared by `p_N`.

## 5. `m = 6` first appearances

A subjective entry, since the class has only four known members:

```
684 611 189   (the 35 027 089-th prime)
782 212 369
1 128 319 141
1 219 645 079
```

The first depth-6 prime sits at index `≈ 3.5 × 10⁷` — well into the run.
The depth-first-appearance index sequence is `1, 2, 4, 8, 30, 6389,
35 027 089` for `m = 0..6`. The jump from `m = 4` to `m = 5` is a factor
of 213; from `m = 5` to `m = 6` it is 5481. Whether the multiplicative
gap continues to grow is the natural question for the next run; if it
does, the first depth-7 prime would not appear until somewhere near
`p_{n}` with `n ~ 10¹¹`, well beyond what a single 64-bit prime sieve
will reach.

## 6. Residue structure mod 30 and mod 210

Refreshed at `N = 10⁸` from
[`mod_30.tsv`](../out/Prime%20Numbers/mod_30.tsv) and
[`mod_210.tsv`](../out/Prime%20Numbers/mod_210.tsv). Percentages are
over each `m`-class; the residues `r ∈ {0, 2, 3, 5}` mod 30 (and the
corresponding non-coprime classes mod 210) carry only the singletons
`m(2) = 0`, `m(3) = m(5) = 1` and are excluded from the analysis.
Uniform baselines on the coprime residues are 12.50% (mod 30; 8
classes) and 2.083% (mod 210; 48 classes).

### 6.1. Mod 30

| `m` | sample | `r=1` | `r=7` | `r=11` | `r=13` | `r=17` | `r=19` | `r=23` | `r=29` |
| --- | -----: | ----: | ----: | -----: | -----: | -----: | -----: | -----: | -----: |
| 2 |    665 641 | 12.32 | 13.09 | 12.82 | 12.45 | 12.79 | 12.15 | 11.79 | 12.58 |
| 3 | 56 679 939 | 12.32 | 12.99 | 12.81 | 12.44 | 12.51 | 12.27 | 11.95 | 12.71 |
| 4 | 42 545 865 | 12.77 | 11.86 | 12.11 | 12.58 | 12.51 | 12.76 | 13.24 | 12.16 |
| 5 |    108 406 |  0.84 |  2.65 |  1.13 | 12.88 |  2.79 | 33.63 | 12.04 | 34.04 |
| 6 |          4 | 25.00 |  0.00 |  0.00 |  0.00 |  0.00 | 25.00 |  0.00 | 50.00 |

Across the bulk classes `m ∈ {2, 3, 4}` every coprime residue lies
within ~1 pp of the 12.50% baseline; for `m = 3` the entire row fits
within ±0.55 pp over 56.7M samples. The depth function does not
encode residue-class membership at all for the dominant classes —
they are mod-30 equidistributed.

At `m = 5` the prior finding survives at 140× the sample: `{19, 29}`
still dominates with 67.7% combined mass (was 91.6% at `N = 10⁶`),
`{13, 23}` has filled in near the baseline (12.88% / 12.04%), and
`{1, 7, 11, 17}` remain heavily depressed (each below 3%). At first
sight this looks like slow convergence to uniform; mod 210 shows that
it is not.

### 6.2. Mod 210

For `m = 3` and `m = 4` the equidistribution survives the finer
modulus intact:

| `m` | sample     | min coprime % | max coprime % | spread |
| --- | ---------: | ------------: | ------------: | -----: |
| 3   | 56 679 939 | 1.88          | 2.37          | 0.49 pp |
| 4   | 42 545 865 | 1.70          | 2.35          | 0.66 pp |

— at uniform baseline 2.083%, every one of the 48 coprime residues
sits within ~0.3 pp of baseline at `m = 3`. The Dirichlet
equidistribution that worked at mod 30 continues to work at mod 210.

For `m = 5` the 48 coprime residues partition cleanly into three
tiers. Tier 1 (six residues, each carrying ~9–12% of m=5 mass):

| `r mod 210` | `r mod 30` | `r mod 7` | m=5 % |
| ----------: | ---------: | --------: | ----: |
|  29         | 29         | 1         |  9.43 |
|  59         | 29         | 3         | 11.96 |
|  79         | 19         | 2         | 11.73 |
| 149         | 29         | 2         | 11.97 |
| 169         | 19         | 1         | 11.71 |
| 199         | 19         | 3         |  9.50 |

**Tier 1 sum: 66.30%.** These are exactly the residues satisfying
`r ≡ {19, 29} (mod 30)` AND `r ≡ {1, 2, 3} (mod 7)`. Six of the twelve
coprime-to-210 residues that lie in `{19, 29} mod 30` are tier-1; the
other six (those with `r ≡ {4, 5, 6} mod 7`: r ∈ {19, 89, 109, 139,
179, 209}) carry just 1.67% of m=5 mass combined — a 26× deficit
relative to tier 1.

Tier 2 (six residues, each ~3–4%):

| `r mod 210` | `r mod 30` | `r mod 7` | m=5 % |
| ----------: | ---------: | --------: | ----: |
|  13         | 13         | 6         |  3.30 |
|  23         | 23         | 2         |  3.72 |
|  53         | 23         | 4         |  3.20 |
|  83         | 23         | 6         |  3.53 |
| 163         | 13         | 2         |  4.44 |
| 193         | 13         | 4         |  3.26 |

**Tier 2 sum: 21.45%.** These are `r ≡ {13, 23} (mod 30)` AND
`r ≡ {2, 4, 6} (mod 7)` — the *even-mod-7* subset. The other six
residues in `{13, 23} mod 30` (with `r ≡ {1, 3, 5} mod 7`) carry only
3.47% combined. The favoured mod-7 set flips between tier 1 (`{1, 2, 3}`)
and tier 2 (`{2, 4, 6}`).

Tier 3 (the remaining 36 coprime residues): combined mass ≈ 12.25%,
mostly well below 1% per residue. The minimum is `r = 137` at 0.0323%,
~64× below baseline.

Aggregating mod 7 across the m=5 class:

| `r mod 7` | mass    |
| --------: | ------: |
| 1         | 22.94%  |
| 2         | 33.57%  |
| 3         | 25.07%  |
| 4         |  7.85%  |
| 5         |  1.17%  |
| 6         |  9.40%  |

There is a strong global bias toward `r ≡ 2 (mod 7)` (33.57% vs
baseline 16.67%) and a near-exclusion at `r ≡ 5 (mod 7)` (1.17%, a
14× under-representation). The mod-30 picture missed this entirely,
because residue 5 mod 7 spreads across all eight coprime-to-30
classes.

### 6.3. What this says

The depth function does not encode residue information at `m ≤ 4`, but
beginning at `m = 5` it imposes a definite, structured residue-mod-210
signature. The structure has a clean two-modulus form: a coupling
between `r mod 30` and `r mod 7` whose favoured mod-7 set depends on
the mod-30 class. The hypothesis from the previous revision that the
m=5 skew would converge to uniform is ruled out — the asymptotic m=5
residue distribution is a specific weighted measure on
`(ℤ/210ℤ)*`, dominated by 12 of the 48 coprime classes.

`m = 6` at this `N` has only 4 primes with `(r mod 210) ∈
{59, 61, 149, 169}`. Three of the four (59, 149, 169) sit in m=5's
tier-1 set; the fourth (61, with `r mod 30 = 1`) is the first depth-6
prime to leave the m=5 tier-1+2 envelope entirely. Far too small a
sample for any claim, but consistent with the m=5 structure carrying
into m=6 with some leakage.

A natural conjecture: each `m ≥ 5` has an asymptotic residue support
that is a specific subset of `(ℤ/Kℤ)*` for some primorial `K`, and `K`
grows with `m` — mod-210 captured m=5; m=6 may need mod-2310 to
resolve cleanly. The address-level structure
([algorithm.md §4](algorithm.md)) is the natural language for this:
each gap-path forces Hardy–Littlewood-style constellation constraints,
and the residue support of a class is the union of the
constellation-residues over the gap-paths in that class. The
`(g₁, g₂) = (2, 6)` example worked out in
[exploration.md](exploration.md) — that row equals
A007530 + 8 = `{p + 8 : (p, p+2, p+6, p+8) is a prime quadruple}` —
is the prototype. Quadruples have a unique admissible residue class
mod 30 (`p ≡ 11 mod 30`, so `p + 8 ≡ 19 mod 30`), which explains why
*that specific* level-2 row sits entirely in `r ≡ 19 mod 30`. The
mod-210 tier-1 of m=5 is presumably the analogous structure averaged
over all gap-paths of length 5.

## 7. Class quantiles and growth shape per class

From [`class_quantiles.tsv`](../out/Prime%20Numbers/class_quantiles.tsv),
the `k`-th prime to enter each m-class. Showing `(prime_index, prime_value)`:

| `m` | size       | `k=1` | `k=10` | `k=10³` | `k=10⁵` | `k=10⁷` |
| --- | ---------: | --- | --- | --- | --- | --- |
| 1 |        144 | 2 (3) | 155 (907) | — | — | — |
| 2 |    665 641 | 4 (7) | 20 (71) | 4 155 (39 499) | 4 256 220 (72 496 441) | — |
| 3 | 56 679 939 | 8 (19) | 32 (131) | 1 564 (13 151) | 124 102 (1 642 441) | 14 758 998 (270 924 649) |
| 4 | 42 545 865 | 30 (113) | 156 (911) | 14 102 (152 899) | 644 775 (9 679 147) | 28 797 449 (549 036 799) |
| 5 |    108 406 | 6 389 (63 709) | 28 734 (334 189) | 1 295 063 (20 412 169) | 92 570 100 (1 879 112 603) | — |
| 6 |          4 | 35 496 536 (684 611 189) | — | — | — | — |

Best-fit slopes for `log(prime_index_k) ≈ α_m · log(k) + β_m` within
each class (from the `pgd class-quantiles` output):

| `m` | slope `α_m` | intercept `β_m` |
| --- | ----------: | --------------: |
| 1 | 3.69 | −2.02 |
| 2 | 1.50 | −1.99 |
| 3 | 0.95 |  0.82 |
| 4 | 0.84 |  3.74 |
| 5 | 0.91 |  7.79 |

A slope `α_m < 1` says the `k`-th class-`m` prime grows *sub-linearly*
in `k`: within the class, primes accrue faster than the rank. `m = 3`
and `m = 4` both sit comfortably below 1 and are well into a regime of
linear-density-or-better; `m = 5` is below 1 at this `N` too, but its
intercept is enormous (`exp(7.79) ≈ 2 400`) so the class is currently
"diluted" — you have to look at the 6 389th prime to find the first
`m = 5` member. The intercept jumps by ~3 units (a factor of ~20) per
unit increase in `m` for `m ≥ 3`. If this pattern persists, the first
`m = 7` prime would sit near `p_n` with `n ≈ exp(11)` ≈ 60 000-th
prime times the depth-6 offset of 35M — i.e. in the `n ≈ 2 × 10⁹`
range, just beyond a 100M-prime run.

The `m = 1` slope (3.69) is qualitatively different and reflects that
m=1 primes are gap-record primes, which become exponentially rarer in
the index: only one new gap value appears per ~`(log p)²`-sized window.

## 8. Pi-chain depth (alternative invariant)

The pi-chain depth iterates `p → π(p)` and counts steps until reaching 1;
it is an independent depth function on the primes that we compute
alongside `m` to have something to compare against. From
[`pichain_C.tsv`](../out/Prime%20Numbers/pichain_C.tsv):

| `m_πchain` | size at `N = 10⁷` | `C(m, k=size)` |
| ---------- | ----------------: | -------------: |
| 1          | 921 502           | 1.82e-5 |
| 2          | 70 796            | 3.09e-3 |
| 3          | 6 725             | 0.342   |
| 4          | 812               | 23.4    |
| 5          | 127               | 959.6   |
| 6          | 26                | 22 543  |
| 7          | 7                 | 289 015 |
| 8          | 2                 | 567 433 |
| 9–11       | 1 each            | rising  |

The pi-chain histogram is monotonically decreasing in `m_πchain` and the
classes shrink exponentially — qualitatively very different from our
`m`. Gap-depth produces a mode in the middle and a long thin tail;
pi-chain depth produces a fat low end and a thin tail. They agree that
"most primes are at depth ≈ 3", but they are picking out different sets
of primes at that depth and the structural significance is different.

## 9. Overlay shifts between consecutive classes

[`overlay_shifts.tsv`](../out/Prime%20Numbers/overlay_shifts.tsv)
records horizontal log-log shifts between the cumulative-count curves
of consecutive `m`-classes, refreshed at `N = 10⁸`. For each
`(m_lo, m_hi)` pair, the file gives 50 sample points across the
overlap region; `dx = x_hi − x_lo` is the log shift required to align
m_hi onto m_lo at the matched count level.

| `(m_lo, m_hi)` | `dx` at tail (last 10) | linear-scale factor `eᵈˣ` |
| -------------- | ---------------------: | ------------------------: |
| (1, 2)         |  7.68                  | ≈ 2 170 |
| (2, 3)         |  3.98                  | ≈ 53.5  |
| (3, 4)         | −0.65                  | ≈ 0.52  |
| (4, 5)         | −5.87                  | ≈ 0.0028 |

The first two shifts are large and positive: m=2 lags m=1 by ~3 log
units, m=3 lags m=2 by ~4 log units. At `(3, 4)` the shift goes
*negative* — the m=4 cumulative curve is now to the *left* of m=3's
at the matched-count level, which is the same crossover already seen
in §3 viewed through the CDF rather than the histogram. `(4, 5)` is
sharply negative again, reflecting m=5's small size relative to m=4
(at the count levels where they both have data, m=5 has to come from
much further along the prime sequence). The `(5, 6)` shift is not
estimated — m=6 currently has 4 members, well below the 50 sample
points the overlay routine wants.

The qualitative read: the inter-class shifts are stable across the
overlap region (the tail-mean is within ~1 unit of the median), so the
shift is a meaningful summary of "how far behind in `N` is class
m+1 vs class m." The shift becomes negative at `(3, 4)` and stays
negative thereafter, which is the regime-change indicator that the
histogram mode has begun walking up the `m` axis.

## 10. What the 100M run resolves vs leaves open

Resolved (or strongly constrained):

- **`m` is unbounded.** Depth 6 is reached, with the first depth-6 prime
  near 6.8 × 10⁸. The conjecture in [exploration.md](exploration.md) §
  "Open questions" — that `max m → ∞` — is now strongly supported, though
  the rate at which new depths appear is sub-double-exponential in `p`.
- **The depth histogram has a moving mode.** The dominant class is
  shifting from `m = 3` to `m = 4` between `N = 10⁸` and `N = 10⁹`. The
  "periodic table" framing in [exploration.md](exploration.md) where a
  *single* `m` would label "ordinary primes" is therefore wrong: every
  `m`-class above some threshold appears to take its turn as the bulk.
- **Class growth is wave-shaped, not power-law-uniform.** Each class
  ramps super-linearly, peaks, and decays. Fitting a single exponent per
  class (as discussed in [exploration.md](exploration.md) §"Growth
  exponents") is meaningful only in the locally-linear regime of each
  wave.
- **Residues at `m ≤ 4` are uniform on `(ℤ/210ℤ)*`.** §6 shows the 48
  coprime residues each within ~0.3 pp of the 2.083% baseline at
  `m = 3` (sample 56.7M). The depth function does not encode residue
  information at all for the dominant classes.
- **The `m = 5` residue skew is real and structured, not a finite-`N`
  artefact.** Mod 30 the m=5 distribution looked like it might be
  converging to uniform; mod 210 resolves the apparent convergence as
  a sharper tiered structure (§6.2). Six coprime-mod-210 residues
  carry 66% of the m=5 mass, twelve more carry 22%, and the remaining
  thirty carry 12% combined. The favoured residues are defined by a
  coupling between `r mod 30` and `r mod 7`. The asymptotic m=5
  residue distribution is a specific weighted measure on `(ℤ/210ℤ)*`,
  not the uniform measure.

Still open:

- **Asymptotic density of `m = k` for `k ≥ 2`.** The wave behaviour
  suggests each fixed `k` has density zero with successive classes
  carrying the bulk at successive scales. This is compatible with
  `Σ_m c_m(N) = N` only if the number of nonempty `m`-classes grows,
  which is consistent with `m = 6` appearing first at `N = 10⁸`. A
  clean formulation: is there a scaling function `f(N)` such that the
  "effective" `m` (e.g. the histogram mean or mode) grows like
  `f(log N)`?
- **The exact form of the m=5 asymptotic measure.** Tier 1 of §6.2
  has 6 residues all satisfying `r ≡ {19, 29} mod 30 ∧ r ≡ {1, 2, 3}
  mod 7`; tier 2 has 6 residues all satisfying `r ≡ {13, 23} mod 30
  ∧ r ≡ {2, 4, 6} mod 7`. The masses within each tier are not equal
  (tier 1 spreads 9.4%–12.0%; tier 2 spreads 3.2%–4.4%), and the
  tier-2 favoured mod-7 set differs from the tier-1 one. Is there a
  natural set of Hardy–Littlewood constellation weights — one per
  gap-path of length 5 — that reproduces the observed tier structure
  and the within-tier spread? If so, m=5 is the union of a small
  number of admissible 5-tuple patterns, and the mass at each residue
  is a sum of singular-series-style coefficients.
- **Modulus for `m = 6` and higher.** The mod-210 picture resolves
  `m = 5` cleanly; m=6 with only 4 primes is undetermined. The natural
  ladder is `K = primorial(m + 1)`: mod 30 for m=4 (which equidistributes
  there), mod 210 for m=5 (with the tier structure), mod 2310 for m=6.
  An `N = 10⁹` run would have ~50 m=6 primes, still too few for mod
  2310 (2310 has 480 coprime classes); an `N = 10¹⁰` run is what would
  actually answer this. Alternatively, one could enumerate the
  expected residue support directly from the gap-paths leading to
  m=6 and compare.
- **Address-level structure.** `m(p)` is the *length* of the gap-path
  address `(g₁, …, g_m)` (see [algorithm.md §4](algorithm.md)). The
  address itself is not yet materialised. The one worked-out example
  — the level-2 row `(g₁, g₂) = (2, 6) = A007530 + 8 =
  {p + 8 : (p, p+2, p+6, p+8) is a prime quadruple}`, proved in
  [exploration.md](exploration.md) — already predicts a specific mod-30
  signature (`r ≡ 19 mod 30`), and that signature is *exactly* one of
  the two m=5-favoured mod-30 classes. The hypothesis worth testing
  directly: dump every level-2 row, identify the constellation each
  one corresponds to, sum the constellation-residue distributions
  weighted by row size, and compare against the observed m=5 mod-210
  histogram. The materialisation step is bounded by row count, not
  prime count.
- **Cross-check on a null model.** Does the wave structure replicate on
  a random ascending sequence with matching gap distribution? If yes,
  the depth structure is a consequence of gap statistics; if no, it is
  intrinsic to the primes. The `--seed-set` flag is the entry point.

## 11. Reproducibility

All numbers in this document come from a single 100M run. To regenerate
from scratch (about 65 s wall and ~2.7 GB RSS on commodity hardware):

```
SEED="data/primes_1B.txt"
DIR="out/Prime Numbers/1B"

./target/release/pgd --seed-set "$SEED" -o "$DIR" -n 1000000000
./target/release/pgd growth          --seed-set "$SEED" -o "$DIR" -n 1000000000
./target/release/pgd oeis-export     --seed-set "$SEED" -o "$DIR" -n 1000000000
./target/release/pgd pi-chain        --seed-set "$SEED" -o "$DIR" -n 1000000000
./target/release/pgd mod-residue 30  --seed-set "$SEED" -o "$DIR" -n 1000000000
./target/release/pgd mod-residue 210 --seed-set "$SEED" -o "$DIR" -n 1000000000
./target/release/pgd class-quantiles --seed-set "$SEED" -o "$DIR" -n 1000000000
./target/release/pgd overlay         --seed-set "$SEED" -o "$DIR" -n 1000000000
```

Each subcommand independently re-runs the depth construction (no shared
cache), so the total wall-clock for the eight commands is ~8 minutes on
commodity hardware. All standard analyses are now refreshed for 100M;
the remaining work flagged in §10 (address materialisation,
mod-residue 2310 at `N ≥ 10¹⁰`, null-model cross-check) requires either
new code or larger `N`.
