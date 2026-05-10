# Prime Gap Depth — algorithm

Formal companion to [exploration.md](exploration.md). That document gives the
motivation, the empirical observations, and the open questions; this one
defines the construction precisely, states what we claim about it, and
matches the pseudocode to [main.rs:159](../src/main.rs#L159).

## 1. Setting

Let `S = (s_1, s_2, s_3, …)` be a strictly increasing sequence of positive
integers. The canonical case is `S = P`, the primes in ascending order, but
the algorithm makes sense for any such `S` and is implemented that way (see
`--seed-set` in the CLI).

For a finite prefix length `N ≥ 1`, write `S↾N = (s_1, …, s_N)`.

A **row** is any non-empty finite subsequence of `S` listed in ascending
order. We identify a row with the strictly increasing tuple of its
**positions** (indices into `S`), so a row `R = (i_1 < i_2 < … < i_k)`
represents the subsequence `(s_{i_1}, s_{i_2}, …, s_{i_k})`. Its
**leader** is `s_{i_1}`, the smallest element. Its **internal gaps** are

```
g_R(j) := s_{i_{j}} − s_{i_(j-1)},     for 2 ≤ j ≤ k.
```

These are gaps *within the row*, not gaps in `S`. They coincide with the
predecessor-gaps in `S` iff `R` consists of consecutive positions.

## 2. The construction

Fix a prefix length `N`. Define a sequence of multisets of rows
`L_0, L_1, L_2, …` by:

- **Base case.** `L_0 = { (1, 2, …, N) }` — a single row containing all
  positions in `S↾N`, in order.
- **Recursive step.** Given `L_ℓ`, form `L_{ℓ+1}` by replacing each row
  `R = (i_1, …, i_k) ∈ L_ℓ` with one child row per distinct value of
  `g_R(j)`, as follows. For each gap value `g > 0` such that
  `g_R(j) = g` for some `j`, the child row is

  ```
  R_g := (i_j : 2 ≤ j ≤ k, g_R(j) = g)
  ```

  i.e. the positions of those elements of `R` (other than the leader) that
  are preceded *within `R`* by gap `g`. The leader `i_1` does not appear in
  any child — `g_R(1)` is undefined, since `i_1` has no in-row predecessor,
  so there is no gap value to bucket it under. If `k = 1`, `R` has no
  children.

The recursion terminates: each row of length `k` produces children whose
total length is `k − 1`, so the total length across all rows at level `ℓ`
strictly decreases as `ℓ` increases (whenever any row has length ≥ 2).
Equivalently, every position appears in at most one row at each level, and
each position becomes a leader exactly once.

**Definition (depth, prefix version).** For each position `i ∈ {1, …, N}`,
let `m_N(s_i)` be the unique level `ℓ` such that `s_i` is the leader of
some row in `L_ℓ`. Equivalently, `m_N(s_i)` is the number of recursive
regroupings needed before `s_i` is the smallest element of its group.

## 3. Independence from the prefix length

**Theorem (well-definedness).** For any positions `i ≤ N ≤ N'`,

```
m_N(s_i) = m_{N'}(s_i).
```

Hence `m(s_i) := lim_{N → ∞} m_N(s_i)` is well-defined as a function
`m : S → ℕ`, and can be computed from any prefix that contains `s_i`.

*Proof sketch.* By induction on the level `ℓ`. The level-0 row in both
constructions is `(1, 2, …, N)` and `(1, 2, …, N')` respectively; their
internal gaps agree on the common prefix `(1, …, N)`. Hence the
buckets-by-gap-value at level 1 agree on positions `≤ N`: each bucket
present in the `N` construction is the truncation to positions `≤ N` of
the corresponding bucket in the `N'` construction (and the `N'`
construction may additionally have buckets for gap values that first
appear past position `N`, but those do not contain any position `≤ N`).

Within each bucket, the order is inherited from `S`, so the leader (the
position with smallest index) is identical in both constructions. The
inductive step is the same argument applied inside each bucket: the
sub-row in the `N` construction is a prefix of the sub-row in the `N'`
construction, so its internal gaps agree on the common prefix and the
bucketing-by-gap commutes with truncation.

The conclusion is that for every position `i ≤ N`, the row at level `ℓ`
containing `i` in the `N` construction is a prefix of the corresponding
row in the `N'` construction, and in particular has the same leader.
Therefore `i` becomes a leader at the same level in both constructions.
∎

The empirical sanity check is `pgd stability`
([main.rs:261](../src/main.rs#L261)), which computes `m_N` at
`N ∈ {10³, 10⁴, 10⁵, 10⁶}` and verifies that the values agree on the
common prefix.

## 4. The gap-path address

Each position `i` is the leader of exactly one row, reached by a unique
descent through the recursion tree. Record the gap value used at each
descent: this gives the **gap-path address**

```
addr(s_i) := (g_1, g_2, …, g_m)     where m = m(s_i).
```

Concretely, `g_ℓ` is the gap value that selects the bucket containing `i`
when the level-`(ℓ−1)` row containing `i` is split. Equivalently, if
`R^{(ℓ)}` is the level-`ℓ` row containing `i` and `R^{(ℓ)}` was produced
from `R^{(ℓ−1)}` as the bucket for gap `g_ℓ`, then `addr(s_i)` is
`(g_1, …, g_m)`.

The depth `m(s_i)` is the length of `addr(s_i)`. The address determines
the row containing `s_i` at every level, and conversely the leader of any
row uniquely determines that row's address. Two elements with the same
address at every level except the last share a row up to that level; this
is one candidate for the "group" axis discussed in
[exploration.md §"Where the periodic-table dream stands"](exploration.md).

The current implementation computes `m` only — the address itself is not
materialized.

## 5. Pseudocode

The reference implementation is [`compute_m` in main.rs:159](../src/main.rs#L159).
In pseudocode, with positions as `usize` indices into `S↾N`:

```
function compute_m(S):
    N := length(S)
    m := array of length N, initialized to ⊥
    queue := stack containing the single entry (level=0, row=[0, 1, …, N−1])

    while queue is non-empty:
        (level, row) := queue.pop()
        if row is empty: continue
        m[row[0]] := level                       # leader assignment
        if length(row) = 1: continue

        buckets := empty map gap_value → list of positions, keyed sorted ascending
        for j in 1 .. length(row) − 1:           # 0-indexed; corresponds to math index j' = j+1, 2..k
            g := S[row[j]] − S[row[j−1]]         # in-row gap g_R(j') = s_{i_{j'}} − s_{i_{j'-1}}
            buckets[g].append(row[j])            # bucket the destination i_{j'} of that gap

        for each (g, bucket) in buckets:
            queue.push((level + 1, bucket))

    return m
```

Notes on the implementation:

- Positions are stored, not values, so each row is a `Vec<usize>` rather
  than a `Vec<u64>`. This keeps every row a slice of indices into the
  original sequence and avoids copying values.
- A `BTreeMap<u64, Vec<usize>>` keys buckets by gap value; the order of
  iteration (smallest gap first) determines stack push order. Because the
  algorithm uses LIFO (`Vec::pop`), child rows for *larger* gaps are
  processed before child rows for *smaller* gaps. This affects only the
  order in which `m[·]` slots are written, not the values themselves.
- The recursion is structured as a worklist rather than actual recursion,
  so depth in `m` does not consume call stack.

## 6. Complexity

Let `N = |S↾N|` and `D = max m_N(s) + 1` be the depth observed at this
prefix. At each level `ℓ`, the rows in `L_ℓ` are pairwise disjoint and
their total size is at most `N`, so the work at level `ℓ` is
`O(N · log G_ℓ)` where `G_ℓ` is the maximum number of distinct gap
values within any single row at level `ℓ` (the `log` is the BTreeMap
insert cost). Bounding `G_ℓ ≤ N`, the total time is `O(N · D · log N)`.

In practice `D` grows extremely slowly — empirically `D = 7` at
`N = 10⁸` (i.e. max m = 6) — so the runtime is effectively linear up to
log factors. The N=10⁸ run took 65 s wall-clock with peak RSS 2.7 GB
on commodity hardware.

Space is `O(N)` for the output `m` array plus `O(N)` for all live rows
combined (each position appears in at most one row in the worklist at
any time). The 2.7 GB RSS at `N = 10⁸` is dominated by the input prime
array (`u64 × 10⁸ ≈ 800 MB`) plus the `Vec<usize>` row arrays.

## 7. Generalization beyond primes

The construction depends on `S` only through the partial order and the
gaps `s_j − s_{j-1}` (for `j ≥ 2`). Any strictly increasing sequence of
positive integers is a valid input. The implementation accepts arbitrary input
sequences via `--seed-set FILE` ([main.rs:1017](../src/main.rs#L1017)),
so the depth construction can be applied to:

- a sieved prime set (the default),
- a residue class of the primes (e.g. primes `≡ 1 (mod 4)`),
- a non-prime sequence with a chosen gap distribution (null-model
  control for whether the depth shape is intrinsic to the primes or
  generic to "irregular ascending integer sequences"),

and the same depth function `m` and address `addr` are defined on the
seed set.
