#!/usr/bin/env python3
"""Generate the first N lucky numbers (OEIS A000959) and write them
one-per-line to a file. Lucky numbers are produced by the sieve:

    1. Start from positive integers 1, 2, 3, 4, 5, ...
    2. The 2nd surviving number is 2 -> remove every 2nd surviving number.
    3. The next surviving k is 3 -> remove every 3rd surviving number.
    4. Continue: each step uses the next surviving number not yet used as
       a sieve key, removing every k-th remaining entry.

Usage: gen_lucky.py N OUTFILE
"""
import sys
import numpy as np


def lucky_numbers(n_target: int) -> np.ndarray:
    # Asymptotic L(n) ~ n log n; pad generously.
    upper = max(200, int(n_target * (np.log(max(n_target, 2)) + np.log(np.log(max(n_target, 3)))) * 1.5) + 100)
    print(f"  initial upper bound: {upper:,}", file=sys.stderr)
    # Start with odd numbers (post k=2 sieve step).
    arr = np.arange(1, upper + 1, 2, dtype=np.int64)
    print(f"  after k=2: {arr.size:,} entries", file=sys.stderr)
    # i is the 0-based index of the next sieve key to apply.
    i = 1
    while i < arr.size:
        k = int(arr[i])
        if k > arr.size:
            break
        # Drop every k-th element (1-indexed positions). Use slice deletion
        # via boolean mask but apply in-place via np.delete for clarity.
        # Indices to delete: k-1, 2k-1, 3k-1, ... < arr.size
        del_idx = np.arange(k - 1, arr.size, k, dtype=np.int64)
        arr = np.delete(arr, del_idx)
        i += 1
        if i % 50 == 0:
            print(f"  step i={i}, k={k}, arr.size={arr.size:,}", file=sys.stderr)
    print(f"  done sieving after {i} steps; arr.size={arr.size:,}", file=sys.stderr)
    return arr[:n_target]


def main() -> None:
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    n = int(sys.argv[1])
    out = sys.argv[2]
    print(f"Sieving for first {n:,} lucky numbers...", file=sys.stderr)
    lucky = lucky_numbers(n)
    print(f"Got {lucky.size:,} lucky numbers; max = {int(lucky[-1]):,}", file=sys.stderr)
    if lucky.size < n:
        print(f"WARNING: only produced {lucky.size}, less than requested {n}", file=sys.stderr)
    with open(out, "w") as f:
        for v in lucky:
            f.write(f"{int(v)}\n")
    print(f"Wrote {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
