"""
convergence_rate.py

Measures how fast m-class primes converge to their PNT-like asymptote
(p_m(k) ~ C * k^alpha * ln k) as k grows, and compares to the known
convergence rate of the all-primes PNT (p(n) ~ n*ln n).

Method:
  1. Load plot_data.tsv
  2. For each m-class, fit the PNT model in log space
  3. Bin data into log-spaced k-bins, compute RMS |residual| per bin
  4. Fit a power law to RMS vs k (slope = convergence exponent)
  5. Overlay the all-primes theoretical correction term ln(ln n)/ln n

Usage:
  python scripts/convergence_rate.py ".\out\Prime Numbers\10^9\plot_data.tsv"
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── Load data ────────────────────────────────────────────────────────────────

path = sys.argv[1] if len(sys.argv) > 1 else r".\out\Prime Numbers\10^9\plot_data.tsv"

print(f"Loading {path} ...")
data_by_m: dict[int, tuple[np.ndarray, np.ndarray]] = {}

with open(path, "r") as fh:
    header = fh.readline()
    k_lists: dict[int, list] = {}
    p_lists: dict[int, list] = {}
    for line in fh:
        parts = line.split("\t")
        m, approx_rank, class_size, prime_index, prime_value = (
            int(parts[0]), int(parts[1]), int(parts[2]),
            int(parts[3]), int(parts[4]),
        )
        k_lists.setdefault(m, []).append(approx_rank)
        p_lists.setdefault(m, []).append(prime_value)

for m in k_lists:
    k = np.array(k_lists[m], dtype=np.float64)
    p = np.array(p_lists[m], dtype=np.float64)
    order = np.argsort(k)
    data_by_m[m] = (k[order], p[order])

print(f"  Loaded {len(data_by_m)} m-classes: {sorted(data_by_m)}")

# ── Fit PNT model per m-class ─────────────────────────────────────────────────
# Model: ln p = alpha*ln k + ln(ln k) + ln C
# i.e.   ln p - ln(ln k) = alpha * ln k + ln C
# Linear regression in (ln k, ln p - ln(ln k))

MIN_K = 20  # exclude very small k
N_BINS = 12  # log-spaced bins across the k range

COLORS = {2:'#e41a1c', 3:'#ff7f00', 4:'#4daf4a', 5:'#377eb8',
          6:'#984ea3', 7:'#a65628', 8:'#f781bf'}

fits: dict[int, dict] = {}

for m, (k_all, p_all) in sorted(data_by_m.items()):
    mask = k_all >= MIN_K
    k, p = k_all[mask], p_all[mask]
    if len(k) < 50:
        continue

    lk   = np.log(k)
    lp   = np.log(p)
    llk  = np.log(lk)           # ln(ln k) — small correction
    y    = lp - llk             # = alpha*ln k + ln C  (ideally)

    A    = np.column_stack([lk, np.ones_like(lk)])
    sol, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    alpha, lnC = sol
    C    = math.exp(lnC)

    # Residuals in log space
    resid = lp - (alpha*lk + llk + lnC)   # ln(actual/fitted)

    fits[m] = dict(alpha=alpha, C=C, lnC=lnC, k=k, p=p,
                   resid=resid, lk=lk)
    print(f"  m={m}: alpha={alpha:.4f}, C={C:.4f},  "
          f"RMS_resid={math.sqrt(np.mean(resid**2)):.4f}")

# ── Bin residuals and measure convergence rate ────────────────────────────────

print("\nConvergence rate (power law fit to RMS-residual vs k):")
print(f"  {'m':>2}  {'slope β':>8}  {'=> RMS ~ k^(-β)':>18}  half-decades to halve RMS")

fig = plt.figure(figsize=(14, 10))
gs  = gridspec.GridSpec(2, 1, hspace=0.35)
ax1 = fig.add_subplot(gs[0])   # RMS residual vs k
ax2 = fig.add_subplot(gs[1])   # all-primes comparison

conv_results = []

for m, fd in sorted(fits.items()):
    k, resid, lk = fd['k'], fd['resid'], fd['lk']

    # Log-spaced bins by k
    k_min, k_max = k.min(), k.max()
    edges = np.logspace(np.log10(k_min), np.log10(k_max), N_BINS + 1)
    bin_k   = []
    bin_rms = []
    for i in range(N_BINS):
        mask = (k >= edges[i]) & (k < edges[i+1])
        if mask.sum() < 5:
            continue
        r = resid[mask]
        bin_k.append(np.sqrt(edges[i] * edges[i+1]))   # geometric midpoint
        bin_rms.append(math.sqrt(np.mean(r**2)))

    if len(bin_k) < 4:
        continue

    bk  = np.array(bin_k)
    brms= np.array(bin_rms)

    # Fit log(RMS) = -beta*log(k) + const
    A    = np.column_stack([np.log(bk), np.ones(len(bk))])
    sol, _, _, _ = np.linalg.lstsq(A, np.log(brms), rcond=None)
    neg_beta, log_A0 = sol
    beta = -neg_beta

    # half-decades to halve: solve k^beta = 2  => log10(k_halve) = log10(2)/beta
    half_decades = math.log10(2) / beta if beta > 0 else float('inf')

    print(f"  m={m}: β = {beta:+.3f}   RMS ~ k^(-{beta:.3f})   "
          f"halves every {half_decades:.2f} decades")
    conv_results.append((m, beta, half_decades))

    col = COLORS.get(m, '#333333')
    ax1.scatter(bk, brms, color=col, s=30, zorder=3, label=f"m={m}")
    # Fitted power law
    kk = np.logspace(np.log10(bk[0]), np.log10(bk[-1]), 200)
    ax1.plot(kk, math.exp(log_A0) * kk**(-beta),
             color=col, linewidth=1.4, linestyle='--')

# ── All-primes theoretical comparison ────────────────────────────────────────
# For p(n) ~ n*ln n, next term is n*ln(ln n), so relative error ~ ln(ln n)/ln n
# In log space: |ln p(n) - ln(n ln n)| ≈ ln(ln n)/ln n  (for large n)
# As a function of n, this decreases like 1/ln(n) -- sub-logarithmic in n.

n_vals = np.logspace(1, 9, 300)
all_prime_err = np.log(np.log(n_vals)) / np.log(n_vals)  # dimensionless relative error

ax1.plot(n_vals, all_prime_err, color='black', linewidth=2.0, linestyle='-',
         label='All primes: ln(ln n)/ln n (theory)')

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('k  (rank within m-class, or n for all-primes)', fontsize=11)
ax1.set_ylabel('RMS log-residual  |ln(actual) − ln(fitted)|', fontsize=11)
ax1.set_title('Convergence to PNT asymptote: m-class primes vs all primes', fontsize=12)
ax1.legend(loc='upper right', fontsize=9)
ax1.grid(True, which='both', alpha=0.3)

# ── Summary panel ────────────────────────────────────────────────────────────
# Show convergence exponents: m-class vs all primes
ax2.axhline(0, color='black', linewidth=0.6, linestyle=':')

ms   = [r[0] for r in conv_results]
betas= [r[1] for r in conv_results]

ax2.bar([str(m) for m in ms], betas,
        color=[COLORS.get(m,'#999') for m in ms],
        edgecolor='black', linewidth=0.7, label='m-class β (power law decay)')

# All-primes effective exponent: ln(ln n)/ln n ~ const/ln n.
# At k=10^6 (mid-range), "effective β" from finite difference:
# d log(1/ln n)/d log n = -1/ln n ≈ -1/14 ≈ -0.071 at n=10^6
# But it's not a power law — the "exponent" changes with n.
# Show the effective local slope at a few n values for context.
ns_ref = [1e4, 1e6, 1e8]
eff_beta_allp = [1.0/math.log(n) for n in ns_ref]  # d(-ln(ln n/ln n))/d(ln n) ≈ 1/ln n
for n_ref, eb in zip(ns_ref, eff_beta_allp):
    ax2.axhline(eb, color='black', linewidth=1.2, linestyle='--', alpha=0.7)
    ax2.text(len(ms) - 0.4, eb + 0.005,
             f"All primes eff. β ≈ {eb:.3f} @ n=10^{int(round(math.log10(n_ref)))}",
             fontsize=8, color='black', ha='right')

ax2.set_xlabel('m-class', fontsize=11)
ax2.set_ylabel('β  (RMS residual ~ k^(−β))', fontsize=11)
ax2.set_title('Convergence exponent β per m-class  vs  all-primes effective rate', fontsize=11)
ax2.legend(fontsize=9)
ax2.grid(True, axis='y', alpha=0.3)

plt.tight_layout()
plt.show()
