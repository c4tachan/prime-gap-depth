"""
rh_consistency.py

Tests whether our m-class data is consistent with the Riemann Hypothesis prediction:

    alpha(k)  -->  1  at rate  k^{-1/2}

i.e.  |alpha_local(k) - 1|  ~  A * k^{-0.5}

Method:
  1. Load plot_data.tsv
  2. For each m-class, compute the local log-log slope d(ln p)/d(ln k) in sliding
     log-spaced bins (using linear regression within each bin).
  3. Plot |alpha_local(k) - 1| vs k on log-log axes.
  4. Fit a power law decay, report the exponent gamma.
  5. Compare gamma to the RH-predicted 0.5.

Note on what this test actually tells us:
  - If gamma ~ 0.5 : data CONSISTENT with RH (not a proof)
  - If gamma >> 0.5 : converging faster than RH requires (also consistent)
  - If gamma << 0.5 : slower convergence, would suggest RH-like regularity is absent
  - The test is weak because our k-range spans only ~5-8 decades, and the
    approach rate itself changes slowly (logarithmically), so we can only
    distinguish gamma=0.5 from gamma=0.0 or gamma=1.0, not finer differences.

Usage:
    python scripts/rh_consistency.py ".\out\Prime Numbers\10^9\plot_data.tsv"
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── Config ────────────────────────────────────────────────────────────────────

path = sys.argv[1] if len(sys.argv) > 1 else r".\out\Prime Numbers\10^9\plot_data.tsv"
N_BINS = 12          # log-spaced k bins per m-class
MIN_PTS_PER_BIN = 6  # minimum samples required to compute a local slope
MIN_K = 50           # discard very small k (non-asymptotic regime)

COLORS = {2:'#e41a1c', 3:'#ff7f00', 4:'#4daf4a', 5:'#377eb8',
          6:'#984ea3', 7:'#a65628'}

# ── Load data ─────────────────────────────────────────────────────────────────

print(f"Loading {path} ...")
k_lists: dict[int, list] = {}
p_lists: dict[int, list] = {}
with open(path, "r") as fh:
    fh.readline()
    for line in fh:
        parts = line.split("\t")
        m            = int(parts[0])
        approx_rank  = int(parts[1])
        prime_value  = int(parts[4])
        k_lists.setdefault(m, []).append(approx_rank)
        p_lists.setdefault(m, []).append(prime_value)

print(f"  m-classes: {sorted(k_lists.keys())}")

# ── Per-class: compute local slope alpha(k) in bins ───────────────────────────

fig = plt.figure(figsize=(14, 10))
gs  = gridspec.GridSpec(2, 1, hspace=0.40)
ax_alpha  = fig.add_subplot(gs[0])   # |alpha(k) - 1| vs k
ax_gamma  = fig.add_subplot(gs[1])   # fitted gamma per m-class

print()
print(f"{'m':>2}  {'global α':>9}  {'fit γ':>7}  {'RH pred':>8}  "
      f"{'consistent?':>12}  {'bins used':>9}")
print("─" * 60)

gamma_results = []

for m in sorted(k_lists.keys()):
    k_all = np.array(k_lists[m], dtype=np.float64)
    p_all = np.array(p_lists[m], dtype=np.float64)
    order = np.argsort(k_all)
    k_all, p_all = k_all[order], p_all[order]

    mask = k_all >= MIN_K
    k_all, p_all = k_all[mask], p_all[mask]
    if len(k_all) < 30:
        continue

    lk_all = np.log(k_all)
    lp_all = np.log(p_all)

    # Global slope (for reference)
    A_g = np.column_stack([lk_all, np.ones(len(lk_all))])
    sol_g, _, _, _ = np.linalg.lstsq(A_g, lp_all, rcond=None)
    alpha_global = sol_g[0]

    # Log-spaced bins
    edges = np.logspace(np.log10(k_all.min()), np.log10(k_all.max()), N_BINS + 1)
    bin_k_mid  = []
    bin_alpha  = []

    for i in range(N_BINS):
        mask_b = (k_all >= edges[i]) & (k_all < edges[i+1])
        if mask_b.sum() < MIN_PTS_PER_BIN:
            continue
        lk_b = lk_all[mask_b]
        lp_b = lp_all[mask_b]
        A_b  = np.column_stack([lk_b, np.ones(len(lk_b))])
        sol_b, _, _, _ = np.linalg.lstsq(A_b, lp_b, rcond=None)
        alpha_b = sol_b[0]
        bin_k_mid.append(math.sqrt(edges[i] * edges[i+1]))
        bin_alpha.append(alpha_b)

    if len(bin_k_mid) < 4:
        continue

    bk    = np.array(bin_k_mid)
    ba    = np.array(bin_alpha)
    dev   = np.abs(ba - 1.0)

    # Only use bins where deviation is positive and measurable
    valid = dev > 1e-6
    if valid.sum() < 3:
        continue
    bk_v, dev_v = bk[valid], dev[valid]

    # Fit: log|alpha - 1| = -gamma * log(k) + const
    A_fit = np.column_stack([np.log(bk_v), np.ones(len(bk_v))])
    sol_f, _, _, _ = np.linalg.lstsq(A_fit, np.log(dev_v), rcond=None)
    neg_gamma, log_A0 = sol_f
    gamma = -neg_gamma

    consistent = "YES" if 0.3 <= gamma <= 0.7 else ("faster" if gamma > 0.7 else "SLOWER")

    print(f"{m:>2}  {alpha_global:>9.4f}  {gamma:>7.3f}  {'0.500':>8}  "
          f"{consistent:>12}  {len(bk_v):>9}")
    gamma_results.append((m, gamma, alpha_global, bk_v, dev_v, log_A0))

    col = COLORS.get(m, '#888888')

    # Plot |alpha - 1| data points
    ax_alpha.scatter(bk_v, dev_v, color=col, s=40, zorder=3,
                     label=f"m={m}  (γ={gamma:.3f})")

    # Fitted power law
    kk = np.logspace(np.log10(bk_v[0]), np.log10(bk_v[-1]), 200)
    ax_alpha.plot(kk, math.exp(log_A0) * kk**(-gamma),
                  color=col, linewidth=1.4, linestyle='--')

# RH reference line slope -0.5
k_ref = np.logspace(1, 9, 300)
# Anchor the RH line to be visually useful: pass through a middle point
# Use m=3 first bin as anchor if available
if gamma_results:
    m0, g0, ag0, bk0, dev0, la0 = gamma_results[1] if len(gamma_results) > 1 else gamma_results[0]
    anchor_k = math.sqrt(bk0[0] * bk0[-1])
    anchor_v = math.exp(la0) * anchor_k**(-g0)    # actual fit at midpoint
    rh_A     = anchor_v * anchor_k**0.5            # scale RH line to same midpoint
    rh_y     = rh_A * k_ref**(-0.5)
    ax_alpha.plot(k_ref, rh_y, color='black', linewidth=2.0, linestyle='-',
                  label='RH prediction: γ = 0.5', zorder=2)

ax_alpha.set_xscale('log')
ax_alpha.set_yscale('log')
ax_alpha.set_xlabel('k  (rank within m-class)', fontsize=11)
ax_alpha.set_ylabel('|α_local(k) − 1|', fontsize=11)
ax_alpha.set_title(
    'Approach of local log-log slope to α=1\n'
    'RH predicts: |α(k)−1| ~ k^{−1/2}   (dashed = fitted, solid black = RH)',
    fontsize=11)
ax_alpha.legend(fontsize=9, loc='upper right')
ax_alpha.grid(True, which='both', alpha=0.3)

# ── Bottom panel: measured gamma vs RH prediction ────────────────────────────

ms     = [r[0] for r in gamma_results]
gammas = [r[1] for r in gamma_results]
cols   = [COLORS.get(m, '#888') for m in ms]

bars = ax_gamma.bar([f"m={m}" for m in ms], gammas,
                    color=cols, edgecolor='black', linewidth=0.8)
ax_gamma.axhline(0.5, color='black', linewidth=2.0, linestyle='-',
                 label='RH prediction γ = 0.5')
ax_gamma.axhspan(0.3, 0.7, color='black', alpha=0.08,
                 label='Rough consistency band ±0.2')

# Annotate bars
for bar, g in zip(bars, gammas):
    ax_gamma.text(bar.get_x() + bar.get_width()/2, g + 0.01,
                  f'{g:.3f}', ha='center', va='bottom', fontsize=9)

ax_gamma.set_ylim(0, max(gammas) * 1.3)
ax_gamma.set_ylabel('Fitted decay exponent γ', fontsize=11)
ax_gamma.set_title(
    'Measured γ per m-class vs RH-predicted 0.5\n'
    'γ > 0.5 means faster approach (still RH-consistent); γ < 0.3 would be anomalous',
    fontsize=11)
ax_gamma.legend(fontsize=9)
ax_gamma.grid(True, axis='y', alpha=0.3)

print()
print("Interpretation:")
print("  RH predicts |alpha(k) - 1| ~ k^{-0.5} for all m-classes.")
print("  gamma > 0.5  =>  converging FASTER than RH requires (consistent)")
print("  gamma ~ 0.5  =>  consistent with RH")
print("  gamma < 0.3  =>  converging slower than RH implies (anomalous)")
print()
print("  Note: this is a CONSISTENCY CHECK, not a proof.")
print("  Our k-range (~5-8 decades) is too short to distinguish gamma=0.4 from")
print("  gamma=0.6 reliably. The test can only rule out grossly wrong behaviour.")

plt.tight_layout()
plt.show()
