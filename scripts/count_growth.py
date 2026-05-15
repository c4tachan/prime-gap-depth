"""
count_growth.py

Reads class_quantiles.tsv from each decade directory and extracts the class
size (count of m-class primes below N) for each m-class.  Then fits three
candidate growth models to count_m(N) and shows which fits best:

  1. ln(ln(N))          -- gap-record / extremal classes (m=1, m=6)
  2. N / ln(N)          -- constant-fraction-of-primes (pi(N)-like)
  3. N^alpha / ln(N)    -- sub-linear power (intermediate classes)

Usage:
    python scripts/count_growth.py
"""

import math
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

BASE = r"f:\src\prime-gap-depth\out\Prime Numbers"
DECADES = [1, 2, 3, 4, 5, 6, 7, 8, 9]

COLORS = {0:'#888888', 1:'#e41a1c', 2:'#ff7f00', 3:'#4daf4a',
          4:'#377eb8', 5:'#984ea3', 6:'#a65628', 7:'#f781bf'}

# ── Collect class sizes from each decade ──────────────────────────────────────

# counts[m][decade] = class_size
counts: dict[int, dict[int, int]] = {}

for dec in DECADES:
    path = os.path.join(BASE, f"10^{dec}", "class_quantiles.tsv")
    if not os.path.exists(path):
        continue
    with open(path) as fh:
        fh.readline()  # header
        seen = set()
        for line in fh:
            parts = line.split("\t")
            m    = int(parts[0])
            size = int(parts[1])
            if m not in seen:
                counts.setdefault(m, {})[dec] = size
                seen.add(m)

# Print raw table
all_m = sorted(counts.keys())
print(f"{'N':>12}  " + "  ".join(f"m={m:>2}" for m in all_m))
print("─" * (14 + 7 * len(all_m)))
for dec in DECADES:
    N = 10 ** dec
    row = f"10^{dec:>2} ({N:>12,})  "
    for m in all_m:
        val = counts.get(m, {}).get(dec, None)
        row += f"{val:>7,}" if val is not None else f"{'—':>7}"
    print(row)

# ── Fit models for each m-class ───────────────────────────────────────────────

print()
print("Growth model fits (count_m vs N):")
print(f"{'m':>2}  {'model':>20}  {'alpha':>7}  {'R²':>6}  description")
print("─" * 70)

fig = plt.figure(figsize=(15, 11))
gs  = gridspec.GridSpec(2, 2, hspace=0.38, wspace=0.32)
ax_raw    = fig.add_subplot(gs[0, :])   # raw counts, full width top
ax_lnln   = fig.add_subplot(gs[1, 0])  # test: count vs ln(ln(N))
ax_pow    = fig.add_subplot(gs[1, 1])  # test: log(count) vs log(N)

fit_results = []

for m in all_m:
    if m not in counts:
        continue
    dec_vals = sorted(counts[m].items())
    if len(dec_vals) < 3:
        continue

    Ns   = np.array([10.0 ** d for d, _ in dec_vals])
    cnts = np.array([float(c)  for _, c in dec_vals])
    col  = COLORS.get(m, '#333333')

    # Raw plot
    ax_raw.plot(Ns, cnts, 'o-', color=col, linewidth=1.5, markersize=5,
                label=f"m={m}")

    # ── Model 1: count ~ a * ln(ln(N)) + b ───────────────────────────────
    x_lnln = np.log(np.log(Ns))
    A1 = np.column_stack([x_lnln, np.ones(len(Ns))])
    sol1, _, _, _ = np.linalg.lstsq(A1, cnts, rcond=None)
    pred1 = A1 @ sol1
    ss_res1 = np.sum((cnts - pred1)**2)
    ss_tot  = np.sum((cnts - cnts.mean())**2)
    r2_1 = 1 - ss_res1/ss_tot if ss_tot > 0 else 0

    # ── Model 2: count ~ a * N/ln(N) ─────────────────────────────────────
    x_pnt = Ns / np.log(Ns)
    A2 = np.column_stack([x_pnt, np.ones(len(Ns))])
    sol2, _, _, _ = np.linalg.lstsq(A2, cnts, rcond=None)
    pred2 = A2 @ sol2
    ss_res2 = np.sum((cnts - pred2)**2)
    r2_2 = 1 - ss_res2/ss_tot if ss_tot > 0 else 0

    # ── Model 3: log(count) ~ alpha*log(N) + c  (pure power law) ─────────
    log_N   = np.log(Ns)
    log_cnt = np.log(cnts)
    A3 = np.column_stack([log_N, np.ones(len(Ns))])
    sol3, _, _, _ = np.linalg.lstsq(A3, log_cnt, rcond=None)
    alpha3, logC3 = sol3
    pred3_log = A3 @ sol3
    ss_res3 = np.sum((log_cnt - pred3_log)**2)
    ss_tot3 = np.sum((log_cnt - log_cnt.mean())**2)
    r2_3 = 1 - ss_res3/ss_tot3 if ss_tot3 > 0 else 0

    best = max([(r2_1,'ln(ln N)'), (r2_2,'N/ln N'), (r2_3,f'N^{alpha3:.3f}')],
               key=lambda x: x[0])
    print(f"{m:>2}  {best[1]:>20}  {alpha3:>7.4f}  {best[0]:>6.4f}  "
          f"lnln R²={r2_1:.4f}  pnt R²={r2_2:.4f}  pow R²={r2_3:.4f}")
    fit_results.append((m, r2_1, r2_2, r2_3, alpha3, col))

    # ln(ln N) test plot
    ax_lnln.plot(np.log(np.log(Ns)), cnts, 'o-', color=col,
                 markersize=5, linewidth=1.2, label=f"m={m}")

    # log-log plot for power law
    ax_pow.plot(np.log10(Ns), np.log10(cnts), 'o-', color=col,
                markersize=5, linewidth=1.2, label=f"m={m}")

# ── Format axes ───────────────────────────────────────────────────────────────

ax_raw.set_xscale('log')
ax_raw.set_yscale('log')
ax_raw.set_xlabel('N', fontsize=11)
ax_raw.set_ylabel('count_m(N)', fontsize=11)
ax_raw.set_title('m-class counts vs N  (log-log)', fontsize=11)
ax_raw.legend(fontsize=9, ncol=4)
ax_raw.grid(True, which='both', alpha=0.3)

ax_lnln.set_xlabel('ln(ln N)', fontsize=11)
ax_lnln.set_ylabel('count_m(N)', fontsize=11)
ax_lnln.set_title('count vs ln(ln N)\nstraight line => ln(ln N) growth', fontsize=10)
ax_lnln.legend(fontsize=8)
ax_lnln.grid(True, alpha=0.3)

ax_pow.set_xlabel('log₁₀(N)', fontsize=11)
ax_pow.set_ylabel('log₁₀(count)', fontsize=11)
ax_pow.set_title('log-log count vs N\nslope = power law exponent α', fontsize=10)
ax_pow.legend(fontsize=8)
ax_pow.grid(True, alpha=0.3)

# Add slope reference lines on log-log panel
for alpha_ref, label in [(1.0, 'α=1 (∝N)'), (0.5, 'α=0.5'), (0.0, 'α=0 (const)')]:
    x0, x1 = 1.0, 9.0
    # anchor at x=5 (N=10^5) through a neutral point
    y0 = 5.0 * alpha_ref
    ax_pow.plot([x0, x1], [y0 + (x0-5)*alpha_ref, y0 + (x1-5)*alpha_ref],
                color='black', linewidth=0.8, linestyle=':', alpha=0.5)
    ax_pow.text(x1 + 0.05, y0 + (x1-5)*alpha_ref, label, fontsize=7, va='center')

plt.suptitle("How do m-class prime counts grow with N?\n"
             "ln(ln N) vs N/ln(N) vs power law", fontsize=12, y=1.01)
plt.tight_layout()
plt.show()
