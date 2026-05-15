import sys, math
import numpy as np

path = r'.\out\Prime Numbers\10^9\plot_data.tsv'
print(f'Loading {path} ...')
k_lists = {}
p_lists = {}
with open(path, 'r') as fh:
    fh.readline()
    for line in fh:
        parts = line.split('\t')
        m = int(parts[0])
        approx_rank = int(parts[1])
        prime_value = int(parts[4])
        k_lists.setdefault(m, []).append(approx_rank)
        p_lists.setdefault(m, []).append(prime_value)

print(f'Loaded m-classes: {sorted(k_lists.keys())}')

MIN_K = 20
N_BINS = 10

for m in sorted(k_lists.keys()):
    k = np.array(k_lists[m], dtype=np.float64)
    p = np.array(p_lists[m], dtype=np.float64)
    order = np.argsort(k)
    k, p = k[order], p[order]
    mask = k >= MIN_K
    k, p = k[mask], p[mask]
    if len(k) < 50:
        continue
    lk = np.log(k)
    lp = np.log(p)
    llk = np.log(lk)
    y = lp - llk
    A = np.column_stack([lk, np.ones_like(lk)])
    sol, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    alpha, lnC = sol
    resid = lp - (alpha*lk + llk + lnC)
    overall_rms = math.sqrt(np.mean(resid**2))

    # Bin by k
    k_min, k_max = k.min(), k.max()
    edges = np.logspace(np.log10(max(k_min, 1)), np.log10(k_max), N_BINS + 1)
    bin_k, bin_rms = [], []
    for i in range(N_BINS):
        mask2 = (k >= edges[i]) & (k < edges[i+1])
        if mask2.sum() < 5:
            continue
        bin_k.append(np.sqrt(edges[i] * edges[i+1]))
        bin_rms.append(math.sqrt(np.mean(resid[mask2]**2)))

    if len(bin_k) < 4:
        print(f'm={m}: not enough bins ({len(bin_k)})')
        continue
    bk  = np.array(bin_k)
    brms= np.array(bin_rms)
    Afit = np.column_stack([np.log(bk), np.ones(len(bk))])
    sol2, _, _, _ = np.linalg.lstsq(Afit, np.log(brms), rcond=None)
    neg_beta, log_A0 = sol2
    beta = -neg_beta
    half_decades = math.log10(2) / beta if beta > 0 else float('inf')

    print(f'm={m:2d}: alpha={alpha:.4f}  overall_RMS={overall_rms:.4f}  '
          f'beta={beta:+.3f}  (RMS~k^(-{beta:.3f}), halves every {half_decades:.2f} decades)')
    for bki, brmsi in zip(bk, brms):
        print(f'         k~{bki:>12,.0f}  RMS={brmsi:.4f}')

print()
print('All-primes (theory): effective beta = 1/ln(n)')
for n in [1e4, 1e5, 1e6, 1e7, 1e8, 1e9]:
    eff = 1.0 / math.log(n)
    lnlnn_lnn = math.log(math.log(n)) / math.log(n)
    print(f'  n=10^{math.log10(n):.0f}: eff_beta={eff:.4f}, ln(ln n)/ln n = {lnlnn_lnn:.4f}')
