# -*- coding: utf-8 -*-
import math

# m=3 class_quantiles: (k, element_value)
data = [
    (1,         19),
    (2,         23),
    (5,         67),
    (10,        131),
    (100,       1283),
    (1000,      13151),
    (10000,     143833),
    (100000,    1642441),
    (1000000,   20190613),
    (10000000,  270924649),
]

# Global power-law fit (least squares on log-log, skip k=1)
pts = [(math.log(k), math.log(p)) for k,p in data if k > 1]
n = len(pts)
sx  = sum(x for x,y in pts)
sy  = sum(y for x,y in pts)
sxx = sum(x*x for x,y in pts)
sxy = sum(x*y for x,y in pts)
slope = (n*sxy - sx*sy) / (n*sxx - sx*sx)
intercept = (sy - slope*sx) / n

print(f"Global log-log fit  ln p(k) = slope*ln(k) + intercept")
print(f"  slope={slope:.4f}  intercept={intercept:.4f}")
print()
print(f"  {'k':>10}  {'p(k)':>14}  {'pred_p':>14}  {'resid(ln)':>10}  {'ratio':>6}  local_slope")
print("  " + "-"*73)
prev_lnk = None
prev_lnp = None
for k, p in data:
    ln_k = math.log(k) if k > 1 else None
    if ln_k is not None:
        pred     = math.exp(slope * ln_k + intercept)
        resid    = math.log(p) - (slope * ln_k + intercept)
        ratio    = math.exp(resid)
        pred_s   = f"{pred:>14,.0f}"
        resid_s  = f"{resid:>+10.3f}"
        ratio_s  = f"{ratio:>5.2f}x"
    else:
        pred_s  = f"{'n/a':>14}"
        resid_s = f"{'n/a':>10}"
        ratio_s = f"{'n/a':>6}"

    if prev_lnk is not None and ln_k is not None:
        local_s = (math.log(p) - prev_lnp) / (ln_k - prev_lnk)
        local_str = f"{local_s:.3f}"
    else:
        local_str = "n/a"

    print(f"  {k:>10,}  {p:>14,}  {pred_s}  {resid_s}  {ratio_s}  {local_str}")

    if ln_k is not None:
        prev_lnk = ln_k
        prev_lnp = math.log(p)

print()
print("pi(k) slope (from index fit): 0.9475")
print(f"p(k)  slope (this fit):        {slope:.4f}")
print()
print("Local slope trend (shows curvature from PNT ln factor):")
rows = list(zip(data, data[1:]))
for (k0,p0),(k1,p1) in rows:
    if k0 >= 10:
        s = (math.log(p1)-math.log(p0))/(math.log(k1)-math.log(k0))
        print(f"  k={k0:>8,} to {k1:>10,}: {s:.3f}")
