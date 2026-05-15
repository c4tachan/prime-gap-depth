import math

# The theoretically motivated model:
#   pi(k) ~ k^alpha  (observed alpha = 0.9475 from index fit)
#   p_n   ~ n * ln(n)   (PNT)
#   => p(k) ~ pi(k) * ln(pi(k)) ~ C * k^alpha * ln(k)
#
# ln p(k) = alpha*ln(k) + ln(ln(k)) + ln(C)
#
# Note: ln(pi(k)) ~ alpha*ln(k) for large k, so the coefficient on ln(ln k) is ~1.

alpha = 0.9475  # from index fit

data = [(10,131),(100,1283),(1000,13151),(10000,143833),
        (100000,1642441),(1000000,20190613),(10000000,270924649)]

# Fit ln(C) as average residual
lnC_list = [math.log(p) - alpha*math.log(k) - math.log(math.log(k)) for k,p in data]
lnC = sum(lnC_list) / len(lnC_list)
C   = math.exp(lnC)

print("Theoretically motivated model (alpha pinned from pi(k) fit):")
print(f"  p(k) = C * k^alpha * ln(k)")
print(f"  alpha = {alpha}   C = {C:.4f}   ln(C) = {lnC:.4f}")
print()
print(f"  {'k':>10}  {'p(k)':>14}  {'pred_p':>14}  {'ratio':>6}  resid(ln)")
print("  " + "-"*58)
residuals = []
for k,p in data:
    pred  = C * (k**alpha) * math.log(k)
    resid = math.log(p) - math.log(pred)
    residuals.append(resid)
    print(f"  {k:>10,}  {p:>14,}  {pred:>14,.0f}  {p/pred:>6.3f}x  {resid:+.3f}")

rms = math.sqrt(sum(r*r for r in residuals)/len(residuals))
print()
print(f"  RMS residual (ln-scale): {rms:.4f}  => typical ratio: {math.exp(rms):.3f}x")
print()

# Extrapolate to full class (k=428,819,600) and compare with pure power law
k_full = 428_819_600
pred_natural = C * (k_full**alpha) * math.log(k_full)
pred_powerlaw = math.exp(1.0439 * math.log(k_full) + 2.3970)
print(f"Extrapolation to k = {k_full:,}  (last m=3 prime at N=10^9):")
print(f"  Natural model:  p ~ {pred_natural:>18,.0f}")
print(f"  Pure power law: p ~ {pred_powerlaw:>18,.0f}")
print()
print(f"  The two agree within {abs(pred_natural/pred_powerlaw - 1)*100:.1f}% at this extrapolation distance.")
print(f"  But the natural model has a theoretical justification; the power law does not.")
print()

# Show the uncertainty band: +/- 1 RMS residual
low  = pred_natural * math.exp(-rms)
high = pred_natural * math.exp(+rms)
print(f"  +/- 1 RMS band:  [{low:>18,.0f},  {high:>18,.0f}]")
print(f"  That is a factor of {high/low:.2f}x width.")
