import math

# Fit p(k) ~ C * k^a * (ln k)^b  (two-term log model, no pure power law)
# The local slope of ln p vs ln k = a + b/ln(k), linearly drifting in 1/ln(k)
# Endpoints of the observed local slope drift:
lnk1, s1 = 4.605, 0.991   # k=10^2
lnk2, s2 = 16.118, 1.128  # k=10^7
b = (s1 - s2) / (1/lnk1 - 1/lnk2)
a = s1 - b / lnk1

data = [(5,67),(10,131),(100,1283),(1000,13151),(10000,143833),
        (100000,1642441),(1000000,20190613),(10000000,270924649)]

# Fit intercept C
vals  = [(math.log(k), math.log(math.log(k)), math.log(p)) for k,p in data]
C     = sum(lnp - a*lnk - b*lnlnk for lnk,lnlnk,lnp in vals) / len(vals)

print(f"Model:  p(k) = C * k^a * (ln k)^b")
print(f"  a = {a:.4f}   b = {b:.4f}   C = {math.exp(C):.4f}")
print()
print(f"  {'k':>10}  {'p(k)':>14}  {'pred_p':>14}  {'ratio':>6}")
print("  " + "-"*50)
for k,p in data:
    lnk   = math.log(k)
    pred  = math.exp(a*lnk + b*math.log(lnk) + C)
    ratio = p / pred
    print(f"  {k:>10,}  {p:>14,}  {pred:>14,.0f}  {ratio:>6.3f}x")

print()
print("For comparison, plain power law (slope=1.0439) residuals:")
slope2, ic2 = 1.0439, 2.3970
for k,p in data:
    pred  = math.exp(slope2*math.log(k) + ic2)
    ratio = p / pred
    print(f"  k={k:>10,}  ratio={ratio:.3f}x")
