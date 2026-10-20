#!/usr/bin/env python3
"""The inflation companion of Separate Ways v6, section 3.2: the mass when the expansion is radiation-like from the end of
inflation (instant reheating, or a quartic inflaton minimum), expansion alone, for several H_inf.

De Sitter joins radiation at H = H_inf; in the variables of Separate Ways v6, section 3, the radiation era has m(t) = t
and H/M = 1/t^2, so the junction is at t_e = sqrt(M/H_inf). The mass solves M^(5/2) I(M) = M_*^(5/2) I_*.
Check: as H_inf grows the result approaches the sudden limit, 257.6 PeV (v6_state_at_bang.py).
"""
import math, sys
from scipy.optimize import brentq
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from v6_reheating_production import Profile, production_integral
fails = 0
def check(ok, what):
    global fails
    print("   %-66s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)
M_STAR, I_STAR = 4.916e8, 0.0127596673634

def mass(H_inf):
    f = lambda lnM: 2.5*lnM + math.log(production_integral(Profile(t_rh=math.sqrt(math.exp(lnM)/H_inf)))[0]) \
        - math.log(M_STAR**2.5*I_STAR)
    return math.exp(brentq(f, math.log(1e8), math.log(1e9), xtol=1e-5))

if __name__ == "__main__":
    res = {}
    for H in (1e11, 1e12, 1e13, 1e15):
        res[H] = mass(H)/1e6
        print("   H_inf = %.0e GeV: radiation-like from the end of inflation, M = %.1f PeV" % (H, res[H]))
    check(abs(res[1e15] - 257.6) < 1.0, "large H_inf approaches the sudden limit, 257.6 PeV")
    check(res[1e11] > res[1e12] > res[1e13], "a lower H_inf gives a heavier N1")
    sys.exit(1 if fails else 0)
