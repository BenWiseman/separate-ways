#!/usr/bin/env python3
"""Separate Ways v6, section 4.3 and Appendix A: without inflation the ceiling has no floor.

CPT invariance fixes each mode half in each component at the bang and leaves one phase per mode.
Over that phase the late occupation runs between n_min = (1 - sqrt(1 - P))/2 and n_max = 1 - n_min,
P = exp(-pi x^2) (x = p/sqrt(gamma)). A state at the band edge n_max out to x_c and at n_min beyond
is CPT-invariant and Hadamard (it departs from the half-angle state only below x_c), and its
production integral I = pi^(-1/2) Int x^2 n dx grows without limit with x_c, so M = M_* (I_*/I)^(2/5)
falls without limit. The table reproduces v5.6 Table 1 (there in the variable sqrt(pi) x).
"""
import math, sys
from scipy.integrate import quad
fails = 0
def check(ok, what):
    global fails
    print("   %-66s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)

M_STAR = 491.6
def nmin(x):
    u = math.exp(-math.pi*x*x)
    return u/(2*(1 + math.sqrt(1 - u)))          # (1 - sqrt(1-u))/2 without cancellation
def I_state(xc):
    lo = quad(lambda x: x*x*(1 - nmin(x)), 0, xc, limit=200)[0] if xc > 0 else 0.0
    hi = quad(lambda x: x*x*nmin(x), xc, 12, limit=400)[0]
    return (lo + hi)/math.sqrt(math.pi)

if __name__ == "__main__":
    I0 = I_state(0.0)
    check(abs(I0 - 0.0127596673634) < 1e-10, "half-angle state: I_* = %.10f" % I0)
    print("   x_c (this paper)   x_c (v5.6)   I/I_*    M (PeV)   line (PeV)")
    rows = []
    for xs5 in (0.10, 0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 3.00):
        xc = xs5/math.sqrt(math.pi)
        r = I_state(xc)/I0; M = M_STAR*r**-0.4
        rows.append((xs5, r, M))
        print("   %8.3f           %5.2f      %7.3f   %6.1f    %6.1f" % (xc, xs5, r, M, M/2))
    check(abs(rows[5][1] - 8.312) < 2e-3 and abs(rows[5][2] - 210.7) < 0.1, "reproduces v5.6 Table 1 at x_c = 1.50: I/I_* = 8.312, 210.7 PeV")
    check(all(rows[i+1][2] < rows[i][2] for i in range(len(rows)-1)), "the mass falls monotonically as the excursion widens")
    check(rows[-1][2] < 100, "plant: a wide enough excursion takes the mass below 100 PeV (no floor)")
    # a step smoothed over a width w: n = n_min + (1 - 2 n_min)/(1 + exp((x - x_c)/w)); the change in I vanishes with w
    xc = 1/math.sqrt(math.pi)
    def I_smooth(w):
        n = lambda x: nmin(x) + (1 - 2*nmin(x))/(1 + math.exp(min((x - xc)/w, 700)))
        return quad(lambda x: x*x*n(x), 0, 12, points=[xc], limit=400)[0]/math.sqrt(math.pi)
    sharp = M_STAR*(I_state(xc)/I0)**-0.4
    for w in (0.05, 0.01):
        print("   step at x_c = %.3f smoothed over w = %.2f: M = %.2f PeV, against %.2f PeV for the sharp step"
              % (xc, w, M_STAR*(I_smooth(w)/I0)**-0.4, sharp))
    d05, d01 = abs(M_STAR*(I_smooth(0.05)/I0)**-0.4 - sharp), abs(M_STAR*(I_smooth(0.01)/I0)**-0.4 - sharp)
    print("   the change shrinks by %.1f for a width five times smaller (25 for a change that goes as w^2); at w = 0.01 it is"
          " %.2f per cent" % (d05/d01, 100*d01/sharp))
    check(15 < d05/d01 < 35 and d01 < 0.002*sharp, "the change from smoothing vanishes as the square of the width")
    # where physics sets a floor: the dark matter's comoving momentum is p = x sqrt(M a_1) and a_1 = pi (aT)^2
    # sqrt(g*/15)/mu (section 4.2), so p over the plasma temperature is x sqrt(pi M sqrt(g*/15)/mu); over a filled band
    # out to x_c the x^2-weighted mean of x is 3x_c/4, and for x_c >> 1, I -> x_c^3/(3 sqrt(pi))
    MU, GS = 5.965e18, 106.75
    from scipy.optimize import brentq
    M_big = lambda xc: M_STAR*1e6*((xc**3/(3*math.sqrt(math.pi)))/I0)**-0.4          # GeV
    p_over_T = lambda xc: 0.75*xc*math.sqrt(math.pi*M_big(xc)*math.sqrt(GS/15)/MU)
    lx = brentq(lambda l: math.log(p_over_T(math.exp(l))), 5, 80)
    xc1 = math.exp(lx); M1 = M_big(xc1)
    lx5 = brentq(lambda l: math.log(M_big(math.exp(l))/5e-6), 5, 80)
    print("   the mean momentum reaches the plasma temperature at x_c = %.1e, M = %.2f keV; at M = 5 keV it is %.2f of it"
          % (xc1, M1*1e6, p_over_T(math.exp(lx5))))
    check(1e-7 < M1 < 1e-5, "the momenta reach the plasma temperature near a keV, where warm-dark-matter bounds apply")
    print("   after the fall of g* from %.2f to 3.91 the dark matter's momenta fall against the photons' temperature by %.2f"
          % (GS, (GS/3.91)**(1/3)))
    # the energy these states hold at production, against the radiation's: with the abundance fixed, n/s = rho_DM0/(s0 M),
    # and mean momentum T, rho_N/rho_rad = (4/3) n T/s/T... = (4/3) Y
    Y1 = (9.7e-48/2.2215e-38)/M1
    print("   at that mass the dark matter holds %.1e of the radiation's energy at production" % (4/3*Y1))
    sys.exit(1 if fails else 0)
