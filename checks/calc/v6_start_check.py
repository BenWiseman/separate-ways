#!/usr/bin/env python3
"""The inflation companion of Separate Ways v6, Appendix A: where the mode solver may start each mode.

Each mode starts in the instantaneous ground state of m(t) sigma_z + x sigma_x deep in inflation.
That state misses the evolved one by an amplitude of about m'/(2 x^2), and the miss interferes
with the production amplitude later on, so a start where m is small but m' is not can bias I.
v6_reheating_production.py therefore starts where t <= -5, m <= 0.02 x and m' <= 1e-6 x^2. This
script compares, for a sudden end of reheating at T_rh = 3e12 GeV, the production integral from
  a start with no condition on m'   (t <= -5 and m <= 0.02 x only),
  the solver's start               (t <= -5, m <= 0.02 x, m' <= 1e-6 x^2),
  a deeper start                   (m <= 1e-6 x and m' <= 1e-9 x^2),
and reports the mass shifts, M ~ I^(-2/5). Plant: the start with no condition on m' must differ from the deep one by
more than 0.1 per cent at H_inf = 1e11 GeV, so the comparison can fail.
"""
import math, sys
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import v6_reheating_production as R
fails = 0
def check(ok, what):
    global fails
    print("   %-72s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)

MP, GSTAR = 2.435e18, 106.75
H_of = lambda T: math.sqrt(math.pi**2*GSTAR/90)*T*T/MP
solver = R.n_exact
def start_at(cond):
    def n(x, prof, t_start=None, t_end=None):
        t0 = min(-5.0, prof.t_e - 1.0)
        while cond(prof, t0, x) and t0 > -1e8: t0 *= 2.0
        return solver(x, prof, t_start=t0, t_end=t_end)
    return n
loose = start_at(lambda p, t, x: p.m(t) > 0.02*x)
deep = start_at(lambda p, t, x: p.m(t) > 1e-6*x or p.mdot(t) > 1e-9*x*x)

def I_with(n, prof):
    R.n_exact = n
    try: return R.production_integral(prof)[0]
    finally: R.n_exact = solver

if __name__ == "__main__":
    shifts = {}
    for H_inf, M in ((1e11, 1.08e9), (1e12, 1.076e9)):
        T = 3e12
        t_rh = math.sqrt(M/H_of(T)); s = (8*M*t_rh/H_inf)**(1/3.)
        prof = R.Profile(t_rh=t_rh, t_e=min(s - t_rh, t_rh))
        I_loose, I_new, I_deep = I_with(loose, prof), I_with(solver, prof), I_with(deep, prof)
        d_loose = 100*((I_deep/I_loose)**0.4 - 1); d_new = abs(I_new/I_deep - 1)
        shifts[H_inf] = (d_loose, d_new)
        print("   H_inf %.0e GeV, T_rh 3e12 GeV: I with no condition on m' %.6e, solver %.6e, deep %.6e" % (H_inf, I_loose, I_new, I_deep))
        print("      with no condition on m' the mass moves by %.2f per cent; the solver's start changes I by %.1e" % (d_loose, d_new))
    check(max(v[1] for v in shifts.values()) < 1e-5, "the solver's start agrees with the deep start to 1e-5 in I")
    check(abs(shifts[1e11][0]) > 0.1, "plant: with no condition on m' it differs by over 0.1 per cent, so the test can fail")
    sys.exit(1 if fails else 0)
