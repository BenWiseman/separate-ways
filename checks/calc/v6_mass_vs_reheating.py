#!/usr/bin/env python3
"""The inflation companion of Separate Ways v6, section 3.2: the N1 mass as a function of the reheating temperature, from
gravitational (Bogoliubov) production alone, for a CPT-symmetric bang followed by inflation.

The state entering the radiation era is the massless vacuum (the least-occupied member of the CPT
family once inflation intervenes; see v6_reheating_production.py). The mass solves
M^(5/2) I(M) = M_c^(5/2) I_min, with M_c = 491.6 PeV the ceiling of Separate Ways v6, section 4 (taken here
as the normalisation) and I(M) the production integral, which depends on M
through t_rh = sqrt(M/H_rh) = T_M/T_rh and through the length of the reheating phase.
"""
import math, sys
import numpy as np
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from v6_reheating_production import Profile, production_integral

MP = 2.435e18            # reduced Planck mass, GeV
GSTAR = 106.75
CH = math.sqrt(math.pi**2*GSTAR/90.0)          # H = CH T^2 / MP in radiation
M_C = 4.916e8            # GeV, the ceiling (the normalisation)
I_MIN = 0.0127596673634

def H_of_T(T): return CH*T*T/MP
def T_of_H(H): return math.sqrt(H*MP/CH)

def solve_mass(T_rh, H_inf, iters=6, verbose=False):
    M = M_C
    for _ in range(iters):
        H_rh = H_of_T(T_rh)
        if H_inf < H_rh*(1 - 1e-9):
            return None
        t_rh = math.sqrt(M/H_rh)
        t_e = (8.0*t_rh*M/H_inf)**(1.0/3.0) - t_rh
        t_e = min(t_e, t_rh)
        I, parts = production_integral(Profile(t_rh=t_rh, t_e=t_e))
        M_new = M_C*(I_MIN/I)**0.4
        if abs(M_new/M - 1) < 1e-5: M = M_new; break
        M = M_new
    return M, I, t_rh, t_e

def _job(args):
    H_inf, T_rh = args
    return (H_inf, T_rh, solve_mass(T_rh, H_inf))

if __name__ == "__main__":
    from multiprocessing import Pool
    T_M = T_of_H(M_C)
    print("T at which H = M_c in radiation: %.4e GeV (H there %.4e GeV)" % (T_M, H_of_T(T_M)), flush=True)
    jobs = [(H, T) for H in [1e10, 1e11, 1e12, 1e13]
            for T in [8e12, 1.2e13, 1.6e13, 1.87e13, 2.5e13, 3.5e13, 5e13, 7e13, 1e14, 1.5e14, 2e14, 3e14, 5e14, 1e15, 2e15]
            if H >= H_of_T(T)*(1 - 1e-9)]
    with Pool(min(30, len(jobs) + 1)) as pool:
        sudden = pool.apply_async(production_integral, (Profile(t_rh=1e-4),))
        res = pool.map(_job, jobs)
        I0, parts = sudden.get()
    print("\n=== sudden limit (t_rh -> 0): the inflationary ceiling ===")
    print("   I_vac = %.6f  (I_min = %.6f, ratio %.4f)  ->  M = %.1f PeV" % (I0, I_MIN, I0/I_MIN, M_C*(I_MIN/I0)**0.4/1e6))
    print("   parts: body %.5f  tail %.5f  rest %.2e (before the 1/sqrt(pi))" % parts)
    print("\n=== M against T_rh, Bogoliubov production only ===")
    print("   H_inf (GeV)   T_rh (GeV)    t_rh      t_e       I          I/I_min    M (PeV)   E_line (PeV)")
    last = None
    for H_inf, T_rh, r in res:
        if last is not None and H_inf != last: print()
        last = H_inf
        M, I, t_rh, t_e = r
        print("   %.0e      %.2e   %8.4f  %8.4f  %.6f   %7.4f   %7.1f     %6.1f" % (H_inf, T_rh, t_rh, t_e, I, I/I_MIN, M/1e6, M/2e6))
