#!/usr/bin/env python3
"""The inflation companion of Separate Ways v6, section 3.2: the N1 mass against the reheating temperature, with both sources.

Expansion: Y_exp = K M^(3/2) I(M, T_rh), with K fixed by the half-angle normalisation
           M_*^(5/2) K I_* = rho/s (Separate Ways v6, section 4), and I from v6_reheating_production.py.
Bath:      Y_bath = b alpha (T_rh/M_P)^3 / C (v6_bath.py), b the reheating boost (1 to 2).
The mass solves M (Y_exp + Y_bath) = rho/s; the left side increases with M, so the root is unique.
"""
import math, sys
from multiprocessing import Pool
from scipy.optimize import brentq
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from v6_reheating_production import Profile, production_integral, n_star_avg
from v6_bath import Y_grav, MP, CH, RHO_OVER_S

M_STAR, I_STAR = 4.916e8, 0.0127596673634
K = RHO_OVER_S/(M_STAR**2.5*I_STAR)

CEILING = "--ceiling" in sys.argv     # every CPT-invariant state: the bound of Separate Ways v6, section 4, with two half-crossings

def I_of(M, T_rh, H_inf):
    H_rh = CH*T_rh**2/MP
    t_rh = math.sqrt(M/H_rh)
    t_e = min((8.0*t_rh*M/H_inf)**(1/3) - t_rh, t_rh)
    return production_integral(Profile(t_rh=t_rh, t_e=t_e), transform=n_star_avg if CEILING else None)[0]

def mass(T_rh, H_inf, alpha, boost):
    Yb = boost*Y_grav(T_rh, alpha)
    f = lambda lnM: math.log(math.exp(lnM)*(K*math.exp(lnM)**1.5*I_of(math.exp(lnM), T_rh, H_inf) + Yb)/RHO_OVER_S)
    lnM = brentq(f, math.log(1e3), math.log(5e9), xtol=1e-5)
    M = math.exp(lnM)
    Ye = K*M**1.5*I_of(M, T_rh, H_inf)
    return M, Yb/(Ye + Yb)

def job(a):
    T_rh, H_inf, alpha, boost, tag = a
    try:
        M, fb = mass(T_rh, H_inf, alpha, boost)
    except ValueError:
        M, fb = float('nan'), float('nan')
    return (tag, T_rh, H_inf, M, fb)

if __name__ == "__main__":
    import numpy as np
    Ts = sorted(set([3e12, 5e12, 8e12, 1.2e13, 1.87e13, 2.5e13, 3.5e13, 5e13, 7e13, 1e14, 1.4e14, 2e14, 3e14]
                    + [float("%.3g" % t) for t in np.geomspace(2e12, 4e14, 25)]))
    variants = [("expansion only", 1e12, 0.0, 1.0),
                ("central: alpha 1.2e-2, no boost", 1e12, 1.2e-2, 1.0),
                ("high: alpha 1.7e-2, boost 2", 1e12, 1.7e-2, 2.0),
                ("central, H_inf 1e13", 1e13, 1.2e-2, 1.0)]
    if CEILING:
        variants = [("ceiling, expansion only", 1e12, 0.0, 1.0), ("ceiling, central", 1e12, 1.2e-2, 1.0)]
    jobs = [(T, H, a, b, tag) for tag, H, a, b in variants for T in Ts if H >= CH*T*T/MP]
    with Pool(min(30, len(jobs))) as pool:
        res = pool.map(job, jobs)
    out = __file__.rsplit("/", 1)[0] + ("/v6_mass_ceiling.tsv" if CEILING else "/v6_mass_with_bath.tsv")
    with open(out, "w") as f:
        f.write("variant\tT_rh_GeV\tH_inf_GeV\tM_GeV\tbath_fraction\n")
        for r in res: f.write("%s\t%.4e\t%.1e\t%.6e\t%.6f\n" % r)
    print("wrote %s" % out)
    for tag, H, a, b in variants:
        print("\n=== %s (H_inf = %.0e GeV) ===" % (tag, H))
        print("   T_rh (GeV)    M (PeV)   line (PeV)   from the bath")
        for r in res:
            if r[0] == tag:
                print("   %.2e     %7.1f     %6.1f      %5.1f%%" % (r[1], r[3]/1e6, r[3]/2e6, 100*r[4]))
