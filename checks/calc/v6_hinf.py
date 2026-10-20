#!/usr/bin/env python3
"""The inflation companion of Separate Ways v6, section 3.2: how much the N1 mass depends on the Hubble rate of inflation.

The particles are made near H = M, after inflation; before that the modes are relativistic and
conformal, so inflation only hands them the massless vacuum. That holds when inflation ends well
before H = M. This solves the mass (central graviton-exchange rate, as in v6_mass_with_bath.py) at
H_inf = 1e10, 1e11, 1e12 and 1e13 GeV for several reheating temperatures and reports the spread.
"""
import sys
from multiprocessing import Pool
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import v6_mass_with_bath as mb
fails = 0
def check(ok, what):
    global fails
    print("   %-66s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)

def job(a):
    T, H = a
    M, fb = mb.mass(T, H, 1.2e-2, 1.0)
    return T, H, M/1e6, fb

if __name__ == "__main__":
    Ts, Hs = (3e12, 1.2e13, 1.87e13, 5e13, 1e14), (1e10, 1e11, 1e12, 1e13)
    jobs = [(T, H) for T in Ts for H in Hs if H >= mb.CH*T*T/mb.MP]
    with Pool(min(20, len(jobs))) as p:
        res = p.map(job, jobs)
    M = {(T, H): m for T, H, m, fb in res}
    print("   T_rh (GeV)   M at H_inf = 1e10   1e11      1e12      1e13 GeV (PeV)")
    for T in Ts:
        print("   %.2e     " % T + "  ".join("%8.1f" % M[(T, H)] if (T, H) in M else "     -  " for H in Hs))
    hi = max(abs(M[(T, H)]/M[(T, 1e12)] - 1) for T in Ts for H in (1e11, 1e13) if (T, H) in M)
    lo = max(M[(T, 1e10)]/M[(T, 1e12)] - 1 for T in Ts if (T, 1e10) in M)
    print("   largest shift from 1e12 GeV: %.2f%% for H_inf = 1e11 or 1e13; %.1f%% for H_inf = 1e10" % (100*hi, 100*lo))
    # with the solver's deep start (v6_start_check.py) the mass moves by at most 0.5 per cent between
    # H_inf = 1e11 and 1e13 GeV
    check(hi < 0.006, "between 1e11 and 1e13 GeV the mass moves by under 0.6 per cent")
    check(lo > 0.01, "plant: at 1e10 GeV it moves by more than 1 per cent, so the test can fail")
    sys.exit(1 if fails else 0)
