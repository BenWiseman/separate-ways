#!/usr/bin/env python3
"""The inflation companion of Separate Ways v6, section 3.2: the plasma's yield for a quartic inflaton minimum when the decay
rate follows the inflaton's falling mass.

In a quartic minimum the oscillating inflaton redshifts like radiation, so the expansion, and with it
the expansion's share of N1, does not depend on how the inflaton decays. Its effective mass falls as
1/a, though, so a constant decay rate is not physical: a decay to fermions through a Yukawa coupling
has Gamma ~ m_eff ~ 1/a, and one to bosons through a trilinear coupling Gamma ~ 1/m_eff ~ a. This
integrates rho_phi' = -4 H rho_phi - Gamma rho_phi and rho_r' = -4 H rho_r + Gamma rho_phi with the
graviton-exchange yield along the history, for Gamma = g0 a^p with p = -1, 0, +1, and T_rh defined as
in v6_mass_smooth.py: the radiation-era temperature at which H equals Gamma where the two cross.
At the reheating temperature above which the plasma makes most of the dark matter, the mass is
2^(-2/5) of the expansion-only mass, 261.8 PeV (v6_mass_smooth.py --w13), and the plasma's yield
equals the expansion's; the script finds that temperature for each decay law, in metric gravity and
with the minimal Einstein-Cartan term.
Plant: p = 0 must reproduce the constant-rate crossover, 8.97e13 GeV (v6_figures.py), to 2 per cent.
"""
import math, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
fails = 0
def check(ok, what):
    global fails
    print("   %-72s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)

MP, GSTAR = 2.435e18, 106.75
CH = math.sqrt(math.pi**2*GSTAR/90)
RHO_OVER_S = 9.7e-48/2.2215e-38
ALPHA, ALPHA_EC = 1.2027e-2, 1.2027e-2 + 7.4185e-2
H_INF, M_EXP = 1e12, 2.618e8                 # GeV; expansion-only mass for a quartic minimum
q = 30.0/(math.pi**2*GSTAR)

def history(g0, p):
    """Units a_e = 1, time 1/H_inf, rho in M_P^2 H_inf^2. Returns (T_rh, plasma yield per unit alpha)."""
    def f(t, y):
        a, rp, rr, na3 = y
        H = math.sqrt(max(rp + rr, 0.0)/3.0)
        G = g0*a**p
        return [a*H, -4*H*rp - G*rp, -4*H*rr + G*rp, (q*rr)**2*a**3]
    done = lambda t, y: y[1]/(y[2] + 1e-300) - 1e-9
    done.terminal = True
    sol = solve_ivp(f, (0.0, 1e12), [1.0, 3.0, 0.0, 0.0], method="LSODA", rtol=1e-10, atol=1e-40,
                    events=done, dense_output=True)
    a, rp, rr, na3 = sol.y[:, -1]
    # where Gamma and H cross, in units of H_inf, and the radiation-era temperature with that H
    ts = np.geomspace(1e-6, sol.t[-1], 4000)
    ys = sol.sol(ts)
    H = np.sqrt((ys[1] + ys[2])/3.0); G = g0*ys[0]**p
    k = int(np.argmax(np.sign(G - H) != np.sign(G[0] - H[0])))
    if k == 0: raise ValueError("Gamma and H do not cross after the end of inflation")
    # Gamma/H rises through 1 for every law used here, so log(Gamma/H) is increasing across the crossing
    Hc = math.exp(np.interp(0.0, np.log(G[k-1:k+1]/H[k-1:k+1]), np.log(H[k-1:k+1])))
    T_rh = math.sqrt(Hc*H_INF*MP/CH)
    T_late = (q*rr*MP**2*H_INF**2)**0.25
    s_a3 = (2*math.pi**2/45)*GSTAR*T_late**3*a**3
    return T_rh, na3*H_INF**3/s_a3

def yield_at(T_rh, p):
    """plasma yield per unit alpha for a history reheating at T_rh. In radiation-like expansion
    H = a^-2, so Gamma = g0 a^p meets H where Gamma_c = g0^(2/(p+2)); bracket g0 around that."""
    Gc = CH*T_rh**2/(MP*H_INF)
    g_est = Gc**((p + 2)/2.0)
    f = lambda lg: math.log(history(math.exp(lg), p)[0]/T_rh)
    lg = brentq(f, math.log(g_est/20), math.log(min(0.9, 20*g_est)), xtol=1e-6)   # g0 < 1: Gamma starts below H
    return history(math.exp(lg), p)[1]

if __name__ == "__main__":
    M_half = M_EXP*2**-0.4
    Y_half = RHO_OVER_S/(2*M_half)
    names = {-1: "fermions (Gamma ~ 1/a)", 0: "constant rate", 1: "bosons, trilinear (Gamma ~ a)"}
    cross = {}
    for p in (-1, 0, 1):
        Ts = np.geomspace(1e13, 4e14, 12)
        ys = np.array([yield_at(T, p) for T in Ts])
        for lab, al in (("metric", ALPHA), ("Einstein-Cartan", ALPHA_EC)):
            T = math.exp(np.interp(math.log(Y_half/al), np.log(ys), np.log(Ts)))
            cross[(p, lab)] = T
        r = yield_at(1.87e13, p)/yield_at(1.87e13, 0)
        print("   %-32s plasma makes half at %.2e GeV (metric), %.2e GeV (Einstein-Cartan); yield at 1.87e13 GeV %.2f x constant rate"
              % (names[p], cross[(p, "metric")], cross[(p, "Einstein-Cartan")], r))
    check(abs(cross[(0, "metric")]/8.97e13 - 1) < 0.02, "plant: a constant rate reproduces v6_mass_smooth.py --w13 (8.97e13 GeV)")
    check(cross[(-1, "metric")] < cross[(0, "metric")] < cross[(1, "metric")], "a rate that is faster early heats the plasma more: fermions lowest, bosons highest")
    sys.exit(1 if fails else 0)
