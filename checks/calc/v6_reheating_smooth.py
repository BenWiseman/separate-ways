#!/usr/bin/env python3
"""The inflation companion of Separate Ways v6, section 3.2: how much the N1 mass depends on the sudden end of reheating.

The production calculation (v6_reheating_production.py) treats reheating as matter domination that
turns into radiation at one instant, t_rh. Here the inflaton decays at a constant rate Gamma instead:

    d rho_phi/dt = -3 H rho_phi - Gamma rho_phi,   d rho_r/dt = -4 H rho_r + Gamma rho_phi,
    H^2 = (rho_phi + rho_r)/(3 M_P^2),

starting at the end of inflation (rho_phi = 3 H_inf^2 M_P^2, rho_r = 0) and preceded by de Sitter,
matched in a and a'. The reheating temperature is defined the usual way, Gamma = H at T_rh in the
radiation era, which is also where the sudden model switches to radiation, so the two models are
compared at the same T_rh. The late radiation era fixes a_1 and eta_0, and the mode variable is
t = sqrt(M a_1) (eta - eta_0), m(t) = sqrt(M/a_1) a. The mass solves M^(5/2) I(M) = M_*^(5/2) I_*
(expansion alone, as the plasma's yield is the same function of T_rh in both models).
Plant: with Gamma -> large (decay within a fraction of a Hubble time) the smooth model must approach
the sudden one.
"""
import math, sys
import numpy as np
from multiprocessing import Pool
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from v6_reheating_production import Profile, n_exact, n_pert
fails = 0
def check(ok, what):
    global fails
    print("   %-70s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)

MP, GSTAR = 2.435e18, 106.75
CH = math.sqrt(math.pi**2*GSTAR/90)
M_STAR, I_STAR = 4.916e8, 0.0127596673634

def background(H_inf, Gamma):
    """a(eta) after inflation in units a_e = 1, eta in 1/(a_e H_inf); returns spline and late a_1, eta_0."""
    g = Gamma/H_inf
    def f(t, y):
        a, rp, rr, eta = y
        H = math.sqrt(max(rp + rr, 0.0)/3.0)
        return [a*H, -3*H*rp - g*rp, -4*H*rr + g*rp, 1.0/a]
    t_end = 400.0/g
    ts = np.concatenate([[0.0], np.geomspace(1e-4, t_end, 6000)])
    sol = solve_ivp(f, (0.0, t_end), [1.0, 3.0, 0.0, 0.0], t_eval=ts, method="LSODA", rtol=1e-10, atol=1e-14)
    a, rp, rr, eta = sol.y
    H = np.sqrt((rp + rr)/3.0)
    a1 = a[-1]**2*H[-1]                      # da/deta in the radiation era
    eta0 = eta[-1] - a[-1]/a1
    drift = abs(a[-200]**2*H[-200]/a1 - 1)
    return CubicSpline(eta, a), eta[-1], a1, eta0, drift

class Smooth:
    """Mode-profile interface for the solver: m(t), mdot(t), t_e, t_rh."""
    def __init__(self, M_over_H, spl, eta_last, a1, eta0, Gamma_over_H):
        self.k = math.sqrt(M_over_H*a1); self.c = math.sqrt(M_over_H/a1)
        self.spl, self.eta_last, self.a1, self.eta0 = spl, eta_last, a1, eta0
        self.t_e = self.k*(0.0 - eta0)
        # t at which H = Gamma in the extrapolated radiation era: a = a1 (eta - eta0), H = a1/a^2
        self.t_rh = self.k*math.sqrt(1.0/Gamma_over_H)/math.sqrt(a1)*1.0
    def a(self, eta):
        if eta <= 0.0: return 1.0/(1.0 - eta)                     # de Sitter, a' matched at eta = 0
        if eta >= self.eta_last: return self.a1*(eta - self.eta0)
        return float(self.spl(eta))
    def ap(self, eta):
        if eta <= 0.0: return 1.0/(1.0 - eta)**2
        if eta >= self.eta_last: return self.a1
        return float(self.spl(eta, 1))
    def m(self, t): return self.c*self.a(self.eta0 + t/self.k)
    def mdot(self, t): return self.c*self.ap(self.eta0 + t/self.k)/self.k

def I_smooth(M, H_inf, Gamma, bg):
    spl, eta_last, a1, eta0, _ = bg
    prof = Smooth(M/H_inf, spl, eta_last, a1, eta0, Gamma/H_inf)
    xs = np.linspace(1e-2, 6.0, 97)
    ns = np.array([n_exact(x, prof) for x in xs])
    I = np.trapezoid(xs**2*ns, xs)
    xt = np.geomspace(6.0, 20.0, 12)
    nt = np.array([n_pert(x, prof) for x in xt])
    return (I + np.trapezoid(xt**2*nt, xt) + xt[-1]**3*nt[-1]/3.0)/math.sqrt(math.pi)

def I_sudden(M, T_rh, H_inf):
    from v6_reheating_production import production_integral
    H_rh = CH*T_rh**2/MP; t_rh = math.sqrt(M/H_rh); t_e = min((8*t_rh*M/H_inf)**(1/3) - t_rh, t_rh)
    return production_integral(Profile(t_rh=t_rh, t_e=t_e))[0]

def solve(job):
    model, T_rh, H_inf = job
    Gamma = CH*T_rh**2/MP
    if model == "smooth":
        bg = background(H_inf, Gamma)
        f = lambda lnM: 2.5*lnM + math.log(I_smooth(math.exp(lnM), H_inf, Gamma, bg)) - math.log(M_STAR**2.5*I_STAR)
    else:
        f = lambda lnM: 2.5*lnM + math.log(I_sudden(math.exp(lnM), T_rh, H_inf)) - math.log(M_STAR**2.5*I_STAR)
    lnM = brentq(f, math.log(5e7), math.log(5e9), xtol=1e-4)
    return model, T_rh, H_inf, math.exp(lnM)/1e6

if __name__ == "__main__":
    bg = background(1e12, CH*(1.87e13)**2/MP)
    print("background at T_rh = 1.87e13 GeV, H_inf = 1e12 GeV: radiation-era slope settled to %.1e" % bg[4])
    check(bg[4] < 1e-3, "the late-time scale factor is linear in conformal time")
    jobs = [(m, T, 1e12) for T in (1.2e13, 1.87e13, 3e13, 5e13) for m in ("smooth", "sudden")]
    jobs += [("smooth", 1.87e13, 1e13), ("smooth", 5.9e14, 1e12), ("sudden", 5.9e14, 1e12)]
    with Pool(min(len(jobs), 10)) as p:
        res = p.map(solve, jobs)
    R = {(m, T, H): M for m, T, H, M in res}
    print("\n   T_rh (GeV)   sudden (PeV)   smooth decay (PeV)   ratio")
    worst = 0.0
    for T in (1.2e13, 1.87e13, 3e13, 5e13):
        a, b = R[("sudden", T, 1e12)], R[("smooth", T, 1e12)]
        worst = max(worst, abs(b/a - 1))
        print("   %.2e     %8.1f        %8.1f           %.3f" % (T, a, b, b/a))
    print("   smooth decay at H_inf = 1e13 GeV, T_rh = 1.87e13: %.1f PeV" % R[("smooth", 1.87e13, 1e13)])
    print("   largest difference between the two reheating models: %.1f%%" % (100*worst))
    fast = R[("smooth", 5.9e14, 1e12)]/R[("sudden", 5.9e14, 1e12)]
    print("   plant, decay within two Hubble times of the end of inflation (T_rh = 5.9e14 GeV): ratio %.3f" % fast)
    check(abs(fast - 1) < 0.03, "plant: a fast decay reproduces the sudden model")
    check(worst < 0.5, "the two reheating models agree to within a factor 1.5")
    sys.exit(1 if fails else 0)
