#!/usr/bin/env python3
"""The inflation companion of Separate Ways v6, section 3.2: the N1 mass against the reheating temperature for perturbative
reheating, the inflaton decaying at a constant rate Gamma = H(T_rh) (v6_reheating_smooth.py).

The graviton-exchange yield is integrated along the actual temperature history, which rises to a
maximum above T_rh during reheating and then falls:
    d(n a^3)/dt = alpha T^8 a^3 / M_P^4,  T^4 = 30 rho_r/(pi^2 g*),
and divided by the late-time entropy. So no reheating "boost" is put in by hand. The mass solves
    M (K M^(3/2) I(M) + Y_plasma) = rho_DM/s0,
with K fixed by the half-angle normalisation as in v6_mass_with_bath.py. Variants: the expansion
alone, with the plasma, and the two-crossing bound over all CPT-invariant states (with the plasma).
Writes v6_mass_smooth.tsv. Checks: for a near-instant decay the yield approaches the sudden formula
alpha (T/M_P)^3/C at the highest temperature reached (with Gamma = H(T_rh) and a decay of order a
Hubble time, half the inflaton is still undecayed at H = Gamma, so T_max < T_rh there), and the
expansion-only mass matches v6_reheating_smooth.py. For slow decay the yield is 0.53 of the sudden
formula at T_rh, consistent with Garcia et al. 2017 once their T_RH (40 instead of 90 in the
Gamma-T relation, T_RH 0.816 of ours) is converted.
"""
import math, sys
import numpy as np
from multiprocessing import Pool
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from v6_reheating_production import n_exact, n_pert, n_star_avg
from v6_reheating_smooth import Smooth
from v6_bath import Y_grav, RHO_OVER_S, C as C_SUDDEN
fails = 0
def check(ok, what):
    global fails
    print("   %-70s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)

MP, GSTAR = 2.435e18, 106.75
CH = math.sqrt(math.pi**2*GSTAR/90)
ALPHA = 1.2027e-2
# --ec: add the contact term of Einstein-Cartan gravity with minimally coupled fermions, 7.4185e-2
# (v6_grav_rate.py), to the graviton-exchange rate; only the "with plasma" variant is solved
EC = "--ec" in sys.argv
if EC: ALPHA = 1.2027e-2 + 7.4185e-2
M_STAR, I_STAR = 4.916e8, 0.0127596673634
K = RHO_OVER_S/(M_STAR**2.5*I_STAR)

W = 1.0/3.0 if "--w13" in sys.argv else 0.0     # inflaton fluid: 0 for a quadratic minimum, 1/3 for a quartic

def background(H_inf, Gamma, w=None):
    """a(eta), late a_1, eta_0, and the graviton-exchange yield along the history. Units a_e = 1,
    time in 1/H_inf, energy densities in M_P^2 H_inf^2. The oscillating inflaton is a fluid with
    pressure w rho (w = 0 for a quadratic minimum, 1/3 for a quartic one)."""
    w = W if w is None else w
    g = Gamma/H_inf
    q = 30.0/(math.pi**2*GSTAR)
    def f(t, y):
        a, rp, rr, eta, na3 = y
        H = math.sqrt(max(rp + rr, 0.0)/3.0)
        return [a*H, -3*(1 + w)*H*rp - g*rp, -4*H*rr + g*rp, 1.0/a, (q*rr)**2*a**3]
    t_end = 400.0/g
    ts = np.concatenate([[0.0], np.geomspace(1e-4, t_end, 6000)])
    sol = solve_ivp(f, (0.0, t_end), [1.0, 3.0, 0.0, 0.0, 0.0], t_eval=ts, method="LSODA", rtol=1e-10, atol=1e-30)
    a, rp, rr, eta, na3 = sol.y
    H = np.sqrt((rp + rr)/3.0)
    a1 = a[-1]**2*H[-1]
    eta0 = eta[-1] - a[-1]/a1
    # yield: n a^3 carries alpha H_inf^3 (q rho_r)^2 a^3 dt in these units; s a^3 = (2 pi^2/45) g* T^3 a^3
    na3_phys = ALPHA*na3[-1]*H_inf**3      # T^8/M_P^4 = (q rho_r)^2 H_inf^4 in these units, and dt = dtau/H_inf
    T_late = (q*rr[-1]*MP**2*H_inf**2)**0.25
    s_a3 = (2*math.pi**2/45)*GSTAR*T_late**3*a[-1]**3
    Tmax = (q*rr.max()*MP**2*H_inf**2)**0.25
    return CubicSpline(eta, a), eta[-1], a1, eta0, na3_phys/s_a3, Tmax

def I_of(M, H_inf, Gamma, bg, transform=None):
    tf = (lambda v: v) if transform is None else transform
    spl, eta_last, a1, eta0 = bg[:4]
    prof = Smooth(M/H_inf, spl, eta_last, a1, eta0, Gamma/H_inf)
    xs = np.linspace(1e-2, 6.0, 97)
    ns = tf(np.array([n_exact(x, prof) for x in xs]))
    xt = np.geomspace(6.0, 20.0, 12)
    nt = tf(np.array([n_pert(x, prof) for x in xt]))
    return (np.trapezoid(xs**2*ns, xs) + np.trapezoid(xt**2*nt, xt) + xt[-1]**3*nt[-1]/3.0)/math.sqrt(math.pi)

def job(a):
    tag, T_rh, H_inf = a
    Gamma = CH*T_rh**2/MP
    bg = background(H_inf, Gamma)
    Yb = 0.0 if tag == "expansion only" else bg[4]
    tf = n_star_avg if tag == "bound, with plasma" else None
    f = lambda lnM: math.log(math.exp(lnM)*(K*math.exp(lnM)**1.5*I_of(math.exp(lnM), H_inf, Gamma, bg, tf) + Yb)/RHO_OVER_S)
    try:
        lnM = brentq(f, math.log(1e3), math.log(5e9), xtol=1e-4)
    except ValueError:
        return tag, T_rh, H_inf, float("nan"), float("nan"), bg[4], bg[5]
    M = math.exp(lnM)
    Ye = K*M**1.5*I_of(M, H_inf, Gamma, bg, tf)
    return tag, T_rh, H_inf, M, Yb/(Ye + Yb), bg[4], bg[5]

if __name__ == "__main__":
    # the yield along the history, against the sudden formula
    for T in (1.87e13, 5e13, 5.9e14):
        bg = background(1e12, CH*T**2/MP)
        print("   T_rh %.2e GeV: plasma yield %.3e, sudden formula %.3e, ratio %.2f; T_max = %.2e GeV"
              % (T, bg[4], Y_grav(T, ALPHA), bg[4]/Y_grav(T, ALPHA), bg[5]))
    bgf = background(1e12, CH*5.9e14**2/MP)
    ratio_fast = bgf[4]/Y_grav(bgf[5], ALPHA)    # near-instant decay: compare with the sudden formula at T_max
    Ts = [3e12, 5e12, 8e12, 1.2e13, 1.87e13, 2.5e13, 3.5e13, 5e13, 7e13, 1e14, 1.4e14, 2e14]
    tags = ["with plasma"] if EC else ["expansion only", "with plasma", "bound, with plasma"]
    jobs = [(tag, T, 1e12) for tag in tags for T in Ts]
    with Pool(min(30, len(jobs))) as p:
        res = p.map(job, jobs)
    out = __file__.rsplit("/", 1)[0] + "/v6_mass_smooth%s%s.tsv" % ("_w13" if W else "", "_ec" if EC else "")
    with open(out, "w") as fh:
        fh.write("variant\tT_rh_GeV\tH_inf_GeV\tM_GeV\tplasma_fraction\tY_plasma\tT_max_GeV\n")
        for r in res: fh.write("%s\t%.4e\t%.1e\t%.6e\t%.6f\t%.4e\t%.4e\n" % r)
    print("wrote %s" % out)
    for tag in tags:
        print("\n=== %s (perturbative reheating, w = %.2f, H_inf = 1e12 GeV) ===" % (tag, W))
        print("   T_rh (GeV)    M (PeV)   line (PeV)   from the plasma   T_max (GeV)")
        for r in res:
            if r[0] == tag:
                print("   %.2e     %7.1f     %6.1f      %5.1f%%         %.2e" % (r[1], r[3]/1e6, r[3]/2e6, 100*r[4], r[6]))
    R = {(r[0], r[1]): r[3] for r in res}
    if EC:
        print("   rate %.4e T^8/M_P^4 (graviton exchange plus the minimal Einstein-Cartan contact term)" % ALPHA)
    elif not W:
        check(abs(R[("expansion only", 1.87e13)]/4.508e8 - 1) < 0.01, "expansion-only mass at 1.87e13 GeV matches v6_reheating_smooth.py (450.8 PeV)")
    else:
        e = [R[("expansion only", T)] for T in Ts]
        print("   radiation-like reheating: expansion-only mass from %.1f to %.1f PeV over the grid" % (min(e)/1e6, max(e)/1e6))
        check(max(e)/min(e) - 1 < 0.02, "w = 1/3: the expansion-only mass does not depend on T_rh (within 2 per cent)")
    print("   near-instant decay (T_rh = 5.9e14 GeV): yield / sudden formula at T_max = %.2f" % ratio_fast)
    check(0.7 < ratio_fast <= 1.0, "plant: near-instant decay gives a yield close to the sudden formula at T_max")
    sys.exit(1 if fails else 0)
