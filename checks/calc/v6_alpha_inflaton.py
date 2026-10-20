#!/usr/bin/env python3
"""The inflation companion of Separate Ways v6, Section 4: a worked inflaton.

Sections 2 and 3 of that paper assume an inflaton that sits at rest on its plateau at the bang, stays there through the
radiation era before inflation, supplies the perturbations (H_inf = 1e12 GeV, eps = 1.0e-6 at 58
e-folds) and reheats afterwards. This script builds one and reports where it stops working.

1. A T-model alpha-attractor (Kallosh and Linde 2013), V = V0 tanh^2 u, u = phi/(sqrt(6 alpha) M_P),
   with alpha fixed by eps = 1.0e-6 at N = 58. Slow roll in u: eps = 4/(3 alpha sinh^2 2u) and
   N = (3 alpha/4)(cosh 2u - cosh 2u_end), with eps(u_end) = 1.
2. Matter couples to the bounded field chi = chi_p tanh u, chi_p = sqrt(6 alpha) M_P, as in the
   superconformal construction (Kallosh, Linde and Roest 2013): a Yukawa y chi psibar psi to a heavy
   Dirac fermion psi that decays to Standard Model particles. The plasma before inflation pushes the
   inflaton with dF/dphi = (T^2/pi^2) K_F(m/T) 2 m dm/dphi, m = y chi, K_F(z) = int x^2/(E(e^E + 1)),
   and dm/dphi = y sech^2 u is small on the plateau. The field starts at rest at the bang, at the
   value where today's Hubble scale leaves the horizon (the start nearest the minimum, where the
   suppression is weakest), and is followed from T = M_P with the integrator of v6_trilinear_bath.py,
   once with the plasma and once without; the difference is the plasma's doing.
3. Reheating: Gamma = y^2 m_phi/(8 pi) for a Dirac pair, T_rh where H(T_rh) = Gamma. The fermion's
   one-loop (Coleman-Weinberg) energy, -(4/64 pi^2) m^4 [ln(m^2/mu^2) - 3/2], with mu = y chi_p, is
   3 (y chi_p)^4/(32 pi^2) on the plateau and zero at the minimum. It rises with the field, so it adds
   to the plateau, and the tree-level height must make up the rest; the coupling is bounded by asking
   that this energy stay below a tenth of the plateau's, and it equals the plateau at the largest y.
4. The N1 mass at those reheating temperatures, from the perturbative-reheating table for a quadratic
   minimum (v6_mass_smooth.tsv), which is the T-model's minimum.
Plants: the slow-roll solution must match the large-N forms r = 12 alpha/N^2, n_s = 1 - 2/N; K_F(0)
must equal pi^2/12; without the plateau's suppression (chi = phi, fermions light) the same coupling
must move the field by more than M_P/2; and at 1.01 times the largest y the tree-level height must
turn negative.
"""
import math, os, sys
from scipy.integrate import quad
from scipy.optimize import brentq
from v6_bath import MP, CH
from v6_trilinear_bath import follow, RHO, V_PLAT, T_MEET

fails = 0
def check(ok, what):
    global fails
    print("   %-82s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)

H_INF, N_STAR, EPS_STAR = 1e12, 58.0, 1.0e-6

def at_N(alpha, N=N_STAR):
    c_end = math.sqrt(1 + 4/(3*alpha))
    u = math.acosh(4*N/(3*alpha) + c_end)/2
    eps = 4/(3*alpha*math.sinh(2*u)**2)
    t, s2 = math.tanh(u), 1/math.cosh(u)**2
    eta = (2/(6*alpha))*(s2**2 - 2*t*t*s2)/(t*t)
    return u, eps, eta

def K_F(z):
    if z < 1e-6: return math.pi**2/12
    if z > 600: return 0.0
    def f(x):
        E = math.sqrt(x*x + z*z)
        return x*x*math.exp(-E)/(E*(1 + math.exp(-E)))
    return quad(f, 0, math.inf, limit=200)[0]

if __name__ == "__main__":
    alpha = brentq(lambda a: at_N(a)[1] - EPS_STAR, 1e-5, 1.0)
    u58, eps, eta = at_N(alpha)
    r, ns = 16*eps, 1 - 6*eps + 2*eta
    chi_p = math.sqrt(6*alpha)*MP
    V0 = V_PLAT/math.tanh(u58)**2                  # the plateau at N = 58 carries H_inf = 1e12 GeV
    m_phi = math.sqrt(V0/(3*alpha))/MP
    phi58 = u58*chi_p
    sup58 = 1/math.cosh(u58)**2
    print("=== the T-model alpha-attractor with eps = 1.0e-6 at N = 58, H_inf = 1e12 GeV ===")
    print("   alpha = %.4f, r = %.2e, n_s = %.4f, inflaton mass at the minimum %.2e GeV" % (alpha, r, ns, m_phi))
    print("   sqrt(6 alpha) M_P = %.2f M_P: sech^2 u falls by e^-2 for each such step out along the plateau" % math.sqrt(6*alpha))
    # Section 5 of the paper: in de Sitter, k^3/(6 pi^2) modes per comoving volume lie below k and modes leave
    # at k = aH, so H^4/(2 pi^2) leave per physical volume per time, 2/(3 pi) per Hubble volume per Hubble time
    print("   de Sitter: modes leave the horizon at 2/(3 pi) = %.2f per Hubble volume per Hubble time per degree of freedom"
          % (2/(3*math.pi)))
    u_end = math.asinh(math.sqrt(4/(3*alpha)))/2
    print("   inflation ends (eps = 1) at phi = %.2f M_P" % (u_end*math.sqrt(6*alpha)))
    print("   today's Hubble scale leaves at phi = %.2f M_P (u = %.2f); chi_p = sqrt(6 alpha) M_P = %.2e GeV;"
          " sech^2 u there = %.1e" % (phi58/MP, u58, chi_p, sup58))
    # n_s against Planck 2018 VI (0.965 +- 0.004) and ACT DR6 with lensing and DESI BAO (0.974 +- 0.003), both abstracts
    print("   n_s = %.4f is %.1f sigma from Planck's 0.965 +- 0.004 and %.1f sigma below ACT DR6's 0.974 +- 0.003"
          % (ns, abs(ns - 0.965)/0.004, (0.974 - ns)/0.003))
    check(abs(r/(12*alpha/N_STAR**2) - 1) < 0.02 and abs(ns - (1 - 2/N_STAR)) < 2e-3,
          "plant: the slow-roll solution matches r = 12 alpha/N^2 and n_s = 1 - 2/N")
    check(abs(K_F(0.0) - math.pi**2/12) < 1e-9 and abs(K_F(1e-3)/(math.pi**2/12) - 1) < 1e-3,
          "plant: K_F(0) = pi^2/12, and K_F is continuous there")

    print("\n=== reheating through y chi psibar psi, and the one-loop bound on y ===")
    T_per_y = math.sqrt(m_phi*MP/(8*math.pi*CH))
    # The fermion's one-loop potential, -(4/64 pi^2) m^4 [ln(m^2/mu^2) - 3/2] with m = y chi_p x, x = tanh u, is
    # -c x^4 ln x^2 plus a quartic that a counterterm absorbs. Fixing the counterterms so that the plateau height and
    # the mass at the minimum keep their values leaves V = V0 x^2 + c x^4 |ln x^2|, c = (y chi_p)^4/16 pi^2, which
    # vanishes at both ends. dV/dx = 2x [V0 - c x^2 (2 ln x^2 + 1)] and x^2 (2 ln x^2 + 1) <= 1 on (0, 1], so the
    # potential rises all the way to its plateau only while c < V0.
    cfun = lambda y: (y*chi_p)**4/(16*math.pi**2)
    y_mono = brentq(lambda y: cfun(y) - V0, 1e-6, 1.0)
    y_tenth = brentq(lambda y: cfun(y) - 0.1*V0, 1e-6, 1.0)
    def rises(c, n=20000):
        xs_ = [(i + 1)/n for i in range(n)]
        return all(2*V0*x - c*(4*x**3*math.log(x*x) + 2*x**3) > 0 for x in xs_)
    print("   T_rh = y x %.2e GeV (Gamma = y^2 m_phi/8 pi = H), with m_phi fixed by the counterterms" % T_per_y)
    print("   the loop adds c x^4 |ln x^2|, c = (y chi_p)^4/16 pi^2; c is a tenth of V0 at y = %.4f, T_rh = %.2e GeV"
          % (y_tenth, y_tenth*T_per_y))
    print("   and equals V0 at y = %.4f, T_rh = %.2e GeV, beyond which the potential no longer rises to its plateau"
          % (y_mono, y_mono*T_per_y))
    check(rises(0.99*V0) and not rises(1.01*V0), "plant: the potential rises to the plateau for c = 0.99 V0 and not for 1.01 V0")
    y_equal = y_mono

    print("\n=== the plasma before inflation, at the largest y, from rest at phi(N = 58) ===")
    y = y_equal
    sech2 = lambda p: 1/math.cosh(p/chi_p)**2
    chi = lambda p: chi_p*math.tanh(p/chi_p)
    Vp = lambda p: 2*V0*math.tanh(p/chi_p)*sech2(p)/chi_p
    def F_plasma(p, T, suppressed=True):
        if suppressed:
            m = y*chi(p); dm = y*sech2(p)
            return (T*T/math.pi**2)*K_F(m/T)*2*m*dm
        return (T*T/math.pi**2)*K_F(0.0)*2*(y*p)*y              # chi = phi, fermions taken light
    T_end = T_MEET*1e-3
    both = follow(lambda p, T: Vp(p) + F_plasma(p, T), T_end, plateau=True, phi0=phi58, T0=MP)
    tree = follow(lambda p, T: Vp(p), T_end, plateau=True, phi0=phi58, T0=MP)
    d_pl = both.y[0][-1] - tree.y[0][-1]
    dN_dphi = math.sqrt(6*alpha)*math.sinh(2*u58)/(4*MP)
    curv = 0.0
    for i in range(401):
        N = i*both.t[-1]/400; T = MP*math.exp(-N); p = phi58 + both.sol(N)[0]
        H2 = (RHO(T) + V_PLAT)/(3*MP**2)
        dF = (F_plasma(p*1.0005, T) - F_plasma(p*0.9995, T))/(0.001*p)
        curv = max(curv, abs(dF)/H2)
    print("   the plasma moves the field by %.1e M_P, which changes the e-folds of inflation by %.1e;"
          " largest |dF/dphi|/H^2 %.1e" % (d_pl/MP, d_pl*dN_dphi, curv))
    print("   (followed to T = %.1e GeV, about %.0f e-folds into inflation, where the field's own slope has moved"
          " it %.1e M_P)" % (T_end, math.log(T_MEET/T_end), tree.y[0][-1]/MP))
    # the bang may sit anywhere further out on the plateau, so what matters is a shift far below one e-fold
    pre = follow(lambda p, T: Vp(p), T_MEET, plateau=True, phi0=phi58, T0=MP)
    print("   from the bang to the start of inflation (T = %.1e GeV) the field's own slope moves it %.1e M_P, %.2f e-folds"
          % (T_MEET, pre.y[0][-1]/MP, pre.y[0][-1]*dN_dphi))
    check(abs(d_pl*dN_dphi) < 0.1, "the plasma changes the number of e-folds by less than a tenth")
    check(curv < 1e-2, "and gives the inflaton no curvature comparable with H^2")
    bare = follow(lambda p, T: Vp(p) + F_plasma(p, T, suppressed=False), T_MEET, plateau=True, phi0=phi58, T0=MP)
    print("   without the suppression (chi = phi, fermions light) the field ends at %.2f M_P from %.2f"
          % ((phi58 + bare.y[0][-1])/MP, phi58/MP))
    check(abs(bare.y[0][-1]) > 0.5*MP, "plant: without the plateau's suppression the same coupling moves it more than M_P/2")

    print("\n=== the trilinear Higgs coupling of Section 4.1, mu chi H^dag H, moved onto the bounded field ===")
    # four real scalars of mass^2 mu chi: one-loop energy (4/64 pi^2) m^4 [ln(m^2/mu_R^2) - 3/2], mu_R^2 = mu chi_p,
    # is 6 (mu chi_p)^2/(64 pi^2) on the plateau; its slope there carries the factor sech^2 u
    mu_max = math.sqrt(V_PLAT*64*math.pi**2/6)/chi_p
    T_tri = mu_max*math.sqrt(MP/(8*math.pi*m_phi*CH))
    T_inst = (30*V_PLAT/(math.pi**2*106.75))**0.25
    print("   its one-loop energy reaches the plateau's at mu = %.1e GeV, where it would reheat perturbatively to %.1e GeV,"
          " above instant reheating, %.1e GeV" % (mu_max, T_tri, T_inst))
    # what instant reheating needs, and what the Higgs loop then does on the plateau at phi(N = 58): its energy there,
    # 6 (mu chi_p)^2/64 pi^2, and its slope, -(mu^2 chi_p/8 pi^2) sech^2 u, against the tree-level slope
    mu_inst = T_inst/math.sqrt(MP/(8*math.pi*m_phi*CH))
    e_frac = 6*(mu_inst*chi_p)**2/(64*math.pi**2)/V_PLAT
    s_frac = (mu_inst**2*chi_p/(8*math.pi**2))*sup58/Vp(phi58)
    from v6_trilinear_bath import K as K_B
    def F_tri(p, T):
        m2 = max(mu_inst*chi(p), 0.0)
        return mu_inst*sech2(p)*T*T*K_B(math.sqrt(m2)/T)/math.pi**2
    tri = follow(lambda p, T: Vp(p) + F_tri(p, T), T_end, plateau=True, phi0=phi58, T0=MP)
    d_tri = tri.y[0][-1] - tree.y[0][-1]
    print("   instant reheating needs mu = %.1e GeV; there the Higgs loop's energy on the plateau is %.1f per cent of the"
          " plateau's and its slope %.1f per cent of the tree-level slope at phi(N = 58)" % (mu_inst, 100*e_frac, 100*s_frac))
    print("   and the plasma before inflation moves the field by %.1e M_P, %.1e e-folds" % (d_tri/MP, d_tri*dN_dphi))
    print("   so the 3e13 GeV cap of Section 4.1 does not apply on chi; the Higgs is still tachyonic when chi < 0 after inflation")
    check(e_frac < 0.1 and s_frac < 0.1 and abs(d_tri*dN_dphi) < 0.1,
          "instant reheating through mu chi H^dag H leaves the plateau's energy, slope and e-folds within 10 per cent")
    print("\n=== the N1 mass at those reheating temperatures (quadratic minimum, perturbative reheating) ===")
    import numpy as np
    tsv = os.path.join(os.path.dirname(os.path.abspath(__file__)), "v6_mass_smooth.tsv")
    rows = [l.rstrip("\n").split("\t") for l in open(tsv) if l.startswith("with plasma\t")]
    Ts = np.array([float(x[1]) for x in rows]); Ms = np.array([float(x[3]) for x in rows]); o = np.argsort(Ts)
    Ts, Ms = Ts[o], Ms[o]
    for lab, yy in (("c at a tenth of V0", y_tenth), ("c at V0", y_equal)):
        T = yy*T_per_y
        M = math.exp(np.interp(math.log(T), np.log(Ts), np.log(Ms)))
        print("   %s: T_rh = %.2e GeV, M = %.0f PeV" % (lab, T, M/1e6))
    print("\n%s" % ("all checks pass" if fails == 0 else "%d FAILED" % fails))
    sys.exit(1 if fails else 0)
