#!/usr/bin/env python3
"""The inflation companion of Separate Ways v6, section 3.1 and Appendix A: which state the dark-matter modes bring from the bang.

With inflation between a CPT-symmetric bang and the radiation era, the N1 modes that become the
dark matter cross twice: once in the mirror's radiation era and once in ours, with a long stretch
between (the bang, both inflations) in which they are massless. Each half-crossing, started from the
massless vacuum, leaves occupation n (n_vac, from v6_reheating_production.py). Between the two
halves the massless levels pick up a relative phase 2 Phi, with Phi of order x e^N, so:

1. From the massless vacuum at the bang the occupation is n for every Phi.
2. The full in-vacuum gives P = 4 n (1 - n) sin^2(Phi + phi0) (Stueckelberg oscillation).
3. The lemma of Separate Ways v6, section 4, bounds every CPT-invariant state by n* = (1 - sqrt(1 - P))/2. Averaged
   over the rapidly varying phase, n* -> (1 - (2/pi) E(4 n (1 - n)))/2 < n. The states that reach
   it need a population p in the massless excited level at the bang, tuned to Phi.
4. Any population p at the bang is N1 radiation there. Its energy relative to the plasma is
   (30/(pi^4 g*)) kappa^4 Int x^3 p dx, kappa = k_DM/(aT) at the bang =
   sqrt(pi M sqrt(g*/15)/mu) e^N (rho_end/rho_rh)^(1/12).

Plants: (1) and (2) are checked against transfer matrices from the exact solver with the massless
stretch inserted; the elliptic average against a direct phase average; a planted error (dropping
the (1 - n) factor) must fail.
"""
import math, sys
import numpy as np
from scipy.integrate import solve_ivp, quad
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from v6_reheating_production import Profile, adiab, n_exact, n_fourier, n_star_avg, production_integral, SX, SZ
fails = 0
def check(ok, what):
    global fails
    print("   %-72s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)

MU, GSTAR, MP = 5.9653e18, 106.75, 2.435e18
CH = math.sqrt(math.pi**2*GSTAR/90)

def evolve(mfun, x, t0, t1, psi0):
    f = lambda t, y: -1j*((mfun(t)*SZ + x*SX) @ y)
    return solve_ivp(f, (t0, t1), psi0.astype(complex), method='DOP853', rtol=1e-11, atol=1e-13).y[:, -1]

def transfer(x, prof):
    """U[i, j] = <late adiabatic i| U |massless j>, i, j = 0 ground, 1 excited, for our half."""
    t_s = min(-5.0, prof.t_e - 1.0)
    while prof.m(t_s) > 1e-4*x: t_s *= 2.0
    t_e = max(40.0, 6.0*x, 4.0*prof.t_rh + 20.0)
    g0, e0 = adiab(prof.m(t_s), x)
    gT, eT = adiab(prof.m(t_e), x)
    U = np.zeros((2, 2), complex)
    for j, v in enumerate((g0, e0)):
        out = evolve(prof.m, x, t_s, t_e, v)
        U[0, j], U[1, j] = np.vdot(gT, out), np.vdot(eT, out)
    return U, t_s, t_e

def mirror_transfer(x, prof, t_s, t_e):
    """V[i, j] = <massless i| U |far-past adiabatic j> for the mirror half, m_mirror(t) = -m(-t)."""
    mm = lambda t: -prof.m(-t)
    gP, eP = adiab(mm(-t_e), x)
    g0, e0 = adiab(mm(-t_s), x)
    V = np.zeros((2, 2), complex)
    for j, v in enumerate((gP, eP)):
        out = evolve(mm, x, -t_e, -t_s, v)
        V[0, j], V[1, j] = np.vdot(g0, out), np.vdot(e0, out)
    return V

if __name__ == "__main__":
    print("=== two half-crossings with a massless stretch between them ===")
    prof = Profile(t_rh=0.5)
    for x in (0.6, 1.2):
        U, t_s, t_e = transfer(x, prof)
        V = mirror_transfer(x, prof, t_s, t_e)
        n = abs(U[1, 0])**2
        L = np.linspace(0, math.pi/x, 241)[:-1]            # massless stretch: phases e^{-/+ i x L}
        P = np.array([abs((U @ np.diag([np.exp(1j*x*l), np.exp(-1j*x*l)]) @ V)[1, 0])**2 for l in L])
        nstar = (1 - np.sqrt(np.clip(1 - P, 0, 1)))/2
        print("   x = %.1f: n_vac = %.5f; in-vacuum P from %.2e to %.5f (4n(1-n) = %.5f);"
              % (x, n, P.min(), P.max(), 4*n*(1 - n)))
        print("            phase-averaged n* = %.6f, elliptic formula %.6f, n_vac/2 = %.6f" % (nstar.mean(), n_star_avg(n), n/2))
        check(abs(P.max() - 4*n*(1 - n)) < 1e-4 and P.min() < 1e-4, "x = %.1f: in-vacuum oscillates between 0 and 4n(1-n)" % x)
        check(abs(nstar.mean()/n_star_avg(n) - 1) < 1e-3, "x = %.1f: phase average of n* equals the elliptic formula" % x)
        check(abs(nstar.mean()/(1 - math.sqrt(1 - 4*n))*2 - 1) > 1e-2, "plant: dropping the (1 - n) factor fails")
        check(nstar.max() <= n + 1e-9, "x = %.1f: no phase lets the bound exceed the massless vacuum's n" % x)

    print("\n=== the sudden limit: massless vacuum at the start of the radiation era ===")
    I_s = production_integral(Profile(t_rh=1e-4))[0]
    I_STAR = 0.0127596673634
    print("   I = %.4f, C_M = I/pi^1.5 = %.4f, M = 491.6 (I_*/I)^(2/5) = %.1f PeV, I/I_BNNP(0.06305) = %.3f"
          % (I_s, I_s/math.pi**1.5, 491.6*(I_STAR/I_s)**0.4, I_s/0.06305))
    check(abs(491.6*(I_STAR/I_s)**0.4 - 257.6) < 0.3, "sudden limit gives 257.6 PeV")

    print("\n=== the ceiling over all CPT-invariant states, against the massless vacuum ===")
    rows = []
    for tr in (1e-3, 0.1, 0.374, 1.0, 2.0, 4.0):
        pr = Profile(t_rh=tr)
        Iv = production_integral(pr)[0]
        Ic = production_integral(pr, transform=n_star_avg)[0]
        rows.append((tr, Iv, Ic))
        print("   t_rh = %6.3f: I_vac = %.5f, I_ceiling = %.5f, mass ratio (I_vac/I_ceil)^(2/5) = %.3f" % (tr, Iv, Ic, (Iv/Ic)**0.4))
    check(all(r[2] < r[1] for r in rows), "the ceiling state always makes fewer particles than the massless vacuum")

    print("\n=== what a population at the bang weighs, against the plasma ===")
    def kappa(M, N, H_inf, T_rh):
        rho_end, rho_rh = 3*H_inf**2*MP**2, math.pi**2/30*GSTAR*T_rh**4
        return math.sqrt(math.pi*M*math.sqrt(GSTAR/15)/MU)*math.exp(N)*(rho_end/rho_rh)**(1/12)
    def popt_avg(n):
        """phase-averaged population of the bound-saturating state: p = sin^2(arctan(|c| tan 2a)/2)."""
        a = math.asin(math.sqrt(min(max(n, 0.0), 0.5)))
        th = np.linspace(0, math.pi, 2001)[:-1]
        return float(np.mean(np.sin(0.5*np.arctan(np.abs(np.cos(th))*math.tan(2*a)))**2))
    pref = 30/(math.pi**4*GSTAR)
    from scipy.optimize import brentq

    def least_energy_up(xs, ns, target):
        """Least Int x^3 p dx over populations 0 <= p <= 1 with random phases that raise Int x^2 n dx
        by target. A unit of p at x adds x^2 (1 - 2n) dx to production and x^3 dx to the energy, so
        the optimum fills bins in order of (1 - 2n)/x (fractional knapsack, exact on the grid)."""
        dx = np.gradient(xs)
        gain, energy = xs**2*(1 - 2*ns)*dx, xs**3*dx
        acc = E = 0.0
        for i in np.argsort(-(gain/energy)):
            if gain[i] <= 0: break
            if acc + gain[i] >= target: return E + energy[i]*(target - acc)/gain[i]
            acc += gain[i]; E += energy[i]
        return np.inf

    def least_energy_tuned(xs, ns, target, sign):
        """Least Int x^3 p dx with phases tuned to the history that changes production by target, raising it
        (sign = +1) or lowering it (sign = -1); see least_energy_down for the construction."""
        dx = np.gradient(xs)
        cth = np.abs(np.cos(np.linspace(0, math.pi, 97)[:-1]))[None, :, None]
        pg = np.concatenate([[0.0], np.geomspace(1e-14, 0.5, 400)])[None, None, :]
        n3 = ns[:, None, None]
        chg = 2*cth*np.sqrt(pg*(1 - pg)*n3*(1 - n3)) + sign*pg*(1 - 2*n3)
        X = xs[:, None, None]
        def at(lam):
            k = np.argmax(X**2*chg - lam*X**3*pg, axis=2)
            r = np.take_along_axis(chg, k[..., None], 2)[..., 0].mean(axis=1)
            pp = pg[0, 0, k].mean(axis=1)
            return np.sum(xs**2*r*dx), np.sum(xs**3*pp*dx)
        lo, hi = 1e-8, 1e4
        for _ in range(200):
            mid = math.sqrt(lo*hi)
            if at(mid)[0] > target: lo = mid
            else: hi = mid
        return at(lo)[1]

    def least_energy_down(xs, ns, target):
        """Least Int x^3 p dx over populations with phases tuned to the expansion history that lower
        Int x^2 n dx by target. With the best sign, a population p at a mode whose interference phase
        has |cos| = c gives n - n(p) = 2 c sqrt(p (1-p) n (1-n)) - p (1 - 2n): first order in the
        amplitude. Phases are uniform across any interval of x, so average over c = |cos theta|, and
        choose p per (x, theta) to maximise x^2 (n - n(p)) - lam x^3 p; bisect on lam."""
        dx = np.gradient(xs)
        cth = np.abs(np.cos(np.linspace(0, math.pi, 97)[:-1]))[None, :, None]
        pg = np.concatenate([[0.0], np.geomspace(1e-14, 0.5, 400)])[None, None, :]
        n3 = ns[:, None, None]
        red = 2*cth*np.sqrt(pg*(1 - pg)*n3*(1 - n3)) - pg*(1 - 2*n3)
        X = xs[:, None, None]
        def at(lam):
            k = np.argmax(X**2*red - lam*X**3*pg, axis=2)
            r = np.take_along_axis(red, k[..., None], 2)[..., 0].mean(axis=1)
            pp = pg[0, 0, k].mean(axis=1)
            return np.sum(xs**2*r*dx), np.sum(xs**3*pp*dx)
        lo, hi = 1e-8, 1e4
        for _ in range(200):
            mid = math.sqrt(lo*hi)
            if at(mid)[0] > target: lo = mid
            else: hi = mid
        R, E = at(lo)
        return E, R

    print("   N1 at the bang has energy (30/(pi^4 g*)) kappa^4 Int x^3 p dx relative to the plasma. The e-folds")
    print("   above which that exceeds the plasma, for the state that saturates the bound and for the least")
    print("   population that moves production by ten per cent up (random phases) or down (tuned phases):")
    last = None
    for T_rh, M in ((1.87e13, 4.134e8), (5e13, 2.784e8)):       # masses from v6_mass_with_bath.py (central)
        H_rh = CH*T_rh**2/MP
        for H_inf in (1e10, 1e12, 1e13):
            if H_inf < H_rh: continue
            t_rh = math.sqrt(M/H_rh); t_e = min((8*t_rh*M/H_inf)**(1/3) - t_rh, t_rh)
            pr = Profile(t_rh=t_rh, t_e=t_e)
            xs = np.linspace(0.01, 6.0, 240); ns = np.array([n_exact(x, pr) for x in xs])
            xt = np.geomspace(6.0, 400.0, 120); nt = np.array([n_fourier(x, pr) for x in xt])
            E_ceil = np.trapezoid(xs**3*[popt_avg(v) for v in ns], xs) + np.trapezoid(xt**3*nt/2, xt)
            Int2 = np.trapezoid(xs**2*ns, xs)
            E_up = least_energy_up(xs, ns, 0.1*Int2)
            E_dn, R_dn = least_energy_down(xs, ns, 0.1*Int2)
            E_tu = least_energy_tuned(xs, ns, 0.1*Int2, +1)
            Ns = [brentq(lambda N: pref*kappa(M, N, H_inf, T_rh)**4*E - 1, 0, 80) for E in (E_ceil, E_up, E_dn, E_tu)]
            print("      T_rh %.2e, H_inf %.0e: Int x^3 p = %.1e, %.1e, %.1e -> N = %.1f (bound), %.1f (up), %.1f (down)"
                  % (T_rh, H_inf, E_ceil, E_up, E_dn, Ns[0], Ns[1], Ns[2]))
            print("         tuned phases raising production ten per cent: Int x^3 p = %.1e -> N = %.1f;"
                  " at N = 15 the largest allowed change is %.1f%% (10%% e^{-2(15 - N)})" % (E_tu, Ns[3], 10*math.exp(-2*(15 - max(Ns[2], Ns[3])))))
            last = (M, H_inf, T_rh, E_up, E_dn, R_dn, Int2)
    M, H_inf, T_rh, E_up, E_dn, R_dn, Int2 = last
    check(abs(R_dn/(0.1*Int2) - 1) < 0.02, "the tuned-phase optimum reaches the ten per cent target")
    xb = math.exp(12.6)*((3*1e24*MP**2)/(math.pi**2/30*GSTAR*1.87e13**4))**(1/12)
    print("   at N = 12.6, H_inf = 1e12, T_rh = 1.87e13: x at the bang = %.2e x; P = exp(-pi %.1e); 1/(16 x^4) = %.1e"
          % (xb, xb**2, 1/(16*xb**4)))
    print("   without inflation kappa = sqrt(pi M sqrt(g*/15)/mu) = %.2e (M = 491.6 PeV)" % math.sqrt(math.pi*4.916e8*math.sqrt(GSTAR/15)/MU))
    k60 = kappa(M, 60, 1e12, 1.87e13)
    print("   with 60 e-folds the dark-matter modes sit at k/(aT) = %.1e at the bang" % k60)
    check(pref*kappa(M, 14, H_inf, T_rh)**4*min(E_up, E_dn) > 1, "at 14 e-folds a ten per cent change either way outweighs the plasma (last case)")
    check(pref*kappa(M, 8, H_inf, T_rh)**4*E_up < 1, "plant: at 8 e-folds it does not, so the test can fail")
    sys.exit(1 if fails else 0)
