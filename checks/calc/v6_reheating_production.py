#!/usr/bin/env python3
"""Gravitational production of the stable right-handed neutrino when inflation precedes the
radiation era, as a function of the reheating temperature. Separate Ways v6, section 3, and its inflation companion,
section 3.

WHY. Without inflation the radiation bang reaches any temperature, and graviton exchange in the
bath then makes N1 without limit (Separate Ways v6, section 5). Any period of inflation removes that, but it
also changes the state in which the N1 modes reach the epoch H = M1, where they are made. This
script follows each mode through inflation, an optional matter-dominated reheating phase and the
radiation era, and returns the production integral I that fixes the mass.

THE MODE EQUATION. A Dirac/Majorana mode of conformal momentum k in conformally flat FRW obeys
    i d/deta (u, v) = [ M a(eta) sigma_z + k sigma_x ] (u, v).
In the radiation era a = a1 (eta - eta0). With t = sqrt(M a1) (eta - eta0) and x = k/sqrt(M a1)
this is i d/dt psi = [ m(t) sigma_z + x sigma_x ] psi with m(t) = t: a two-level crossing at unit
sweep rate whose full sweep leaves P = exp(-pi x^2). Boyle, Finn and Turok's variable is x_sw = sqrt(pi) x, so
the production integral I = (1/pi^2) Int x_sw^2 n dx_sw = pi^(-1/2) Int x^2 n dx.

THE BACKGROUND. Radiation for t >= t_rh (m = t). Matter-dominated reheating for t_e <= t <= t_rh,
m = (t + t_rh)^2/(4 t_rh), matched in a and a' at t_rh. De Sitter for t <= t_e,
m = m_e/(1 - (m'_e/m_e)(t - t_e)). In these units H/M = m'/m^2, so H_rh = M/t_rh^2 and
H_inf = 8 M t_rh/(t_e + t_rh)^3. Instantaneous reheating is t_e = t_rh.

THE STATE. During inflation the relevant modes have m << x: the field is effectively massless and
conformal, so inflation neither excites nor relaxes it. A state fixed at a CPT-symmetric bang
before inflation is carried through with its populations in the massless (sigma_x) basis intact
and their relative phase scrambled by the long conformal time. The massless vacuum, p_e = 0, is
the least-occupied choice; any other member of the CPT family adds p_e (1 - 2 n_vac) >= 0.
"""
import sys, math
import numpy as np
from scipy.integrate import solve_ivp, quad

SX = np.array([[0, 1], [1, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)
fails = 0
def check(ok, what):
    global fails
    print("   %-62s %s" % (what, "yes" if ok else "NO  <-- FAILED"))
    if not ok: fails += 1


class Profile:
    """m(t) for inflation -> (matter) -> radiation, or a supplied function."""
    def __init__(self, t_rh, t_e=None):
        self.t_rh = t_rh
        self.t_e = t_rh if t_e is None else t_e
        te, tr = self.t_e, self.t_rh
        self.me = (te + tr)**2/(4*tr) if te < tr else tr
        self.mpe = (te + tr)/(2*tr) if te < tr else 1.0
    def m(self, t):
        te, tr = self.t_e, self.t_rh
        if t >= tr: return t
        if t >= te: return (t + tr)**2/(4*tr)
        return self.me/(1.0 - (self.mpe/self.me)*(t - te))
    def mdot(self, t):
        te, tr = self.t_e, self.t_rh
        if t >= tr: return 1.0
        if t >= te: return (t + tr)/(2*tr)
        return self.mpe*(self.m(t)/self.me)**2


def adiab(m, x):
    w, v = np.linalg.eigh(m*SZ + x*SX)
    return v[:, 0], v[:, 1]


def n_exact(x, prof, t_start=None, t_end=None):
    """Occupation at late times, starting in the instantaneous ground state deep in inflation."""
    if t_start is None:
        # Start in the instantaneous ground state where the evolution is adiabatic to high order. The
        # state then misses the true one by an amplitude of about m'/(2 x^2), which would interfere
        # with the production amplitude later on, so go back until m' < 1e-6 x^2 (and m < 0.02 x).
        # In de Sitter m' falls as 1/t^2, so this costs a few tens of oscillations at most.
        t_start = min(-5.0, prof.t_e - 1.0)
        while (prof.m(t_start) > 0.02*x or prof.mdot(t_start) > 1e-6*x*x) and t_start > -1e7: t_start *= 2.0
    if t_end is None: t_end = max(40.0, 6.0*x, 4.0*prof.t_rh + 20.0)
    g0, _ = adiab(prof.m(t_start), x)
    f = lambda t, y: -1j*((prof.m(t)*SZ + x*SX) @ y)
    # the inflationary stretch is slow and adiabatic; the crossing near t ~ 0..few is where the work is
    sol = solve_ivp(f, (t_start, t_end), g0.astype(complex), method='DOP853', rtol=1e-10, atol=1e-12)
    _, eT = adiab(prof.m(t_end), x)
    return abs(np.vdot(eT, sol.y[:, -1]))**2


def n_pert(x, prof):
    """First-order adiabatic perturbation theory, accurate when n << 1 (large x):
    c_e = Int dt (x m'/(2 E^2)) exp(2 i Phi),  E = sqrt(m^2 + x^2),  Phi = Int E dt."""
    # integrate from deep inflation to late radiation; coupling ~ x m'/(2E^2)
    ts = []
    t0 = min(-5.0, prof.t_e - 1.0)
    while prof.m(t0) > 1e-6*x: t0 *= 2.0
    T = max(60.0, 10*x)
    # dense grid where the phase moves: step << 1/E
    grid = np.concatenate([np.linspace(t0, prof.t_e, 4000), np.linspace(prof.t_e, prof.t_rh, 4000),
                           np.linspace(prof.t_rh, T, int(400*T*max(1.0, x/20)) + 1)])
    grid = np.unique(grid)
    m = np.array([prof.m(t) for t in grid]); md = np.array([prof.mdot(t) for t in grid])
    E = np.sqrt(m*m + x*x)
    Phi = np.concatenate([[0.0], np.cumsum(0.5*(E[1:] + E[:-1])*np.diff(grid))])
    integrand = x*md/(2*E*E)*np.exp(2j*Phi)
    c = np.sum(0.5*(integrand[1:] + integrand[:-1])*np.diff(grid))
    return abs(c)**2


def n_fourier(x, prof):
    """Large-x tail. When x >> m wherever m' changes, E ~ x and the first-order amplitude is the
    Fourier transform of m'(t) at frequency 2x:  c = (1/2x) Int m'(t) exp(2 i x t) dt,  so
    n = |FT[m'](2x)|^2/(4 x^2). For a pure kink (m' = step) this is 1/(16 x^4); the C^1 junctions of
    the real profile make it fall as x^-6 beyond x ~ 1/t_rh."""
    w = 2.0*x
    te, tr = prof.t_e, prof.t_rh
    # radiation: Int_{tr}^inf e^{iwt} dt = -e^{iw tr}/(iw) = i e^{iw tr}/w
    c = 1j*np.exp(1j*w*tr)/w
    # matter: m' = (t+tr)/(2tr) on [te, tr]; Int (A t + B) e^{iwt} dt = e^{iwt}[(A t + B)/(iw) + A/w^2]
    if te < tr:
        A, B = 1.0/(2*tr), 0.5
        F = lambda t: np.exp(1j*w*t)*((A*t + B)/(1j*w) + A/w**2)
        c += F(tr) - F(te)
    # de Sitter: m' = mpe/(1 - b (t - te))^2, b = mpe/me, on (-inf, te]; substitute u = te - t >= 0
    bb = prof.mpe/prof.me
    g = lambda u: prof.mpe/(1.0 + bb*u)**2
    re = quad(g, 0, np.inf, weight='cos', wvar=w, limlst=200)[0]
    im = quad(g, 0, np.inf, weight='sin', wvar=w, limlst=200)[0]
    # Int g(u) e^{iw(te-u)} du = e^{iw te} (C - i S)
    c += np.exp(1j*w*te)*(re - 1j*im)
    return abs(c)**2/(4*x*x)


def n_star_avg(n):
    """Least occupation over CPT-invariant states when two half-crossings (the mirror's and ours)
    are separated by a long massless stretch: n* = (1 - sqrt(1 - P))/2 with P = 4 n (1 - n) sin^2 Phi,
    averaged over the rapidly varying Phi: (1 - (2/pi) E(4 n (1 - n)))/2, E the complete elliptic
    integral of the second kind (parameter convention). Series q/8 + 3q^2/128 + 5q^3/512 for small q."""
    from scipy.special import ellipe
    n = np.clip(np.asarray(n, float), 0.0, 1.0)
    q = 4.0*n*(1.0 - n)
    return np.where(q < 1e-4, q/8 + 3*q*q/128 + 5*q**3/512, 0.5*(1.0 - (2.0/math.pi)*ellipe(q)))


def production_integral(prof, x_join=None, nx=None, transform=None):
    """I = pi^-1/2 Int x^2 n dx: exact two-level evolution for x <= x_join, Fourier tail beyond.
    The tail formula needs m << x wherever m' varies, and m reaches t_rh there, so the join moves
    out with t_rh (validated against the exact solver in the plants). transform maps the
    massless-vacuum occupation to another state's (e.g. n_star_avg) before integrating."""
    tf = (lambda a: a) if transform is None else transform
    # The join stays at x = 6. There the tail formula matches the exact solver to under 1 per cent
    # for fast transitions, where the tail matters; for slow transitions the occupation beyond 6
    # is below 1e-6, near the exact solver's own floor, and carries under 1 per cent of I (plants).
    if x_join is None: x_join = 6.0
    if nx is None: nx = 97
    xs = np.linspace(1e-2, x_join, nx)     # below 0.01 the integrand x^2 n is under 5e-5 of the total
    ns = tf(np.array([n_exact(x, prof) for x in xs]))
    I_body = np.trapezoid(xs**2*ns, xs)
    xmax = max(60.0, 60.0/max(prof.t_rh, 1e-3))
    xt = np.geomspace(x_join, xmax, 160)
    nt = tf(np.array([n_fourier(x, prof) for x in xt]))
    I_tail = np.trapezoid(xt**2*nt, xt)
    I_rest = xt[-1]**3*nt[-1]/3.0          # beyond xmax the tail falls at least as x^-6
    return (I_body + I_tail + I_rest)/math.sqrt(math.pi), (I_body, I_tail, I_rest)


if __name__ == "__main__":
    print("=== plants: the solver must reproduce the two production integrals of Separate Ways v6, section 3 ===")
    # full sweep from the in-vacuum: P = exp(-pi x^2), I = 1/(4 pi^1.5)
    class Full:
        t_e = -1e9; t_rh = -1e9
        def m(self, t): return t
        def mdot(self, t): return 1.0
    xs = np.linspace(0.01, 2.5, 60)
    Ps = []
    for x in xs:
        g, _ = adiab(-60.0, x)
        f = lambda t, y, x=x: -1j*((t*SZ + x*SX) @ y)
        s = solve_ivp(f, (-60.0, 60.0), g.astype(complex), method='DOP853', rtol=1e-10, atol=1e-12)
        _, e = adiab(60.0, x); Ps.append(abs(np.vdot(e, s.y[:, -1]))**2)
    Ps = np.array(Ps)
    check(np.max(np.abs(Ps - np.exp(-math.pi*xs**2))) < 2e-3, "full sweep returns P = exp(-pi x^2) (to the finite-T error)")
    I_in = np.trapezoid(xs**2*np.exp(-math.pi*xs**2), xs)/math.sqrt(math.pi)
    check(abs(I_in - 1/(4*math.pi**1.5)) < 1e-4, "and I_in = 1/(4 pi^1.5) = %.6f" % (1/(4*math.pi**1.5)))
    nmin = lambda x: (1 - math.sqrt(1 - math.exp(-math.pi*x*x)))/2
    I_min = quad(lambda x: x*x*nmin(x), 0, 12, limit=400)[0]/math.sqrt(math.pi)
    check(abs(I_min - 0.0127596673634) < 1e-9, "the half-angle state gives I_min = %.10f" % I_min)

    print("\n=== the sudden limit: radiation from the crossing point, massless vacuum entering ===")
    prof = Profile(t_rh=1e-4)
    for x in [0.5, 1.0, 3.0, 6.0]:
        ne = n_exact(x, prof); npt = n_pert(x, prof)
        print("   x = %.1f   exact %.6e   first order %.6e   1/(16 x^4) %.3e" % (x, ne, npt, 1/(16*x**4)))
    check(abs(n_pert(6.0, prof)/n_exact(6.0, prof) - 1) < 0.05, "first-order tail matches the exact solution at x = 6")
    for tr in [1e-4, 0.1]:
        pr = Profile(t_rh=tr)
        r = n_fourier(6.0, pr)/n_exact(6.0, pr)
        check(abs(r - 1) < 0.03, "fast transition t_rh = %g: tail formula matches exact at x = 6 (ratio %.4f)" % (tr, r))
    for tr, te in [(0.5, None), (1.0, 0.3), (2.0, 0.2)]:
        I, parts = production_integral(Profile(t_rh=tr, t_e=te))
        share = (parts[1] + parts[2])/(parts[0] + parts[1] + parts[2])
        check(share < 0.01, "slow transition t_rh = %g: tail beyond x = 6 carries %.2f%% of I" % (tr, 100*share))
    # plant for the share test: the sudden limit's x^-4 tail must NOT pass it
    I, parts = production_integral(Profile(t_rh=1e-4))
    share = (parts[1] + parts[2])/(parts[0] + parts[1] + parts[2])
    check(share > 0.01, "plant: the sudden limit's tail carries %.1f%% of I, so the share test can fail" % (100*share))
    sys.exit(1 if fails else 0)
