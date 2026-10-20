#!/usr/bin/env python3
"""Separate Ways v6, Appendix B: the rate at which graviton exchange in the bath makes N1.

Each Standard Model pair annihilates through an s-channel graviton into a massless fermion pair.
With L = -(1/M_P) h T (M_P the reduced Planck mass) and the de Donder propagator, the amplitude is
M = T_i P T_f / (M_P^2 s). This script builds T_i and T_f from explicit helicity spinors and
polarisation vectors, sums |M|^2 over spins, integrates over angle, and folds in Bose-Einstein and
Fermi-Dirac distributions. The coefficient alpha in dn_N/dt = alpha T^8/M_P^4 for one Majorana N
equals the Dirac pair-event rate (a Dirac field is two degenerate Majoranas).

Expected (by hand, Appendix B): angle-integrated spin sums pi/15 (real scalar pair), pi/5 (Dirac
pair), 4 pi/5 (one gauge boson, both helicities), in units s^2/M_P^4; alpha = 1.2027e-2.
Plants: the scalar's non-minimal coupling cannot matter (T_f is traceless); a polarisation
replaced by its momentum gives zero for the field-strength vertex and not for the vertex as printed
in Boyle, Turok, Vaibhav (arXiv:2609.11750, App. D), so their Feynman-gauge polarisation sum admits
unphysical states. Their printed spin-1 result (29/12 of ours, shape 15 + 5c^2 - 6c^4) does not
vanish along the beam, where angular momentum requires it to; it was not reproduced from their
printed vertex under either symmetrisation convention (scratchpad btv_variants.py), so only the
transversality failure and the forward value are claimed.
"""
import math, sys
import numpy as np
fails = 0
def check(ok, what):
    global fails
    print("   %-72s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)

ETA = np.diag([1.0, -1.0, -1.0, -1.0])
I2 = np.eye(2); Z2 = np.zeros((2, 2))
SIG = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
G = [np.block([[Z2, I2], [I2, Z2]]).astype(complex)] + [np.block([[Z2, s], [-s, Z2]]) for s in SIG]   # chiral basis

def lower(v): return ETA @ v
def dot(a, b): return a @ ETA @ b
def slash(p): return sum(G[m]*lower(p)[m] for m in range(4))
def bar(w): return w.conj() @ G[0]

def helicity_states(n):
    """Two-component spinors with sigma.n xi = +/- xi."""
    sn = sum(n[i]*SIG[i] for i in range(3))
    w, V = np.linalg.eigh(sn)
    return {+1: V[:, np.argmax(w)], -1: V[:, np.argmin(w)]}

def u_spinor(p, h):
    E = p[0]; xi = helicity_states(p[1:]/E)[h]
    return math.sqrt(2*E)*(np.concatenate([np.zeros(2), xi]) if h > 0 else np.concatenate([xi, np.zeros(2)]))

def v_spinor(p, h):
    E = p[0]; eta = helicity_states(p[1:]/E)[-h]
    return math.sqrt(2*E)*(np.concatenate([eta, np.zeros(2)]) if h > 0 else np.concatenate([np.zeros(2), -eta]))

def T_fermion_annihilate(p1, s1, p2, s2):
    """<0|T^{mu nu}|f(p1) fbar(p2)> = (1/4) vbar(p2) [g^mu k^nu + g^nu k^mu] u(p1), k = p1 - p2."""
    vb, u = bar(v_spinor(p2, s2)), u_spinor(p1, s1); k = p1 - p2
    J = np.array([vb @ G[m] @ u for m in range(4)])
    return 0.25*(np.outer(J, k) + np.outer(k, J))

def T_fermion_create(p3, s3, p4, s4):
    """<chi(p3) chibar(p4)|T^{rho sigma}|0> = (1/4) ubar(p3) [...] v(p4), k = p3 - p4."""
    ub, v = bar(u_spinor(p3, s3)), v_spinor(p4, s4); k = p3 - p4
    J = np.array([ub @ G[m] @ v for m in range(4)])
    return 0.25*(np.outer(J, k) + np.outer(k, J))

def T_scalar(p1, p2, xi=0.0):
    """Two real massless scalars, non-minimal coupling xi (improvement term ~ xi (P^mu P^nu - eta P^2))."""
    P = p1 + p2
    return (np.outer(p1, p2) + np.outer(p2, p1) - ETA*dot(p1, p2)) - 2*xi*(np.outer(P, P) - ETA*dot(P, P))

def fs(p, e):
    """Field strength of one photon, f^{mu alpha} = p^mu e^alpha - p^alpha e^mu."""
    return np.outer(p, e) - np.outer(e, p)

def T_vector(p1, e1, p2, e2):
    f1, f2 = fs(p1, e1), fs(p2, e2)
    A = f1 @ ETA @ f2.T          # f1^{mu a} eta_{ab} f2^{nu b}
    B = f2 @ ETA @ f1.T
    trace = np.einsum("ma,nb,mn,ab->", f1, f2, ETA, ETA)
    return A + B - 0.5*np.linalg.inv(ETA)*trace

def T_vector_btv(p, xi_, q, xit):
    """BTV App. D, contracted vertex, with a_(g b_d) = (a_g b_d + a_d b_g)/2; returned with upper indices."""
    sym = lambda a, b: 0.5*(np.outer(a, b) + np.outer(b, a))
    pq, xx, xq, xtp = dot(p, q), dot(xi_, xit), dot(xi_, q), dot(xit, p)
    Vlow = 0.5*(xx*(ETA*pq - sym(lower(p), lower(q))) - xq*xtp*ETA - 2*pq*sym(lower(xi_), lower(xit))
                + xq*sym(lower(xit), lower(p)) + xtp*sym(lower(xi_), lower(q)))
    return ETA @ Vlow @ ETA

def T_vector_btv_plain(p, xi_, q, xit):
    """the same with a_(g b_d) read as the plain sum a_g b_d + a_d b_g, under which (D3) differs from the Maxwell
    tensor (times their overall 1/2) only in the third term, -2 p.q xi_(g xit_d) against -p.q"""
    sym = lambda a, b: np.outer(a, b) + np.outer(b, a)
    pq, xx, xq, xtp = dot(p, q), dot(xi_, xit), dot(xi_, q), dot(xit, p)
    Vlow = 0.5*(xx*(ETA*pq - sym(lower(p), lower(q))) - xq*xtp*ETA - 2*pq*sym(lower(xi_), lower(xit))
                + xq*sym(lower(xit), lower(p)) + xtp*sym(lower(xi_), lower(q)))
    return ETA @ Vlow @ ETA

def T_vector_maxwell_plain(p, xi_, q, xit):
    """the Maxwell tensor written as (D3) is, plain-sum brackets, with -p.q in the third term and no overall 1/2"""
    sym = lambda a, b: np.outer(a, b) + np.outer(b, a)
    pq, xx, xq, xtp = dot(p, q), dot(xi_, xit), dot(xi_, q), dot(xit, p)
    Vlow = (xx*(ETA*pq - sym(lower(p), lower(q))) - xq*xtp*ETA - pq*sym(lower(xi_), lower(xit))
            + xq*sym(lower(xit), lower(p)) + xtp*sym(lower(xi_), lower(q)))
    return ETA @ Vlow @ ETA

PROP = 0.5*(np.einsum('mr,ns->mnrs', ETA, ETA) + np.einsum('ms,nr->mnrs', ETA, ETA) - np.einsum('mn,rs->mnrs', ETA, ETA))
def amp(Ti, Tf, s): return np.einsum('mn,mnrs,rs->', Ti, PROP, Tf)/s     # in units of 1/M_P^2

def kin(s, c, phi=0.3):
    E = math.sqrt(s)/2; st = math.sqrt(1 - c*c)
    p1 = E*np.array([1, 0, 0, 1.0]); p2 = E*np.array([1, 0, 0, -1.0])
    n = np.array([st*math.cos(phi), st*math.sin(phi), c])
    return p1, p2, E*np.concatenate([[1], n]), E*np.concatenate([[1], -n])

def transverse_pols(p):
    n = p[1:]/p[0]; a = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(n, a); e1 /= np.linalg.norm(e1); e2 = np.cross(n, e1)
    return [np.concatenate([[0], e1]), np.concatenate([[0], e2])]

def finals(p3, p4):
    return [T_fermion_create(p3, a, p4, b) for a in (1, -1) for b in (1, -1)]

def spin_sum(channel, s, c, xi=0.0):
    p1, p2, p3, p4 = kin(s, c)
    Tfs = finals(p3, p4)
    if channel == "scalar":
        Tis = [T_scalar(p1, p2, xi)]
    elif channel == "fermion":
        Tis = [T_fermion_annihilate(p1, a, p2, b) for a in (1, -1) for b in (1, -1)]
    elif channel == "vector":
        Tis = [T_vector(p1, e1, p2, e2) for e1 in transverse_pols(p1) for e2 in transverse_pols(p2)]
    return sum(abs(amp(Ti, Tf, s))**2 for Ti in Tis for Tf in Tfs)

def vector_feynman_sum(s, c, vertex):
    """Sum over polarisations with sum e_mu e_nu -> -eta_mu nu for both photons (all four basis vectors)."""
    p1, p2, p3, p4 = kin(s, c); Tfs = finals(p3, p4)
    basis = [np.eye(4)[k] for k in range(4)]
    tot = 0.0
    for Tf in Tfs:
        A = np.zeros((4, 4), complex)       # M = e1_a e2_b A^{ab} with lower-index polarisations
        for a in range(4):
            for b in range(4):
                e1 = ETA @ basis[a]; e2 = ETA @ basis[b]          # upper-index vectors whose lower components are unit
                A[a, b] = amp(vertex(p1, e1, p2, e2), Tf, s)
        tot += np.einsum('ab,cd,ac,bd->', A, A.conj(), ETA, ETA).real
    return tot

if __name__ == "__main__":
    xs, ws = np.polynomial.legendre.leggauss(24)
    s = 1.0
    def ang(fn): return 2*math.pi*sum(w*fn(x) for x, w in zip(xs, ws))
    expect = {"scalar": math.pi/15, "fermion": math.pi/5, "vector": 4*math.pi/5}
    print("=== angle-integrated spin sums, units of s^2/M_P^4 ===")
    got = {}
    for ch in ("scalar", "fermion", "vector"):
        got[ch] = ang(lambda c: spin_sum(ch, s, c))
        print("   %-8s %.6f   (by hand %.6f)" % (ch, got[ch], expect[ch]))
        check(abs(got[ch]/expect[ch] - 1) < 1e-9, "%s channel matches the hand calculation" % ch)
    # angular shapes, against the helicity selection rules
    for ch, shape in (("scalar", lambda c: c*c*(1 - c*c)/8), ("fermion", lambda c: (4*c**4 - 3*c*c + 1)/16),
                      ("vector", lambda c: (1 - c**4)/4)):
        ok = all(abs(spin_sum(ch, s, c) - shape(c)) < 1e-12 for c in (-0.9, -0.3, 0.2, 0.7))
        check(ok, "%s: |M|^2 has the spin-2 shape" % ch)
    check(spin_sum("vector", s, 0.999999) < 1e-5, "vector -> fermion pair vanishes in the forward direction")
    # plants
    check(abs(ang(lambda c: spin_sum("scalar", s, c, xi=1/6))/got["scalar"] - 1) < 1e-12,
          "plant: the scalar's coupling to curvature does not change the rate")
    p1, p2, p3, p4 = kin(s, 0.4); e2 = transverse_pols(p2)[0]; Tf = finals(p3, p4)[1]
    check(abs(amp(T_vector(p1, p1, p2, e2), Tf, s)) < 1e-12, "plant: field-strength vertex vanishes for e1 -> p1")
    check(abs(amp(T_vector_btv(p1, p1, p2, e2), Tf, s)) > 1e-3, "plant: the BTV vertex does not vanish for e1 -> p1")
    r_ok = vector_feynman_sum(s, 0.4, T_vector)/spin_sum("vector", s, 0.4)
    check(abs(r_ok - 1) < 1e-10, "plant: for the transverse vertex the -eta sum equals the physical sum")
    btv = ang(lambda c: vector_feynman_sum(s, c, T_vector_btv))
    print("   BTV vertex as printed, Feynman-gauge sum: %.4f x physical, forward value %.3g (their printed result:"
          % (btv/got["vector"], vector_feynman_sum(s, 1.0, T_vector_btv)))
    print("   29/12 x physical, shape 15 + 5c^2 - 6c^4; not reproduced from the printed vertex, so not traced further)")
    check(vector_feynman_sum(s, 1.0, T_vector_btv) > 1e-3,
          "the printed BTV vertex gives a non-zero rate along the beam")
    def vector_physical_sum(c, vertex):
        p1, p2, p3, p4 = kin(s, c)
        return sum(abs(amp(vertex(p1, e1, p2, e2), Tf, s))**2
                   for e1 in transverse_pols(p1) for e2 in transverse_pols(p2) for Tf in finals(p3, p4))
    # their vertices are half the canonical ones and their propagator four times ours (their scalar cross-section is
    # exactly canonical), so their photon vertex enters |M|^2 with a factor 4 against ours
    check(abs(ang(lambda c: vector_physical_sum(c, T_vector_maxwell_plain))/got["vector"] - 1) < 1e-12,
          "plant: the Maxwell tensor written as (D3) is, plain-sum brackets, gives our channel")
    four = {}
    for name, V in (("(ab+ba)/2", T_vector_btv), ("ab+ba", T_vector_btv_plain)):
        four[name] = (4*ang(lambda c: vector_physical_sum(c, V))/got["vector"],
                      4*ang(lambda c: vector_feynman_sum(s, c, V))/got["vector"])
        print("   BTV (D3) with brackets read %-9s, over ours: physical polarisations %.3f, Feynman-gauge sum %.3f"
              % (name, four[name][0], four[name][1]))
    check(all(abs(v - 29/12) > 0.1 and abs(v - 1) > 0.03 for pair in four.values() for v in pair),
          "neither reading of (D3), with either polarisation sum, gives their (D4) or our channel")
    check(abs(ang(lambda c: vector_physical_sum(c, T_vector))/got["vector"] - 1) < 1e-12,
          "plant: the field-strength vertex with physical polarisations gives our channel")
    # their (D4), (f2^4/512)(15 + 5c^2 - 6c^4) averaged over polarisations, is (1/128)(...) in Einstein units by their
    # footnote 8 and (1/32)(...) summed over the two photons' polarisations
    d4 = ang(lambda c: (15 + 5*c*c - 6*c**4)/32)
    print("   their (D4) summed over polarisations, (15 + 5c^2 - 6c^4)/32 in our units, integrates to %.4f = 29 pi/15;"
          " ours to 4 pi/5" % d4)
    check(abs(d4/(29*math.pi/15) - 1) < 1e-12 and abs(d4/got["vector"] - 29/12) < 1e-9,
          "their (D4) is 29/12 of our gauge-boson channel")

    # thermal rate: R = (1/1024 pi^6) Int f1 f2 E1 E2 dE1 dE2 dcos12 x K s^2, s = 2 E1 E2 (1 - cos12)
    JB, JF = math.pi**4/15, 7*math.pi**4/120          # Int E^3 f dE / T^4
    pref = (32/3)/(1024*math.pi**6)
    per = {"scalar": 0.5*got["scalar"]*pref*JB**2, "fermion": got["fermion"]*pref*JF**2, "vector": 0.5*got["vector"]*pref*JB**2}
    alpha = 4*per["scalar"] + 22.5*per["fermion"] + 12*per["vector"]
    exact = math.pi**3*(1/162000 + 49/307200 + 1/4500)
    print("\n=== thermal coefficient, dn_N/dt = alpha T^8/M_P^4 (Majorana N; = Dirac pair-event rate) ===")
    print("   per real scalar %.3e, per Dirac fermion %.3e, per gauge boson %.3e" % (per["scalar"], per["fermion"], per["vector"]))
    print("   alpha = 4 x + 22.5 x + 12 x = %.4e   (closed form %.4e)" % (alpha, exact))
    check(abs(alpha/exact - 1) < 1e-9, "alpha matches pi^3 (1/162000 + 49/307200 + 1/4500)")
    # direct numerical check of the thermal folding
    from scipy.integrate import quad
    def fold(stat):
        f = (lambda E: 1/math.expm1(E)) if stat == "B" else (lambda E: 1/(math.exp(E) + 1))
        return quad(lambda E: E**3*f(E), 0, 80)[0]
    check(abs(fold("B")/JB - 1) < 1e-8 and abs(fold("F")/JF - 1) < 1e-8, "plant: Bose and Fermi moments by quadrature")
    clery = 4*per["scalar"] + 45*per["fermion"] + 12*per["vector"]
    print("   with 45 Dirac fermions (Clery et al. arXiv v2 Eq. 90; half their published Eq. A12): %.4e, %.2f times ours" % (clery, clery/alpha))
    # BTV, Einstein limit: dn/dt = 2 sum_i n_i^2 <sigma v>_i, n_i = xi_i zeta3 T^3/pi^2, <sigma v> = 24 A T^2/M^4
    z3 = 1.2020569
    def btv_rate(A0, Ah, A1, pair=2.0, scalar_half=1.0, weyl_states=94):
        S = (4*scalar_half*A0 + weyl_states*(0.75**2)*Ah + 24*A1)/(3840*math.pi)
        return pair*24*(z3/math.pi**2)**2*S
    b0 = btv_rate(4, 3, 29); b1 = btv_rate(4, 3, 12); b2 = btv_rate(4, 3, 12, pair=1.0); b3 = btv_rate(4, 3, 12, pair=1.0, scalar_half=0.5)
    print("   BTV Einstein limit as printed: %.4e (%.2f x ours)" % (b0, b0/alpha))
    # their averaged rate Gbar = C_IR T^5/M_P^4 with dn/dt = 2 n_tot Gbar, n_tot = 98.5 zeta3 T^3/pi^2
    C_IR = b0/(2*98.5*z3/math.pi**2)
    print("   the same as their C_IR = %.3e; the closed form they print, 1393/(25216 pi^3) = %.3e, is %.3f of it (1/zeta3 = %.3f)"
          % (C_IR, 1393/(25216*math.pi**3), 1393/(25216*math.pi**3)/C_IR, 1/z3))
    print("     transverse gauge vertex: %.4e; Majorana pair counting: %.4e; identical real scalars: %.4e (%.3f x ours)"
          % (b1, b2, b3, b3/alpha))
    check(1.0 < b3/alpha < 1.15, "corrected BTV rate within 15% of ours")
    b45 = btv_rate(4, 3, 12, pair=1.0, scalar_half=0.5, weyl_states=90)
    print("     with 45 Weyl fermions in place of their 47: %.4e (%.3f x ours); the rest is their statistics" % (b45, b45/alpha))
    print("   their C_IR needs the 29 of (D4): with our gauge-boson channel it would be %.2e, against the %.2e they quote"
          % (C_IR*b1/b0, C_IR))
    print("   their remaining %.3f = %.3f (1/2 for identical scalars) x %.3f (47 Weyl fermions) x %.3f (their statistics)"
          % (b2/alpha, b2/b3, b3/b45, b45/alpha))
    print("   of their 4.27 the errors make %.2f (vertex x Majorana count x identical scalars); the rest is the bath and statistics"
          % ((b0/b1)*(b1/b2)*(b2/b3)))
    check(abs((b2/b3)*(b3/b45)*(b45/alpha)/(b2/alpha) - 1) < 1e-12, "the three factors multiply to the remainder")
    print("   their 47 Weyl fermions against 45 raise the rate by %.1f per cent in their statistics, %.1f in ours"
          % (100*(b3/b45 - 1), 100*(((4*per["scalar"] + 23.5*per["fermion"] + 12*per["vector"])/alpha) - 1)))
    # the remainder: BTV put quantum number densities into a Maxwell-Boltzmann <sigma v>. Done consistently in
    # Maxwell-Boltzmann statistics (n = T^3/pi^2 per state, Int E^3 e^-E dE = 6) the two calculations agree.
    ours_mb = 4*0.5*got["scalar"]*pref*36 + 22.5*got["fermion"]*pref*36 + 12*0.5*got["vector"]*pref*36
    btv_mb = 24/math.pi**4*(4*4*0.5 + 94*3 + 24*12)/(3840*math.pi)
    print("   Maxwell-Boltzmann throughout: ours %.4e, corrected BTV %.4e (ratio %.3f)" % (ours_mb, btv_mb, btv_mb/ours_mb))
    check(abs(btv_mb/ours_mb - 1) < 0.03, "in consistent Boltzmann statistics the corrected BTV rate matches ours to 3%")
    check(abs(ours_mb - 1.156e-2) < 2e-5, "Boltzmann value equals the 1.16e-2 found from Clery's amplitudes (lit A2)")
    # BTV's bath has 47 Weyl fermions (N2 and N3 thermalised), ours 45; with 45 the two agree exactly
    btv_mb45 = 24/math.pi**4*(4*4*0.5 + 90*3 + 24*12)/(3840*math.pi)
    print("   corrected BTV with 45 Weyl fermions instead of 47: %.6e against ours %.6e" % (btv_mb45, ours_mb))
    check(abs(btv_mb45/ours_mb - 1) < 1e-9, "the 2 per cent Boltzmann residual is BTV's 47 Weyl fermions against 45")
    # Above T* their quadratic-gravity rate has the same channels with sigma ~ f2^4/s: their C_UV is exactly C_IR/768,
    # which is what quarter-size channels and Maxwell-Boltzmann moments give. A rate there weighs the first moment
    # Int E f dE of each species where the Einstein limit weighs the third, so the same corrections move their two
    # coefficients by different factors, and their yield, Eq. (11), goes as C_IR^(1/4) C_UV^(3/4) f2^3.
    check(abs(1393/(19365888*math.pi**3)*768/(1393/(25216*math.pi**3)) - 1) < 1e-12, "their C_UV is exactly C_IR/768")
    def btv_over_ours(k):
        """their rate over ours for sigma ~ s^k: their channels (no 1/2 for scalars, 29/12 for gauge bosons, 47 Weyl
        fermions, two N1 per event) with quantum densities in a Maxwell-Boltzmann average, against ours"""
        JB_, JF_, JM = {1: (math.pi**4/15, 7*math.pi**4/120, 6.0), -1: (math.pi**2/6, math.pi**2/12, 1.0)}[k]
        ours_k = 4*0.5*got["scalar"]*JB_**2 + 22.5*got["fermion"]*JF_**2 + 12*0.5*got["vector"]*JB_**2
        theirs_k = 2*(4*got["scalar"]*(z3*JM)**2 + 23.5*got["fermion"]*(0.75*z3*JM)**2
                      + 12*0.5*(29/12)*got["vector"]*(z3*JM)**2)
        return theirs_k/ours_k
    def btv_steps(k):
        """their coefficient over ours after each correction in turn: as printed, transverse vertex, one N1 per
        event, 1/2 for identical scalars, then our statistics (the last step makes it 1 by construction)"""
        JB_, JF_, JM = {1: (math.pi**4/15, 7*math.pi**4/120, 6.0), -1: (math.pi**2/6, math.pi**2/12, 1.0)}[k]
        ours_k = 4*0.5*got["scalar"]*JB_**2 + 22.5*got["fermion"]*JF_**2 + 12*0.5*got["vector"]*JB_**2
        def theirs(vec, pair, half, ndirac=23.5):
            return pair*(4*half*got["scalar"]*(z3*JM)**2 + ndirac*got["fermion"]*(0.75*z3*JM)**2
                         + 12*0.5*vec*got["vector"]*(z3*JM)**2)/ours_k
        return [theirs(29/12, 2, 1), theirs(1, 2, 1), theirs(1, 1, 1), theirs(1, 1, 0.5), theirs(1, 1, 0.5, 22.5), 1.0]
    st_ir, st_uv = btv_steps(1), btv_steps(-1)
    print("   their Einstein-limit coefficient over ours, step by step: " + ", ".join("%.3f" % x for x in st_ir))
    print("   their coefficient above T* over ours, step by step:      " + ", ".join("%.3f" % x for x in st_uv))
    check(abs(st_ir[0]/(b0/alpha) - 1) < 1e-9 and abs(st_ir[1]/(b1/alpha) - 1) < 1e-9 and abs(st_ir[3]/(b3/alpha) - 1) < 1e-9,
          "the step model reproduces Table B.1's Einstein-limit steps")
    r_ir, r_uv = btv_over_ours(1), btv_over_ours(-1)
    print("   their coefficients over ours: Einstein limit %.3f, above T* %.3f" % (r_ir, r_uv))
    check(abs(r_ir/(b0/alpha) - 1) < 1e-9, "the same moments reproduce their Einstein-limit ratio")
    f2_bound = 7e-5*r_ir**(1/12)*r_uv**(1/4)
    print("   their yield goes as C_IR^(1/4) C_UV^(3/4) f2^3, so their bound f2 < 7e-5 rises to %.2e" % f2_bound)
    print("   by the factor (%.2f^(1/4) %.2f^(3/4))^(1/3) = %.2f" % (r_ir, r_uv, (r_ir**0.25*r_uv**0.75)**(1/3)))
    # Thermal masses. A mass m in the colliding pair changes the angle-integrated rate at relative order m^2 <1/s>,
    # with <1/s> the mean of 1/s over the rate: Int f f E1 E2 s / Int f f E1 E2 s^2 = (3/8)(I2/I3)^2 / T^2, where
    # I_n = Int E^n f dE and s = 2 E1 E2 (1 - cos). The largest is the gluon's asymptotic mass, m^2 = g3^2 T^2,
    # with alpha_s run at one loop with six flavours from 0.1180 at M_Z.
    I2B, I3B, I2F, I3F = 2*z3, math.pi**4/15, 1.5*z3, 7*math.pi**4/120
    invs_b, invs_f = 0.375*(I2B/I3B)**2, 0.375*(I2F/I3F)**2
    def invs_numeric(boson):
        sgn = -1 if boson else 1
        num = quad(lambda E: E**2/(math.exp(E) + sgn), 0, 60)[0]**2*2*2
        den = quad(lambda E: E**3/(math.exp(E) + sgn), 0, 60)[0]**2*4*(8/3)
        return num/den
    check(abs(invs_numeric(True)/invs_b - 1) < 1e-8 and abs(invs_numeric(False)/invs_f - 1) < 1e-8,
          "the closed forms for <1/s> match direct integration")
    # Einstein-Cartan with non-minimal fermion couplings alpha, beta: Shaposhnikov et al.'s Eq. (6) scales the
    # Majorana rate by C_M = (9/4)(24 a^4 + 21 b+^4), a^2 = 1 + alpha^2 - beta^2, b+^2 = 1 - (alpha + beta)^2
    def C_M(al, be):
        return 2.25*(24*(1 + al**2 - be**2)**2 + 21*(1 - (al + be)**2)**2)
    print("   Einstein-Cartan, non-minimal: C_M(0,0) = %.2f; C_M(0,+1) = %.1f, C_M(0,-1) = %.1f; C_M(10,0)/C_M(0,0) = %.0f"
          % (C_M(0, 0), C_M(0, 1), C_M(0, -1), C_M(10, 0)/C_M(0, 0)))
    check(C_M(0, 1) == 0 and C_M(0, -1) == 0 and C_M(10, 0) > 1e3*C_M(0, 0),
          "the contact term vanishes at alpha = 0, beta = +-1 and grows without limit")
    g3sq = 4*math.pi/(1/0.1180 + 7/(2*math.pi)*math.log(1e14/91.1876))
    print("   thermal masses: <1/s> = %.4f/T^2 for bosons, %.4f/T^2 for fermions; at 1e14 GeV g3^2 = %.3f, so the gluon's"
          " m^2 <1/s> = %.3f" % (invs_b, invs_f, g3sq, g3sq*invs_b))
    # every species at 1e14 GeV, one-loop couplings from M_Z (alpha_2^-1 = 29.6, GUT-normalised alpha_1^-1 = 59.0,
    # b = -19/6, 41/10) and y_t = 0.5: asymptotic masses g3^2 T^2 (gluons), (11/12) g^2 T^2 (W, B), g3^2 T^2/3 (quarks),
    # (3 g2^2 + g'^2)/32 T^2 (leptons), (3 g2^2/16 + g'^2/16 + y_t^2/4) T^2 (Higgs); each weighted by its share of alpha
    L = math.log(1e14/91.1876)
    g2sq = 4*math.pi/(29.6 + 19/6/(2*math.pi)*L); g1sq = 4*math.pi/(59.0 - 4.1/(2*math.pi)*L)*3/5
    share_v, share_f, share_s = 12*per["vector"]/alpha, 22.5*per["fermion"]/alpha, 4*per["scalar"]/alpha
    shift = (share_v*(8*g3sq + 3*(11/12)*g2sq + (11/12)*g1sq)/12*invs_b
             + share_f*(36*g3sq/3 + 9*(3*g2sq + g1sq)/32)/45*invs_f
             + share_s*(3*g2sq/16 + g1sq/16 + 0.25/4)*invs_b)
    print("   all species at 1e14 GeV (g2^2 = %.3f, g'^2 = %.3f): the rate shifts at relative order %.4f" % (g2sq, g1sq, shift))
    check(0 < shift < g3sq*invs_b, "the shift of the whole rate is below the gluon channel's")
    # below the bound the plasma's yield goes as f2^3, so its share of the abundance is (f2/bound)^3
    # the bath's two extra Weyl fermions are a choice of model, not an error: the bound with their bath kept
    w47_ir = b3/b45
    print("   keeping their 47-fermion bath moves that bound by %.2f per cent" % (100*(1 - (1/w47_ir)**(1/12))))
    print("   below that bound the plasma makes under one per cent of the dark matter for f2 < %.3g" % (0.01**(1/3)*f2_bound))
    check(f2_bound < 7e-5*(b0/alpha)**(1/3), "plant: a single factor 4.27 on both coefficients would move the bound further")
    # Free evolution through the crossing: the rate at which graviton exchange changes N1 number, per N1, at the
    # equilibrium density, Gamma = R/n_eq, over H. In Einstein gravity it grows as T^3; in quadratic gravity it peaks
    # at T*, where the two rates meet, T* = (C_UV/C_IR)^(1/4) f2 M_P with the corrected coefficients.
    MP = 2.435e18                                      # GeV, the reduced Planck mass
    def gamma_over_H(T):
        neq = 0.75*z3/math.pi**2*2*T**3
        return alpha*T**8/MP**4/neq/(math.sqrt(math.pi**2*106.75/90)*T*T/MP)
    print("   Einstein gravity: Gamma/H = %.1e at the reduced Planck temperature, 1 at %.1f M_P" % (gamma_over_H(MP), gamma_over_H(MP)**(-1/3)))
    T_star_q = (1/768*r_ir/r_uv)**0.25*f2_bound*MP
    print("   quadratic gravity at f2 = %.1e: T* = %.3e GeV, where Gamma/H = %.1e, its largest value" % (f2_bound, T_star_q, gamma_over_H(T_star_q)))
    check(gamma_over_H(T_star_q) < 1e-10, "in quadratic gravity N1 never interacts while it crosses (Gamma/H far below 1)")
    # integrated across the crossing, from the bang to H = M (T_c = 1.87e13 GeV): int Gamma dt = int (Gamma/H) dlnT, with
    # Gamma/H growing as T^3 below T* and falling as 1/T above it
    T_c, g_pk = 1.87e13, gamma_over_H(T_star_q)
    gq = lambda lnT: g_pk*(math.exp(lnT)/T_star_q)**3 if math.exp(lnT) < T_star_q else g_pk*T_star_q/math.exp(lnT)
    n_int = quad(gq, math.log(T_c), math.log(T_star_q))[0] + quad(gq, math.log(T_star_q), math.log(T_star_q) + 60)[0]
    print("   across the crossing, from the bang to H = M, each N1 meets %.1e graviton-exchange interactions" % n_int)
    check(abs(n_int/(g_pk*((1 - (T_c/T_star_q)**3)/3 + 1)) - 1) < 1e-6, "the integral matches its closed form")
    mD2 = 2*g3sq                                      # gluon Debye mass squared over T^2, the largest
    # the rate-weighted mean of s for bosons: Int f f E1 E2 s^3 / Int f f E1 E2 s^2 = 2 (I4/I3)^2 (3/2) T^2
    s_mean = 2*(24*1.0369277551/(math.pi**4/15))**2*1.5
    print("   elastic scattering, forward-peaked and cut at the Debye mass, is faster by about <s>/m_D^2 = %.0f/%.2f = %.0f,"
          " so %.0e per N1 across the crossing" % (s_mean, mD2, s_mean/mD2, s_mean/mD2*n_int))
    kappa = math.sqrt(math.pi*4.916e8*math.sqrt(106.75/15)/5.965e18)
    print("   an N1 the plasma makes with momentum 3T has x = %.1e; the crossing fills x of order 1" % (3/kappa))
    print("   Clery et al.'s 45-Dirac value is 11351 pi^3/20736000 = %.6e; ours with 45 Dirac %.6e" % (11351*math.pi**3/20736000, clery))
    check(abs(clery/(11351*math.pi**3/20736000) - 1) < 1e-9, "with 45 Dirac fermions the rate equals Clery et al.'s closed form")
    # Clery et al.'s published PRD numbers this rate Eq. (A12) and counts the two particles each scattering
    # makes, so it is 11351 pi^3/10368000, twice their arXiv v2 Eq. (90) and twice our 45-Dirac event sum
    print("   Clery et al.'s published Eq. (A12) is 11351 pi^3/10368000 = %.6e, %.3f x our 45-Dirac event sum"
          % (11351*math.pi**3/10368000, 11351*math.pi**3/10368000/clery))
    check(abs(11351*math.pi**3/10368000/clery - 2) < 1e-9, "Clery et al.'s published fermion rate is twice the 45-Dirac event rate")
    # Bernal et al. 2018 (PRD, checked at source): vertex h T/2M_P, Eq. (3), with the canonical propagator, Eq. (A1),
    # so every squared amplitude is 1/16 of the canonical one; their rate is Eq. (B17)
    bernal = 11351*math.pi**3/331776000
    print("   Bernal et al.'s Eq. (B17) is 11351 pi^3/331776000 = %.4e, %.4f x ours, %.4f of our 45-Dirac sum"
          % (bernal, bernal/alpha, bernal/clery))
    check(abs(bernal/clery - 1/16) < 1e-12, "Bernal et al.'s rate is exactly 1/16 of the 45-Dirac sum (a vertex half the canonical one)")
    print("   Boyle, Turok and Vaibhav per event (one pair per collision): %.3f x ours" % (b0/2/alpha))
    # Tang and Wu 2017 tabulate A = kappa^-4 Int dcos sum|M|^2 with kappa^4 = 16/M_P^4 (their Eq. 3.7 and Table I)
    tw = [got[ch]/(2*math.pi)/16 for ch in ("scalar", "fermion", "vector")]
    print("   in Tang and Wu's normalisation our channels are A = s^2/%.0f, s^2/%.0f, s^2/%.0f (their Table I: 480, 160, 40)"
          % tuple(1/x for x in tw))
    check(all(abs(x*d - 1) < 1e-9 for x, d in zip(tw, (480, 160, 40))), "our channels equal Tang and Wu's Table I")
    # Garny et al. 2018 (JCAP, checked at source): their Eq. (3.1) weights the Dirac all-spin sum by 45; their real-scalar
    # and vector stress tensors, Eqs. (3.5) and (3.9), carry a factor 1/2, so those squared amplitudes are a quarter of
    # ours; and their cross-section, Eq. (3.3), integrates over all angles with no 1/2 for an identical pair, which
    # doubles those two channels again. In their rate the channels therefore stand at 1/2, 2 and 1/2 of ours.
    w_garny = {"scalar": (1/4)*2, "fermion": 45/22.5, "vector": (1/4)*2}
    print("   Garny et al.'s channels in their rate, relative to ours: scalar %.2f, fermion %.2f, vector %.2f"
          % (w_garny["scalar"], w_garny["fermion"], w_garny["vector"]))
    garny = 4*per["scalar"]*w_garny["scalar"] + 22.5*per["fermion"]*w_garny["fermion"] + 12*per["vector"]*w_garny["vector"]
    print("   with those weights and Bose-Einstein and Fermi-Dirac statistics: %.4e, %.3f x ours (they print no total)" % (garny, garny/alpha))
    # their own thermal average: production <sigma v>(n_eq)^2, their (2.5), n_eq = T^3/pi^2 (g = 1), and their (3.17) in
    # the massless limit, <sigma v>_0 = (2 pi/5) T^2/m_p^4 and <sigma v>_1/2 = <sigma v>_1 = (24 pi/5) T^2/m_p^4, m_p^4 = 64 pi^2 M_P^4
    garny_mb = (4*2*math.pi/5 + 45*24*math.pi/5 + 12*24*math.pi/5)/math.pi**4/(64*math.pi**2)
    print("   with their own Maxwell-Boltzmann average, (2.5) and (3.17): %.5e = 43/(10 pi^5), %.3f x ours" % (garny_mb, garny_mb/alpha))
    check(abs(garny_mb/(43/(10*math.pi**5)) - 1) < 1e-12, "Garny et al.'s formulas give 43/(10 pi^5) in their average")
    check(w_garny["scalar"] == 0.5 and w_garny["vector"] == 0.5 and w_garny["fermion"] == 2.0,
          "Garny et al.: the 1/2 in their stress tensors and the missing 1/2 for identical pairs leave half; the count doubles")

    # The published channels themselves, massless, in units of s^2/M_P^4 with G = 1/(8 pi M_P^2), so G^2 pi^2 = 1/64:
    # Garny et al. JCAP (3.14)-(3.16), Bernal et al. PRD (B11)-(B13), Clery et al. PRD (A9)-(A11), each as printed
    def published(c):
        t = -(1 - c)/2; u = -t*(1 + t)                    # s = 1, so -t(s+t) = (1 - c^2)/4 and s + 2t = c
        poly = 1 + 10*t + 42*t*t + 64*t**3 + 32*t**4
        return {"Garny":  ((8/64)*u*(1 + 2*t)**2, (8/64)*poly, (32/64)*u*(1 + 2*u*(-1))),
                "Bernal": (u*(1 + 2*t)**2/32, poly/128, u*(1 - 2*u)/8),
                "Clery":  (u*(1 + 2*t)**2/4, poly/8, u*(1 - 2*u))}
    ratios = {}
    for name in ("Garny", "Bernal", "Clery"):
        rs = []
        for c in (-0.7, -0.2, 0.3, 0.8):
            ours = (spin_sum("scalar", s, c), spin_sum("fermion", s, c), spin_sum("vector", s, c))
            rs.append(tuple(p/o for p, o in zip(published(c)[name], ours)))
        spread = max(abs(r[k] - rs[0][k]) for r in rs for k in range(3))
        ratios[name] = (rs[0], spread)
        print("   %s et al.'s channels over ours (scalar, fermion, vector): %.4f, %.4f, %.4f at every angle (spread %.0e)"
              % (name, rs[0][0], rs[0][1], rs[0][2], spread))
    check(all(v[1] < 1e-12 for v in ratios.values()), "each published channel is a fixed multiple of ours at every angle")
    check(max(abs(a_ - b_) for a_, b_ in zip(ratios["Garny"][0], (0.25, 1, 0.25))) < 1e-12, "Garny et al.: 1/4, 1, 1/4")
    check(max(abs(a_ - 1/16) for a_ in ratios["Bernal"][0]) < 1e-12, "Bernal et al.: 1/16 in every channel")
    check(max(abs(a_ - b_) for a_, b_ in zip(ratios["Clery"][0], (0.5, 1, 0.5))) < 1e-12,
          "Clery et al.: ours, with the 1/2 for an identical pair in the scalar and vector channels")
    # the count in energy densities: per pair of single-particle states the rate is A rho_a rho_b/(24 pi^2 M_P^4), with
    # rho the energy density of one state; a Weyl field gives one coupled pair, A = pi/10, rho_f = (7/8)(pi^2/30) T^4
    rho_f = 7/8*math.pi**2/30
    fermion_rho = 45*(math.pi/10)*rho_f**2/(24*math.pi**2)
    print("   fermions from energy densities: 45 (pi/10) rho_f^2/(24 pi^2) = %.6e; from the Dirac count, 22.5 x %.6e = %.6e"
          % (fermion_rho, per["fermion"], 22.5*per["fermion"]))
    check(abs(fermion_rho/(22.5*per["fermion"]) - 1) < 1e-9, "45 coupled Weyl pairs of states give the 22.5-Dirac fermion channel")
    # A Weyl pair is one of the two helicity combinations of a Dirac pair, and carries half its sum
    p1, p2, p3, p4 = kin(s, 0.37); Tfs = finals(p3, p4)
    combos = {(a_, b_): ang(lambda c, a_=a_, b_=b_: sum(abs(amp(T_fermion_annihilate(*kin(s, c)[:1], a_, kin(s, c)[1], b_), Tf_c, s))**2
                                                    for Tf_c in finals(*kin(s, c)[2:])))
              for a_ in (1, -1) for b_ in (1, -1)}
    live = sorted(v for v in combos.values() if v > 1e-12)
    print("   fermion channel by helicity combination: %s" % ", ".join("%.6f" % v for v in sorted(combos.values())))
    check(len(live) == 2 and abs(live[0]/live[1] - 1) < 1e-10 and abs(live[0]/got["fermion"] - 0.5) < 1e-10,
          "a Weyl pair (one helicity combination) carries exactly half the Dirac sum")
    # with the factor 1/2 for an identical pair the channels are the widths of a massive spin-2 particle into the same pairs
    w = (got["scalar"]/2, got["fermion"], got["vector"]/2)
    print("   with 1/2 for identical pairs the channels stand 1 : %.3f : %.3f (real scalar pair, Dirac pair, gauge-boson pair)"
          % (w[1]/w[0], w[2]/w[0]))
    check(abs(w[1]/w[0] - 6) < 1e-9 and abs(w[2]/w[0] - 12) < 1e-9, "the channels stand 1 : 6 : 12, the spin-2 widths of Han, Lykken and Zhang")
    check(abs(got["vector"]/got["scalar"] - 6) > 1, "plant: without the identical-pair 1/2 the scalar-to-vector ratio is not 6")

    # Einstein-Cartan gravity with minimally coupled fermions: eliminating torsion leaves
    # L = -(3/16 M_P^2) A.A, A the total axial current (Shaposhnikov et al., PRL 126 (2021) 161301, Eq. 2,
    # alpha = beta = 0). The cross term between a Standard Model Dirac current and N's is (3/8)/M_P^2;
    # with the Dirac-pair counting above it gives one Majorana N1 per pair event, as for the graviton.
    G5 = 1j*G[0] @ G[1] @ G[2] @ G[3]
    def contact_amps(c):
        p1, p2, p3, p4 = kin(s, c)
        out = []
        for a in (1, -1):
            for b in (1, -1):
                Ji = np.array([bar(v_spinor(p2, b)) @ G[m] @ G5 @ u_spinor(p1, a) for m in range(4)])
                Ti = T_fermion_annihilate(p1, a, p2, b)
                for c3 in (1, -1):
                    for c4 in (1, -1):
                        Jf = np.array([bar(u_spinor(p3, c3)) @ G[m] @ G5 @ v_spinor(p4, c4) for m in range(4)])
                        Mc = 0.375*dot(Ji, Jf)
                        Mg = amp(Ti, T_fermion_create(p3, c3, p4, c4), s)
                        out.append((Mc, Mg))
        return out
    con = lambda c: sum(abs(m[0])**2 for m in contact_amps(c))
    shape = lambda c: (9/8)*(((1 - c)/2)**2 + ((1 + c)/2)**2)          # (9/8)(t^2 + u^2), s = 1
    check(all(abs(con(c) - shape(c)) < 1e-12 for c in (-0.8, -0.1, 0.5, 0.95)), "contact |M|^2 = (9/8)(t^2 + u^2)/M_P^4 per Dirac pair")
    ec_ang = ang(con)
    print("\n=== Einstein-Cartan contact term, minimal coupling ===")
    print("   angle-integrated spin sum %.6f s^2/M_P^4 (3 pi = %.6f); %.3f times the graviton's Dirac channel" % (ec_ang, 3*math.pi, ec_ang/got["fermion"]))
    check(abs(ec_ang/(3*math.pi) - 1) < 1e-9, "contact channel integrates to 3 pi")
    inter = ang(lambda c: sum((m[0]*m[1].conjugate()).real for m in contact_amps(c)))
    inter_fixed = abs(sum((m[0]*m[1].conjugate()).real for m in contact_amps(0.5)))
    print("   interference with graviton exchange: %.1e after angular integration, %.3f at cos = 0.5" % (inter, inter_fixed))
    check(abs(inter) < 1e-10, "graviton (J = 2) and contact (J = 1) amplitudes do not interfere in the rate")
    check(inter_fixed > 1e-3, "plant: at fixed angle the interference is not zero, so the test can fail")
    alpha_ec = 22.5*per["fermion"]*ec_ang/got["fermion"]
    # T_star from v6_bath.py, 8.008e13 GeV, scales as the rate to the power -1/3
    print("   alpha_EC = %.4e = %.2f alpha; total %.4e = %.2f alpha; T_star falls by %.3f, from 8.008e13 to %.2e GeV"
          % (alpha_ec, alpha_ec/alpha, alpha + alpha_ec, (alpha + alpha_ec)/alpha, ((alpha + alpha_ec)/alpha)**(1/3),
             8.008e13/((alpha + alpha_ec)/alpha)**(1/3)))
    # Shaposhnikov et al.'s own normalisation: Omega_N/Omega_DM = 3.6e-2 C_M (M/10 keV)(T/M_P)^3, C_M = (9/4) 45 at
    # alpha = beta = 0, against ours Omega/Omega_DM = M (a/C)(T/M_P)^3 / (rho/s)
    a_sst = 3.6e-2*(9/4)*45/1e-5*(4.3664e-10)*480.6
    print("   from Shaposhnikov et al.'s abundance formula: %.3e (%.3f of ours; their r(y) is approximate)" % (a_sst, a_sst/alpha_ec))
    check(abs(a_sst/alpha_ec - 1) < 0.05, "agrees with Shaposhnikov et al.'s abundance formula to 5 per cent")

    # Coefficients for other dark matter. Event rate R = pi^3 (N_0 r_0 + N_1/2 r_1/2 + N_1 r_1) T^8/M_P^4, with
    # N_0 real scalars, N_1/2 Weyl fermions and N_1 gauge bosons in the bath, r_i per unit of each (a Weyl
    # fermion is half a Dirac spin sum; identical initial bosons and identical final particles carry 1/2).
    from fractions import Fraction
    def Tf_dm(dm, p3, p4):
        if dm == "real scalar, xi = 0": return [T_scalar(p3, p4, 0.0)]
        if dm == "real scalar, xi = 1/6": return [T_scalar(p3, p4, 1/6)]
        if dm == "Dirac fermion": return finals(p3, p4)
        return [T_vector(p3, e3, p4, e4) for e3 in transverse_pols(p3) for e4 in transverse_pols(p4)]
    def Ti_sm(ch, p1, p2):
        if ch == "scalar": return [T_scalar(p1, p2, 0.0)]
        if ch == "fermion": return [T_fermion_annihilate(p1, a, p2, b) for a in (1, -1) for b in (1, -1)]
        return [T_vector(p1, e1, p2, e2) for e1 in transverse_pols(p1) for e2 in transverse_pols(p2)]
    def channel(ch, dm):
        def f(c):
            p1, p2, p3, p4 = kin(s, c)
            return sum(abs(amp(Ti, Tf, s))**2 for Ti in Ti_sm(ch, p1, p2) for Tf in Tf_dm(dm, p3, p4))
        return ang(f)
    print("\n=== graviton-exchange coefficients for other dark matter (event rates, units pi^3 T^8/M_P^4) ===")
    print("   per real scalar, per Weyl fermion, per gauge boson of the bath; SM: 4, 45, 12")
    table = {}
    for dm, ident, ppe in (("real scalar, xi = 0", 0.5, 2), ("real scalar, xi = 1/6", 0.5, 2), ("Dirac fermion", 1.0, 2),
                           ("massless vector", 0.5, 2)):
        A = {ch: channel(ch, dm) for ch in ("scalar", "fermion", "vector")}
        r = {"scalar": ident*0.5*A["scalar"]*pref*JB**2/math.pi**3, "fermion": ident*0.5*A["fermion"]*pref*JF**2/math.pi**3,
             "vector": ident*0.5*A["vector"]*pref*JB**2/math.pi**3}
        R = math.pi**3*(4*r["scalar"] + 45*r["fermion"] + 12*r["vector"])
        R45 = math.pi**3*(4*r["scalar"] + 90*r["fermion"] + 12*r["vector"])
        table[dm] = (r, R, R45, A)
        fr = [str(Fraction(r[k]).limit_denominator(10**9)) for k in ("scalar", "fermion", "vector")]
        print("   %-22s r = %s, %s, %s; R = %.4e; particles per event %d; with 45 Dirac %.4e, %.3f times too high"
              % (dm, fr[0], fr[1], fr[2], R, ppe, R45, R45/R))
    mv = {k: table["massless vector"][0][k] + table["real scalar, xi = 0"][0][k] for k in ("scalar", "fermion", "vector")}
    Rmv = table["massless vector"][1] + table["real scalar, xi = 0"][1]; Rmv45 = table["massless vector"][2] + table["real scalar, xi = 0"][2]
    fr = [str(Fraction(mv[k]).limit_denominator(10**9)) for k in ("scalar", "fermion", "vector")]
    print("   %-22s r = %s, %s, %s; R = %.4e; particles per event 2; with 45 Dirac %.4e, %.3f times too high"
          % ("massive vector, T >> m", fr[0], fr[1], fr[2], Rmv, Rmv45, Rmv45/Rmv))
    Rmv6 = table["massless vector"][1] + table["real scalar, xi = 1/6"][1]
    print("   with a conformally coupled longitudinal mode the massive-vector row falls by %.1f per cent" % (100*(1 - Rmv6/Rmv)))
    rs = table["real scalar, xi = 0"]
    print("   Clery et al.'s scalar coefficient 3997 pi^3/41472000 = %.6e; ours with 45 Dirac %.6e" % (3997*math.pi**3/41472000, rs[2]))
    check(abs(rs[2]/(3997*math.pi**3/41472000) - 1) < 1e-9, "the 45-Dirac scalar rate equals Clery et al.'s arXiv v2 coefficient")
    print("   Clery et al.'s published Eq. (A8) is 3997 pi^3/20736000 = %.6e, %.3f x our 45-Dirac scalar event sum"
          % (3997*math.pi**3/20736000, 3997*math.pi**3/20736000/rs[2]))
    check(abs(3997*math.pi**3/20736000/rs[2] - 2) < 1e-9, "Clery et al.'s published scalar rate is twice the 45-Dirac event rate")
    check(abs(table["Dirac fermion"][1]/alpha - 1) < 1e-9, "the Dirac event rate equals alpha (one Majorana N1 per Dirac-pair event)")
    xi0, xi6 = table["real scalar, xi = 0"][3], table["real scalar, xi = 1/6"][3]
    print("   scalar pair from plasma scalars: xi = 0 integrates to %.6f (2 pi/15 = %.6f), %.3f times xi = 1/6"
          % (xi0["scalar"], 2*math.pi/15, xi0["scalar"]/xi6["scalar"]))
    check(abs(xi0["scalar"]/(2*math.pi/15) - 1) < 1e-9 and abs(xi0["scalar"]/xi6["scalar"] - 6) < 1e-9,
          "minimally coupled scalars: -tu/(M_P^2 s), integral 2 pi/15, six times the conformal value")
    check(abs(xi0["fermion"]/xi6["fermion"] - 1) < 1e-12 and abs(xi0["vector"]/xi6["vector"] - 1) < 1e-12 and abs(xi0["scalar"]/xi6["scalar"] - 1) > 0.5,
          "only the scalar channel depends on the dark matter's xi (fermion and gauge stress tensors are traceless)")
    print("   Majorana fermion: event rate %.4e (half the Dirac one), %.4e particles per unit volume and time" % (alpha/2, alpha))
    b_vertex, b_count, b_rest = b0/b1, b1/b2, b2/alpha
    print("   BTV's factor %.2f = vertex %.3f x Majorana count %.3f x scalar factor and statistics %.3f" % (b0/alpha, b_vertex, b_count, b_rest))
    print("   at fixed abundance a fermion's T_max rises by 1.41^(1/3) = %.3f" % (clery/alpha)**(1/3))
    a47 = 4*per["scalar"] + 23.5*per["fermion"] + 12*per["vector"]
    print("   with N2 and N3 in the bath (47 Weyl fermions) alpha rises by %.1f per cent" % (100*(a47/alpha - 1)))
    # Clery et al.'s published rates count two particles per event, so with the 45 count their thermal abundances are
    # high by the count factor alone. arXiv v2 entered the scalar event rate directly as dn_X/dt, without that factor
    # 2, which made its thermal scalar abundance low by 2/1.381 instead.
    print("   Clery et al.'s published thermal abundances: high by %.3f (fermion) and %.3f (scalar), the 45 count alone"
          % (clery/alpha, rs[2]/rs[1]))
    print("   (their arXiv v2 thermal scalar abundance: x %.3f from the 45 count, /2 from events entered as particles: low by %.2f)"
          % (rs[2]/rs[1], 2/(rs[2]/rs[1])))
    sys.exit(1 if fails else 0)
