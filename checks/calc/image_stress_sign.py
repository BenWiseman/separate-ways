#!/usr/bin/env python3
"""What the image term does to the energy density, and why this model cannot say.

A.18 settles the power of the image divergence at the contact sphere and leaves the sign,
which A.15 says is the half that decides whether the divergence closes the region or only
marks it.  This file goes after the sign in the same model and returns a negative result
with a reason, which is worth more than a number would have been if the number were an
artifact of the model.

THE ASSEMBLY.  Same geometry as A.18, ds^2 = -dt^2 + dr^2 + a^2 dOmega^2 with a = 1, and a
massless minimally coupled scalar, whose energy density for a static observer is

    rho = (1/2) [ (d_t phi)^2 + (d_r phi)^2 + |grad_perp phi|^2 ],

every term positive classically.  Point-split with the fold's image kernel: the fold here
is Theta(t, r, n) = (T - t, r, -n), which is an isometry, an involution, and free, and
which reverses time orientation the way the real one does.  Writing F_ab for the mixed
derivative of G(x, Theta y) at y = x, the three pullbacks are

    d/dt_y  = -d/dt'        (Theta reverses t)
    d/dr_y  = +d/dr'        (Theta fixes r)
    d/dn_y  = -d/dn'        (the antipodal map's differential is minus the identity)

so rho_img = (1/2)(F_tt + F_rr + F_ang) with F_tt = +G_tautau, F_rr = -G_rhorho and
F_ang = -g^AB d_A d_B' G, all evaluated at the image pair: time separation xi, radial
separation zero, transverse angle pi.

THE ANGULAR PIECE.  At the antipode the first term of the mixed derivative of
P_l(cos gamma) drops, because the projection of n' off n vanishes when n' = -n, and what
is left is 2 P_l'(-1) = -(-1)^l l(l+1).  So the angular contraction carries exactly the
mass of the tower, which is what makes the cancellation below happen.

WHAT COMES OUT.  Mode by mode the three pieces are W_tautau, -W_rhorho and +l(l+1)/a^2 W,
and that combination is the wave operator, which annihilates the two-point function.  The
fold's own pullback signs are what assemble it: reverse time and the time piece changes
sign, take the antipodal map and the angular contraction collapses to the Laplacian
eigenvalue.  So the image term carries no energy density at all here, at every separation
and not only at the caustic, and the divergence A.18 found in the two-point function is
invisible in rho.  Checked to 2e-12 against pieces of order 1.9e9.

That is not the fold carrying nothing.  It is this geometry having nothing for the residue
to be, because the transverse sphere supplies the whole of the radial equation's potential
and there is no remainder.  Give the tower any mass the sphere does not supply and exactly
that much survives:

    rho_img = -(1/2) mu^2 sum_l c_l (-1)^l W,

verified against the assembled sum to 6e-8 at mu = 2.  So in this model the sign of the
energy density is the sign of minus the non-geometric part of the potential, times the A.18
sum, and that sum is real and positive.

One step of the transfer is NOT shown here and should not be assumed.  The assembly is done
where the two-dimensional factor is flat, so the point-split radial term is a plain second
derivative.  On a hole g^rr = f varies, the covariant wave operator carries first-derivative
terms this assembly never had to reproduce, and whether the same three pullbacks still
deliver the wave operator exactly is open.  The constant-mu test shows a remainder in the
potential survives untouched; it says nothing about a varying metric.  If it does carry, the
remainder on Schwarzschild is the curvature term f f'/r of the Regge-Wheeler potential and
computing it on the contact orbit is what is left.

The checks below are built so that this can fail.  The derivatives come from Bessel
identities rather than from the wave equation, so the cancellation is evidence and not a
tautology; each of the three pullback signs is flipped in turn and each flip must move the
sum by exactly the piece it flipped; dropping the time reversal from Theta, which leaves a
perfectly good isometry, must leave a non-zero density; and the mass test above must return
the residue it predicts rather than any residue at all.
"""
import numpy as np
from scipy.special import eval_legendre, hankel2

A = 1.0


def mass(l, mu=0.0):
    """2D mass of the tower.  The sphere supplies l(l+1)/a^2; mu^2 is whatever else is
    in the radial equation and is not geometric, which is the point of the test below."""
    return np.sqrt(l * (l + 1.0) / A ** 2 + mu ** 2)


def w2(m, tau, rho, eps):
    """2D Wightman function at timelike separation, as a function of tau and rho."""
    s = np.sqrt((tau - 1j * eps) ** 2 - rho ** 2 + 0j)
    return -0.25j * hankel2(0, m * s)


def d2(f, h):
    """Second derivative from a five-point stencil, for cross-checking the closed forms."""
    fm2, fm1, f0, fp1, fp2 = f
    return (-fm2 + 16 * fm1 - 30 * f0 + 16 * fp1 - fp2) / (12 * h * h)


def derivatives(m, xi, eps):
    """W, d^2W/dtau^2 and d^2W/drho^2 at rho = 0, from Bessel identities and not from
    the wave equation.  With H0' = -H1 and H1' = H0 - H1/z,

        W      = -(i/4) H0(z),
        W_tt   = -(i/4) m^2 ( -H0(z) + H1(z)/z ),
        W_rr   = -(i/4) m H1(z) / tau,          since d^2 u/d rho^2 = -1/tau at rho = 0,

    with z = m tau and tau carrying the i-epsilon.  Nothing here knows the mode satisfies
    a wave equation, which is what makes the cancellation below evidence.
    """
    tau = xi - 1j * eps
    z = m * tau
    h0, h1 = hankel2(0, z), hankel2(1, z)
    w = -0.25j * h0
    w_tt = -0.25j * m ** 2 * (-h0 + h1 / z)
    w_rr = -0.25j * m * h1 / tau
    return w, w_tt, w_rr


def pieces(xi, eps, lmax, signs=(-1.0, +1.0, -1.0), mu=0.0):
    """F_tt, F_rr, F_ang at the image pair.  signs are the three Theta pullbacks."""
    st, sr, sa = signs
    l = np.arange(1, lmax + 1, dtype=float)
    m = mass(l, mu)
    c = (2.0 * l + 1.0) / (4.0 * np.pi)
    par = (-1.0) ** l                       # P_l(cos pi)
    w, w_tt, w_rr = derivatives(m, xi, eps)

    # d_t d_t' G = -W_tt and d_r d_r' G = -W_rr, then the pullback signs
    f_tt = np.sum(c * par * st * (-w_tt))
    f_rr = np.sum(c * par * sr * (-w_rr))
    # g^AB d_A d_B' P_l at gamma = pi is 2 P_l'(-1) = -(-1)^l l(l+1)
    # the sphere supplies only its own eigenvalue, whatever else sits in the radial
    # equation; with mu = 0 that is the whole 2D mass and with mu > 0 it is not
    f_ang = np.sum(c * sa * (-par * l * (l + 1.0) / A ** 2) * w)
    return f_tt, f_rr, f_ang


def density(xi, eps=None, lmax=None, signs=(-1.0, +1.0, -1.0), mu=0.0):
    if eps is None:
        eps = (np.pi - xi) / 40.0
    if lmax is None:
        lmax = int(40.0 / eps)
    f_tt, f_rr, f_ang = pieces(xi, eps, lmax, signs, mu)
    return 0.5 * (f_tt + f_rr + f_ang), (f_tt, f_rr, f_ang)


def bare_sum(xi, eps=None, lmax=None, mu=0.0):
    """The A.18 object itself, sum_l c_l (-1)^l W, with no derivatives on it."""
    if eps is None:
        eps = (np.pi - xi) / 40.0
    if lmax is None:
        lmax = int(40.0 / eps)
    l = np.arange(1, lmax + 1, dtype=float)
    c = (2.0 * l + 1.0) / (4.0 * np.pi)
    w, _, _ = derivatives(mass(l, mu), xi, eps)
    return np.sum(c * (-1.0) ** l * w)


if __name__ == "__main__":
    print(__doc__)
    print("=" * 78)

    print("\n[0] the closed form for the angular contraction, against a stencil")
    ls = np.arange(1, 31, dtype=float)
    hx = 1e-5
    stencil = (eval_legendre(ls, -1.0 + hx) - eval_legendre(ls, -1.0 - hx)) / (2 * hx)
    closed = -0.5 * (-1.0) ** ls * ls * (ls + 1.0)
    err = np.max(np.abs(stencil - closed) / np.abs(closed))
    print(f"    P_l'(-1) = (-1)^(l+1) l(l+1)/2 for l up to 30, worst relative error {err:.2e}")
    assert err < 1e-6

    print("\n[0b] the contraction itself, by differences on the sphere rather than by the")
    print("     projection argument, since it is what the cancellation turns on")

    def unit(v):
        return v / np.linalg.norm(v)

    def contraction(l, h=1e-3):
        n, nb = np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, -1.0])
        basis = [np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])]
        tot = 0.0
        for e in basis:
            pp = eval_legendre(l, np.dot(unit(n + h * e), unit(nb + h * e)))
            pm = eval_legendre(l, np.dot(unit(n + h * e), unit(nb - h * e)))
            mp = eval_legendre(l, np.dot(unit(n - h * e), unit(nb + h * e)))
            mm = eval_legendre(l, np.dot(unit(n - h * e), unit(nb - h * e)))
            tot += (pp - pm - mp + mm) / (4 * h * h)
        return tot

    worst = 0.0
    for l in (1, 2, 3, 5, 8, 12):
        closed = -((-1.0) ** l) * l * (l + 1.0)
        rel = abs(contraction(l) - closed) / abs(closed)
        worst = max(worst, rel)
    print(f"     equals -(-1)^l l(l+1), the sphere's own eigenvalue, for l up to 12,")
    print(f"     worst relative error {worst:.1e}")
    assert worst < 1e-4

    print("\n[1] the three pieces at the fold's own pullback signs")
    for s in (0.02, 0.005, 0.001):
        rho, (a_, b_, c_) = density(np.pi - s)
        scale = max(abs(a_), abs(b_), abs(c_))
        print(f"    offset {s:<7.4f} F_tt {a_.real:+13.2f}  F_rr {b_.real:+13.2f}  "
              f"F_ang {c_.real:+13.2f}   |rho|/scale {abs(rho) / scale:.2e}")
        assert abs(rho) / scale < 1e-6, "the cancellation failed; check the assembly"
    print("    each piece is large and the sum is machine zero: rho_img vanishes")

    print("\n[2] flip one pullback sign at a time, and it must stop vanishing")
    print("    each flip changes the sum by twice the piece it flips, so that is the")
    print("    scale it is judged against; the radial piece is small and still must show")
    for idx, name in ((0, "time"), (1, "radial"), (2, "angular")):
        sg = [-1.0, +1.0, -1.0]
        sg[idx] = -sg[idx]
        rho, parts = density(np.pi - 0.005, signs=tuple(sg))
        flipped = abs(parts[idx])
        # The ratio printed to four places showed one part in 1e4 while the assertion below
        # enforces one part in 1e6, which is what the paper states. The evidence for the
        # stated precision has to be in the output, or the claim cannot be checked against
        # the script that backs it.
        dev = abs(abs(rho) / flipped - 1.0)
        print(f"    {name:<8} flipped: |rho| = {abs(rho):.4g}, the piece itself "
              f"{flipped:.4g}, ratio {abs(rho) / flipped:.4f}, departing from 1 by "
              f"{dev:.2e} against a tolerance of one part in 1.0e+06")
        assert dev < 1e-6, f"flipping {name} did not register"

    print("\n[3] drop the time reversal from Theta and the density comes back")
    print("    (antipodal alone is still an isometry and still an involution)")
    rho, (a_, b_, c_) = density(np.pi - 0.005, signs=(+1.0, +1.0, -1.0))
    print(f"    rho = {rho.real:+.3f}, against pieces of order {max(abs(a_), abs(c_)):.1f}")

    print("\n[4] what the zero is really saying: give the tower a mass the sphere does")
    print("    not supply, and exactly that much survives")
    for mu in (0.0, 0.5, 1.0, 2.0):
        rho, _ = density(np.pi - 0.005, mu=mu)
        predicted = -0.5 * mu ** 2 * bare_sum(np.pi - 0.005, mu=mu)
        if mu == 0.0:
            print(f"    mu = {mu:.1f}:  rho = {rho.real:+13.4f},  nothing to survive, "
                  f"and nothing does")
            assert abs(rho) < 1e-3
            continue
        rel = abs(rho - predicted) / abs(predicted)
        print(f"    mu = {mu:.1f}:  rho = {rho.real:+13.4f},  "
              f"-mu^2/2 times the A.18 sum = {predicted.real:+13.4f},  agree to {rel:.1e}")
        assert rel < 1e-4

    print("\n[5] so the zero is not a property of the fold, it is this model having")
    print("    nothing for the residue to be.  The fold's three pullbacks assemble the")
    print("    energy density into the wave operator, which annihilates the two-point")
    print("    function; what survives in any other geometry is the part of the radial")
    print("    potential the transverse sphere does not supply.  On Schwarzschild that")
    print("    is the curvature term f f'/r of the Regge-Wheeler potential, and the")
    print("    sign of the energy density is the sign of that term times the A.18 sum.")
    print("    That sum is real and positive.  The curvature term is not computed here.")
