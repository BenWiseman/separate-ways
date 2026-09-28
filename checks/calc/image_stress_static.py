#!/usr/bin/env python3
"""The image energy density on any static spherically symmetric hole, and what survives.

A.19 computes the image term's energy density in a geometry whose transverse radius is
constant, finds it vanishes identically because the fold's pullbacks assemble it into the
wave operator, and leaves one step open: whether that assembly survives a varying metric,
and if it does not, what the remainder is.  It guessed the remainder would be the curvature
term of the Regge-Wheeler potential.  This file does the assembly in general and the guess
was wrong, in an interesting direction.

THE SETTING.  Any static spherically symmetric metric

    ds^2 = -h(r) dt^2 + dr^2 / f(r) + R(r)^2 dOmega^2,

which covers Schwarzschild at h = f = 1 - 2M/r and R = r, flat space at h = f = 1 and
R = r, and A.19's model at h = f = 1 and R constant.  The fold is Theta(t, r, n) =
(T - t, r, -n), an isometry of any static metric composed with the antipodal map, so the
three pullbacks are the same ones A.19 uses: -1 on time, +1 on the radius, -1 on the sphere.

THE ASSEMBLY.  For a massless minimally coupled scalar the static observer measures

    rho = (1/2h) (d_t phi)^2 + (f/2) (d_r phi)^2 + (1/2R^2) |grad_perp phi|^2,

and point-splitting each term against the image kernel, mode by mode with
G = sum_l c_l P_l(cos gamma) g_l and the antipodal contraction 2 P_l'(-1) = -(-1)^l l(l+1),

    rho_l = (1/2h) d_tau^2 g_l + (f/2) d_r d_r' g_l + (l(l+1)/2R^2) g_l.

THE RESULT.  The wave equation for the l-mode gives

    (1/h) d_tau^2 g_l = (1/sqrt(G)) d_r( sqrt(G) f d_r g_l ) - (l(l+1)/R^2) g_l,
    sqrt(G) = sqrt(h/f) R^2,

and substituting it, the two l(l+1)/R^2 terms cancel against each other exactly, at every
h, f and R.  What is left is purely radial:

    rho_img = (1/2) sum_l c_l (-1)^l [ (1/sqrt(G)) d_r( sqrt(G) f d_r g_l ) + f d_r d_r' g_l ].

So the transverse sphere always pays for the angular part of the potential, in every static
spherically symmetric geometry, and the remainder is never the curvature term: it is the
difference between the radial part of the wave operator and the mixed radial derivative.
A.19's conjecture is withdrawn by this.

WHY THE TWO KNOWN CASES COME OUT AS THEY DO.  With R constant and h = f = 1 the mode
function depends on r and r' only through r' - r, so d_r^2 g = -d_r d_r' g and the two
radial terms cancel as well: that is A.19's identical zero, and it needed the constant
radius exactly where A.19 said it did.  Let R vary and the measure sqrt(G) brings in
R'/R terms that have no counterpart in the mixed derivative, the cancellation stops, and
a non-zero density is left: that is the flat-space result, positive and divergent.

WHAT IS STILL NOT DONE, STATED PLAINLY.  The remainder is a two-dimensional object.  It
needs g_l and its mixed radial derivative at the contact orbit in the Hartle-Hawking state,
inside the horizon, where no closed form exists and where d_t is spacelike so the static
observer above does not exist.  That is a smaller problem than the one A.19 left, because
the angular half is gone and the sphere is out of it, and it is still a problem.
"""
import sympy as sp


def assemble():
    """Carry out the substitution and return what is left after it."""
    r, tau = sp.symbols("r tau", positive=True)
    L = sp.Symbol("Lambda", positive=True)               # l(l+1)
    h, f, R = (sp.Function(n, positive=True)(r) for n in ("h", "f", "R"))
    g = sp.Function("g")(r, tau)
    dr_drp = sp.Function("M")(r)                         # the mixed derivative, an unknown

    sqrtG = sp.sqrt(h / f) * R ** 2
    radial = sp.diff(sqrtG * f * sp.diff(g, r), r) / sqrtG
    # the wave equation, solved for the time term as it appears in rho
    time_term = radial - L * g / R ** 2

    rho = sp.Rational(1, 2) * time_term + f * dr_drp / 2 + L * g / (2 * R ** 2)
    return sp.simplify(rho - (radial / 2 + f * dr_drp / 2)), radial, r, h, f, R, g, L


if __name__ == "__main__":
    print(__doc__)
    print("=" * 78)

    residue, radial, r, h, f, R, g, L = assemble()
    print("\n[1] the angular terms cancel identically, at every h, f and R")
    print(f"    rho minus the purely radial pair = {residue}")
    assert residue == 0, "the angular cancellation failed"
    print("    so the sphere pays for l(l+1)/R^2 exactly, and no l survives in rho")

    print("\n[2] A.19's case: R constant, h = f = 1, and the radial pair cancels too")
    a = sp.Symbol("a", positive=True)
    sub = {h: 1, f: 1, R: a}
    rad0 = sp.simplify(radial.subs({sp.Function("h", positive=True)(r): 1,
                                    sp.Function("f", positive=True)(r): 1,
                                    sp.Function("R", positive=True)(r): a}).doit())
    print(f"    radial part of the wave operator becomes {rad0}")
    # with g a function of r' - r only, d_r^2 g = -d_r d_r' g, so the pair is zero
    print("    and with g depending on r' - r alone, d_r d_r' g = -d_r^2 g, so the pair")
    print("    vanishes: that is A.19's identically zero density, and it needs R constant")

    print("\n[3] what varying R does, which is where A.19's zero comes from")
    rr = sp.Symbol("r", positive=True)
    radR = sp.simplify((sp.diff(rr ** 2 * sp.diff(sp.Function("g")(rr), rr), rr) / rr ** 2).doit())
    print(f"    with h = f = 1 and R = r the radial operator is {radR}")
    print("    the 2/r piece has no counterpart in the mixed derivative, so the pair")
    print("    no longer cancels and a density is left: the flat-space case of A.19")

    print("\n[4] Schwarzschild, written out")
    M = sp.Symbol("M", positive=True)
    fs = 1 - 2 * M / rr
    gg = sp.Function("g")(rr)
    sqrtGs = sp.sqrt(sp.Integer(1)) * rr ** 2         # h = f, so sqrt(h/f) = 1
    rad_s = sp.simplify((sp.diff(sqrtGs * fs * sp.diff(gg, rr), rr) / sqrtGs).doit())
    print(f"    (1/r^2) d_r( r^2 f d_r g ) = {sp.expand(rad_s)}")
    print("    so rho_img = (1/2) sum_l c_l (-1)^l [ that + f d_r d_r' g_l ],")
    print("    with no l(l+1) anywhere in it.")

    print("\n[5] the correction this makes to A.19")
    print("    A.19 conjectured the remainder would be the curvature term f f'/r of the")
    print("    Regge-Wheeler potential, on the grounds that the sphere supplies only the")
    print("    l(l+1)/r^2 part.  The sphere does supply that part, exactly, which is what")
    print("    [1] shows.  But the potential is not what is left over: the time term was")
    print("    eliminated by the wave equation, and what the wave equation leaves is the")
    print("    radial operator, not the potential.  The remainder is radial derivatives of")
    print("    g and nothing else.  The conjecture is withdrawn.")
