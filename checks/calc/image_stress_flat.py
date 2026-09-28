#!/usr/bin/env python3
"""The same energy density where the transverse radius varies, in closed form.

A.19 finds the image energy density vanishing identically, and identifies the reason: the
fold's three pullbacks assemble it into the wave operator, which annihilates the two-point
function.  It also flags the step it could not take, that the assembly was done where the
transverse radius is a constant, and a geometry where the radius varies brings terms the
flat assembly never had to reproduce.

Flat four-dimensional Minkowski is exactly that geometry.  Written in spherical coordinates
it is -dt^2 + dr^2 + r^2 dOmega^2, the transverse radius is r rather than a constant, and
everything is known in closed form, so the step A.19 left open can be taken here rather than
argued about.

THE FOLD.  Theta(t, x) = (T - t, -x): a time reflection about t = T/2 composed with the
spatial point reflection, which is a PT transformation and an isometry and an involution.
Its differential is minus the identity on every Cartesian component, which in spherical
terms is exactly A.19's three pullbacks: the radial basis vector at -x is -n, so the radial
pullback is +1, and the angular ones are -1.  The two descriptions agree.

One defect to state rather than bury: Theta fixes the single event (T/2, 0), so it is not
free, and the real fold is.  Away from the origin that makes no difference to what is
computed here, and the contact surface sits at r > 0.  A second: a straight line is the only
geodesic joining x to its image, so there is no caustic and no refocusing, which is why the
power below is the generic one and not A.18's.  This model is here for the sign, not the
power, and it delivers the generic power as a by-product that checks it.

WHAT COMES OUT.  With all four pullbacks equal to -1,

    rho_img = -(1/2) [ d_t d_t' + sum_i d_i d_i' ] G = (1/2) [ G_titi + laplacian G ],

and the wave equation in the separation turns that into rho_img = d^2 G / d(Delta t)^2
exactly.  Evaluated at the image pair, where Delta x = -2x and Delta t = T - 2t = xi,

    rho_img = ( 8 r^2 + 6 xi^2 ) / ( 4 pi^2 D^3 ),        D = 4 r^2 - xi^2,

which is positive everywhere outside the contact surface and diverges as D goes to zero.
D vanishes linearly in the offset, so the divergence is the generic third power, and the
sign is positive.  Everything is derived symbolically below rather than quoted.
"""
import sympy as sp

t, tp, T, eps = sp.symbols("t t_prime T epsilon", real=True)
x1, x2, x3 = sp.symbols("x_1 x_2 x_3", real=True)
y1, y2, y3 = sp.symbols("y_1 y_2 y_3", real=True)


def wightman(dt, dx):
    """Massless scalar Wightman function in 4D Minkowski, spacelike separation."""
    return 1 / (4 * sp.pi ** 2 * (sum(c ** 2 for c in dx) - dt ** 2))


def assembled_density():
    """rho_img from the definition: mixed derivatives with the fold's pullbacks."""
    X, Y = (x1, x2, x3), (y1, y2, y3)
    # Theta y = (T - t_y, -y); G is evaluated between x and Theta y
    dt = (T - tp) - t
    dx = [-Y[i] - X[i] for i in range(3)]
    G = wightman(dt, dx)
    # the four mixed derivatives, then set y = x and t_y = t
    sub = {y1: x1, y2: x2, y3: x3, tp: t}
    f_tt = sp.diff(G, t, tp).subs(sub)
    f_ii = [sp.diff(G, X[i], Y[i]).subs(sub) for i in range(3)]
    return sp.simplify(sp.Rational(1, 2) * (f_tt + sum(f_ii)))


def closed_form():
    """The claim: (8 r^2 + 6 xi^2) / (4 pi^2 (4 r^2 - xi^2)^3)."""
    r, xi = sp.symbols("r xi", positive=True)
    return (8 * r ** 2 + 6 * xi ** 2) / (4 * sp.pi ** 2 * (4 * r ** 2 - xi ** 2) ** 3), r, xi


if __name__ == "__main__":
    print(__doc__)
    print("=" * 78)

    rho = assembled_density()
    print("\n[1] the assembled density, straight from the mixed derivatives")
    print(f"    rho_img = {sp.simplify(rho)}")

    claim, r, xi = closed_form()
    # put the assembled expression in terms of r and xi: |x| = r, Delta t = T - 2t = xi
    rho_rxi = rho.subs({x1: r, x2: 0, x3: 0}).subs({T: xi + 2 * t})
    diff = sp.simplify(sp.together(rho_rxi - claim))
    print("\n[2] against the closed form quoted in the docstring")
    print(f"    difference = {diff}")
    assert diff == 0, "the closed form is not what the derivatives give"

    print("\n[3] the wave-equation route must give the same thing")
    dt_s, r_s = sp.symbols("Delta_t r", positive=True)
    G_sep = 1 / (4 * sp.pi ** 2 * (4 * r_s ** 2 - dt_s ** 2))
    route = sp.simplify(sp.diff(G_sep, dt_s, 2) - claim.subs({r: r_s, xi: dt_s}))
    print(f"    d^2 G / d(Delta t)^2 minus the closed form = {route}")
    assert route == 0

    print("\n[4] the sign, and it is not a choice")
    print("    numerator 8 r^2 + 6 xi^2 is positive for every r and xi, and D^3 is")
    print("    positive outside the contact surface, so rho_img > 0 there and the")
    print("    divergence on approach is positive.")
    for rv, xv in ((1.0, 1.9), (1.0, 1.99), (1.0, 1.999), (2.0, 3.999)):
        val = float(claim.subs({r: rv, xi: xv}))
        print(f"    r = {rv:.1f}, xi = {xv:.4f}:  D = {4*rv**2 - xv**2:.6f},  "
              f"rho_img = {val:+.6e}")
        assert val > 0

    print("\n[5] the power, as a check that this is the generic case and not A.18's")
    import math
    offs, vals = [], []
    for k in range(6):
        d = 1e-3 / 2 ** k                     # D = 4 r^2 - xi^2, at r = 1
        xv = math.sqrt(4.0 - d)
        offs.append(d)
        vals.append(float(claim.subs({r: 1.0, xi: xv})))
    slope = ((math.log(vals[-1]) - math.log(vals[0]))
             / (math.log(offs[-1]) - math.log(offs[0])))
    print(f"    fitted power in D over five halvings: {slope:+.6f}")
    print(f"    written as the p in D^-p, which is how the text reads it: {abs(slope):.6f}")
    print(f"    the shortfall from -3 is the numerator moving with xi at finite offset,")
    print(f"    not an error: 8 r^2 + 6 xi^2 runs from 31.994 to 32 across the window")
    print(f"    the generic null pair gives -3, A.18's caustic gave -7/2, and there is")
    print(f"    no caustic here because a straight line is the only connecting geodesic")
    assert abs(slope + 3.0) < 1e-3

    print("\n[6] so where the transverse radius varies the cancellation does not happen,")
    print("    and the image energy density is positive and divergent at contact.  A.19's")
    print("    zero belongs to the constant-radius model.  What this does not settle is")
    print("    the black hole: no caustic here, no horizon, and the fold has a fixed")
    print("    point at the origin where the real one has none.")
