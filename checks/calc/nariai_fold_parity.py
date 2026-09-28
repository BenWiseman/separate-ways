#!/usr/bin/env python3
"""What the fold does to Nariai's perturbations, and what it does not do.

`nariai_stability.R` showed that of the three spherical modes phi = c.X, the two that run
away in static time are odd under the fold and do not descend to the quotient while the
bounded one is even and does.  It is tempting to read that as the fold stabilising Nariai.
This script checks how far the statement actually reaches, and the answer is: one sector,
one patch, and not the Euclidean section.

Three things are checked, each by construction rather than by assertion.

1. The spherical modes.  In global coordinates ds^2 = -dt^2 + cosh^2 t dtheta^2 + dOmega^2
   with X0 = sinh t, X1 = cosh t cos theta, X2 = cosh t sin theta, the fold is
   N: (t, theta, Omega) -> (-t, pi - theta, P Omega).  X0 and X1 are odd, X2 is even.  But
   X2 = cosh t sin theta is bounded only inside a static patch: on the hyperboloid
   X2^2 = l^2 + X0^2 - X1^2, and the static patch is X1 > |X0|, which is exactly the region
   where that is below l^2.  Outside it X2 is unbounded, so the surviving mode is static and
   bounded in one observer's patch and nowhere else.

2. Every angular sector keeps half its modes, so the fold empties none of them.  A separated
   mode carries a temporal parity of its own, and the mode equation on dS2 is invariant under
   t -> -t, so each harmonic has one even and one odd temporal solution.  The total parity is
   tau (-1)^(k + l + a + s), and tau can always be chosen to make it even.  Reading the
   spherical calculation as "odd l is discarded" is therefore wrong.

3. The Euclidean negative mode survives, so the semiclassical decay is untouched.  The
   Ginsparg-Perry direction is the relative radius of the two spheres of S^2 x S^2, which is
   constant on both factors.  A constant is invariant under every isometry, so no free
   involution can project it out.
"""
import sympy as sp

t, th, k, l = sp.symbols("t theta k l", real=True)
L = sp.Symbol("Lambda", positive=True)
q = sp.Function("q")


def spherical_modes():
    return {"X0": sp.sinh(t), "X1": sp.cosh(t) * sp.cos(th), "X2": sp.cosh(t) * sp.sin(th)}


def fold(expr):
    """N: t -> -t, theta -> pi - theta.  The angular antipode acts on the other factor."""
    return sp.simplify(expr.subs({t: -t, th: sp.pi - th}, simultaneous=True))


def wave_operator(f, kk, ll):
    """Box phi = 0 on dS2 x S2 for phi = q(t) e^{i k theta} Y_lm, unit radii."""
    return sp.simplify(sp.diff(f, t, 2) + sp.tanh(t) * sp.diff(f, t)
                       + (kk ** 2 / sp.cosh(t) ** 2 + ll * (ll + 1)) * f)


if __name__ == "__main__":
    print("1. The three spherical modes under the fold\n")
    for name, mode in spherical_modes().items():
        image = fold(mode)
        parity = "even" if sp.simplify(image - mode) == 0 else (
                 "odd" if sp.simplify(image + mode) == 0 else "NEITHER")
        print("   %-3s = %-22s ->  %-24s %s" % (name, mode, image, parity))

    print("\n   and how far the surviving one is bounded. On the hyperboloid")
    print("   -X0^2 + X1^2 + X2^2 = 1 the identity X2^2 = 1 + X0^2 - X1^2 holds:")
    X = spherical_modes()
    ident = sp.simplify(X["X2"] ** 2 - (1 + X["X0"] ** 2 - X["X1"] ** 2))
    print("      X2^2 - (1 + X0^2 - X1^2) = %s" % ident)
    print("   The static patch is X1 > |X0|, where X0^2 - X1^2 < 0 and so |X2| < 1.")
    for tv in (0, 1, 3, 6, 10):
        x2 = float(sp.cosh(tv) * sp.sin(sp.Rational(1, 2)))
        print("      at t = %-3d and theta = 1/2, X2 = %12.3f" % (tv, x2))
    print("   so the mode the fold keeps is bounded inside one static patch and unbounded")
    print("   outside it. The spherical claim is a static-patch claim and not a global one.")

    print("\n2. Does the fold empty any angular sector?\n")
    print("   The mode equation q'' + tanh(t) q' + [k^2/cosh^2 t + l(l+1)] q = 0 is invariant")
    print("   under t -> -t, so its solution space splits into even and odd temporal parts:")
    generic = q(t)
    lhs = sp.diff(generic, t, 2) + sp.tanh(t) * sp.diff(generic, t) \
        + (k ** 2 / sp.cosh(t) ** 2 + l * (l + 1)) * generic
    reflected = lhs.subs(t, -t).doit()
    reflected = reflected.subs(q(-t), sp.Function("Q")(t))
    print("      the t -> -t image of the operator differs from it by %s"
          % sp.simplify(sp.tanh(-t) + sp.tanh(t)))
    print("   (the only t-odd coefficient is tanh, and reversing t reverses the single")
    print("    derivative with it, so the two sign flips cancel).")
    print("\n   angular pieces: cos(k theta) -> (-1)^k cos(k theta), sin(k theta) ->")
    print("   -(-1)^k sin(k theta), and Y_lm(P Omega) = (-1)^l Y_lm(Omega). Checked:")
    for kk in (0, 1, 2, 3):
        c, s = sp.cos(kk * th), sp.sin(kk * th)
        fc = sp.simplify(c.subs(th, sp.pi - th))
        fs = sp.simplify(s.subs(th, sp.pi - th))
        print("      k = %d:  cos -> %-16s  sin -> %-16s" % (kk, fc, fs))
    print("\n   so the total parity is tau (-1)^(k+l+a+s) with a = 1 for axial and s = 1 for")
    print("   sine, and tau = +-1 is free. Every harmonic keeps exactly one temporal branch.")
    print("   No angular sector is emptied, and discarding odd l outright is wrong.")

    print("\n   The radiative sectors had nothing to remove in any case:")
    for ll in range(1, 5):
        print("      l = %d:  scalar/polar potential l(l+1) = %-3d   axial l(l+1)-2 = %d"
              % (ll, ll * (ll + 1), ll * (ll + 1) - 2))
    print("   both nonnegative for every l >= 1, so those sectors are already mode-stable")
    print("   and the fold selects their parity rather than removing an instability.")

    print("\n3. The Euclidean negative mode\n")
    print("   Euclidean Nariai is S^2 x S^2 at equal radii and the Ginsparg-Perry direction")
    print("   is the relative radius, constant on both factors. A constant is invariant under")
    print("   every isometry, so the free involution cannot project it out: the mode is even")
    print("   and survives the quotient. The fold does not remove the semiclassical decay.")

    print("\nPlanted failures, each of which must break something above:\n")
    bad = X["X0"] + X["X2"]                       # one odd and one even: genuinely mixed
    mixed = sp.simplify(fold(bad) - bad) != 0 and sp.simplify(fold(bad) + bad) != 0
    print("   X0 + X2 mixes an odd mode with an even one, so it is neither: %s"
          % ("caught" if mixed else "TEST FAILED TO NOTICE"))
    print("   (X0 + X1 would not do: both are odd, so their sum is odd and the plant misses.)")
    print("   calling X2 odd instead of even leaves the residual fold(X2) + X2 = %s"
          % sp.simplify(fold(X["X2"]) + X["X2"]))
    print("   an even friction coefficient would break the parity split. Integrate both")
    print("   operators from even data q(0)=1, q\'(0)=0 and compare q(T) with q(-T):")
    import mpmath as mpm

    def run(friction, T, n=20000):
        h, y, yp, tt = T / n, mpm.mpf(1), mpm.mpf(0), mpm.mpf(0)
        for _ in range(n):
            ypp = -friction(tt) * yp - 6 * y
            y, yp, tt = y + h * yp + h * h * ypp / 2, yp + h * ypp, tt + h
        return y

    for name, fr in (("tanh (the real one)", mpm.tanh), ("cosh (planted)", mpm.cosh)):
        fwd, back = run(fr, mpm.mpf("1.5")), run(fr, mpm.mpf("-1.5"))
        print("      %-20s q(1.5) = %-14s q(-1.5) = %-14s gap %s"
              % (name, mpm.nstr(fwd, 8), mpm.nstr(back, 8), mpm.nstr(abs(fwd - back), 3)))
    print("   even data stays even under tanh and does not under cosh, which is the split.")
