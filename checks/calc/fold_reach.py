#!/usr/bin/env python3
"""How far is a point from its own fold image?  The field, not just its zero.

The paper uses the separation between a point and its fold image in one place only, to find
where it vanishes.  That is the contact result.  But the separation is a scalar FIELD on the
whole spacetime, determined by the geometry and the fold with nothing chosen, and away from its
zero it has a profile and a gradient.  If the fold has any local effect at all, this field is what
carries it, so it is worth computing rather than leaving implicit.

On the $t=0$ slice of Kruskal the geometry is the Einstein-Rosen bridge, $d\\ell^2 + r(\\ell)^2
d\\Omega^2$ with $\\ell$ the proper radial distance and $\\ell=0$ at the throat $r=r_h$.  The fold
sends a static observer at $(r,\\theta)$ in one exterior to $(r,\\pi-\\theta)$ in the other at the
same Killing time, so the separation is the length of the shortest curve that crosses the bridge
AND turns through $\\pi$ on the sphere.  That is a geodesic on a surface of revolution, so
Clairaut applies: $r\\sin\\psi = L$ is conserved, and

    ds = dl / sqrt(1 - L^2/r^2),      dtheta = L dl / (r^2 sqrt(1 - L^2/r^2)).

For each starting radius, solve for the $L$ that sweeps exactly $\\pi$, then read off the length.

Two things come out.  On this slice the separation never falls below $\\pi r_h$, half the horizon
circumference, because the angular half-turn has to be paid for and the sphere is smallest at the
throat; that is a bound on spacelike separation across the bridge and not on the interior, where
the contact result makes a point and its image causally connected.  And the field is monotone outward with a nonzero radial gradient, so the fold does
define a preferred direction and a scale at every point of the exterior, both set by the hole's
mass alone.  What the field does not do is act: the separation is a geometric fact, and turning it
into a local effect needs a term in an effective action that depends on it, which this paper does
not have.  Recording the field is what makes that a well-posed question instead of a feeling.
"""
import mpmath as mpm

mpm.mp.dps = 20
RH = mpm.mpf(1)                                  # horizon radius, units 2M = 1


def half_integrals(r0, L):
    """Angle swept and length, from the throat out to r0.

    dl = dr / sqrt(1 - rh/r) has an integrable endpoint singularity at the throat; the
    substitution r = rh/(1 - u^2) removes it exactly, since then sqrt(1 - rh/r) = u and
    dl = 2 rh du/(1 - u^2)^2.  Nothing is regulated by hand.
    """
    u0 = mpm.sqrt(1 - RH / r0)

    def pieces(u):
        r = RH / (1 - u ** 2)
        dl = 2 * RH / (1 - u ** 2) ** 2
        w = mpm.sqrt(1 - L ** 2 / r ** 2)
        return dl / w, L * dl / (r ** 2 * w)

    length = mpm.quad(lambda u: pieces(u)[0], [0, u0])
    angle = mpm.quad(lambda u: pieces(u)[1], [0, u0])
    return angle, length


def l_of_r(r0):
    u0 = mpm.sqrt(1 - RH / r0)
    return mpm.quad(lambda u: 2 * RH / (1 - u ** 2) ** 2, [0, u0])


def separation(r0):
    """Length of the geodesic from (r0, 0) on one sheet to (r0, pi) on the other."""
    f = lambda L: 2 * half_integrals(r0, L)[0] - mpm.pi
    lo, hi = mpm.mpf("1e-8"), RH * (1 - mpm.mpf("1e-10"))
    L = mpm.findroot(f, (lo, hi), solver="bisect", tol=mpm.mpf("1e-14"))
    return 2 * half_integrals(r0, L)[1], L


if __name__ == "__main__":
    print("Separation from a static observer to its fold image, units 2M = 1.\n")
    print("  the angular half-turn alone, taken at the throat, costs pi r_h = %s"
          % mpm.nstr(mpm.pi * RH, 8))
    print("  so on this slice nothing is nearer its image than half a horizon circumference.")
    print("  That is a bound on SPACELIKE separation across the bridge. It says nothing about")
    print("  the interior, where the contact result has a point and its image causally")
    print("  connected at r = M, hence null or timelike separation. Different region,")
    print("  different sign, and the two must not be quoted as one statement.\n")
    print("   r        proper depth l(r)   separation s(r)    Clairaut L    s - pi r_h")
    prev, monotone = None, True
    rows = {}
    for label in ("1.001", "1.05", "1.2", "1.5", "2", "3", "5", "10"):
        r0 = mpm.mpf(label)
        s_, L = separation(r0)
        rows[label] = s_
        print("   %-8s %-19s %-18s %-13s %s"
              % (label, mpm.nstr(l_of_r(r0), 8), mpm.nstr(s_, 8),
                 mpm.nstr(L, 6), mpm.nstr(s_ - mpm.pi * RH, 6)))
        if prev is not None and s_ <= prev:
            monotone = False
        prev = s_

    print("\n  Checks that can fail:")
    print("    monotone outward:                              %s" % monotone)
    print("    s(r) >= pi r_h at every radius sampled:        %s"
          % all(v >= mpm.pi * RH - mpm.mpf("1e-9") for v in rows.values()))
    print("    s -> pi r_h as r -> r_h:  s(1.001) - pi = %s"
          % mpm.nstr(rows["1.001"] - mpm.pi * RH, 4))
    straight = 2 * l_of_r(mpm.mpf(3)) + mpm.pi * RH
    print("    the broken path, straight in then a turn at the throat, gives %s at r = 3,"
          % mpm.nstr(straight, 8))
    print("    which must exceed the geodesic %s: %s"
          % (mpm.nstr(rows["3"], 8), straight > rows["3"]))
    print("    planted: demanding a quarter turn instead of half must give a shorter curve:")
    g = lambda L: 2 * half_integrals(mpm.mpf(3), L)[0] - mpm.pi / 2
    Lq = mpm.findroot(g, (mpm.mpf("1e-8"), RH * (1 - mpm.mpf("1e-10"))), solver="bisect")
    sq = 2 * half_integrals(mpm.mpf(3), Lq)[1]
    print("       quarter turn at r = 3 gives %s against %s for the half turn: %s"
          % (mpm.nstr(sq, 8), mpm.nstr(rows["3"], 8), sq < rows["3"]))

    print("\n  And the comparison that decides whether any of this is a local structure at all:")
    print("  in EXACT de Sitter the same field is constant. The fold's invariant there is")
    print("  Z(x, Ax) = X.(-X)/l^2 = -1 at every point of the hyperboloid, checked below at")
    print("  random points, so the separation to one's image is the same everywhere and its")
    print("  gradient vanishes identically. There is no local fold structure in empty de Sitter.")
    import random
    random.seed(11)
    worst = mpm.mpf(0)
    for _ in range(6):
        # a point on -X0^2 + X1^2 + X2^2 + X3^2 + X4^2 = 1
        t_, a_, b_ = [mpm.mpf(random.uniform(-2, 2)) for _ in range(3)]
        sp = mpm.sqrt(1 + t_ ** 2)
        X = [t_, sp * mpm.cos(a_), sp * mpm.sin(a_) * mpm.cos(b_),
             sp * mpm.sin(a_) * mpm.sin(b_), mpm.mpf(0)]
        eta = [-1, 1, 1, 1, 1]
        Z = sum(e * x * (-x) for e, x in zip(eta, X))
        worst = max(worst, abs(Z + 1))
    print("     worst |Z(x, Ax) + 1| over six random points: %s" % mpm.nstr(worst, 4))
    print("     planted: the map X -> -X with one sign left unflipped gives Z = %s, not -1."
          % mpm.nstr(sum(e * x * y for e, x, y in
                         zip([-1, 1, 1, 1, 1], X, [X[0]] + [-v for v in X[1:]])), 6))
    print("\n  So the field is switched on by mass and is identically flat without it, which is")
    print("  the sharpest form of the question: whatever the fold does locally, it does in")
    print("  proportion to how far the geometry departs from the empty case.")

    print("\n  How the reach approaches its floor, which is the clock.")
    print("  For a static observer the redshift factor is sqrt(f) with f = 1 - rh/r, so f is")
    print("  the SQUARE of the local clock rate. Near the horizon the excess separation over")
    print("  the floor is proportional to it, with a pure-number coefficient:")
    u = mpm.pi / (2 * mpm.sqrt(2))
    closed = 2 * mpm.sqrt(2) * mpm.coth(u)
    print("     s - pi rh  ->  C rh f      with C = 2 sqrt(2) coth(pi/(2 sqrt 2)) = %s"
          % mpm.nstr(closed, 12))
    print("     eps        f               (s - pi rh)/f      minus C")
    resid = []
    for e in ("1e-2", "3e-3", "1e-3", "3e-4", "1e-4"):
        eps = mpm.mpf(e)
        r0 = RH * (1 + eps)
        f = 1 - RH / r0
        sv, _ = separation(r0)
        ratio = (sv - mpm.pi * RH) / f
        resid.append(ratio - closed)
        print("     %-10s %-15s %-18s %s"
              % (e, mpm.nstr(f, 6), mpm.nstr(ratio, 12), mpm.nstr(ratio - closed, 4)))
    print("     the residual falls linearly in eps, ratios %s, so the limit is the closed"
          % ", ".join(mpm.nstr(resid[i - 1] / resid[i], 6) for i in range(1, len(resid))))
    print("     form and not a fit: eps drops by 10/3 then 3 alternately and so does it.")
    print("     planted: 9 pi/8 = %s would show residuals that do NOT go to zero, %s at the"
          % (mpm.nstr(9 * mpm.pi / 8, 10),
             mpm.nstr(abs(closed - 9 * mpm.pi / 8), 4)))
    print("     smallest eps, which is what a near-miss closed form looks like.")
    print("\n  So the fold's reach above its floor vanishes as the square of the clock rate.")
    print("  Where the clock stops the reach is at its floor, and that place is the horizon.")
    print("  The two are not separate facts about the geometry; they are one quantity.")

    print("\n  What a coupling that falls with separation would look like across the exterior.")
    print("  A massless two-point function goes as 1/s^2 at short range, so an image")
    print("  correlation evaluated between a point and its own image would carry 1/s(r)^2.")
    print("  This is conditional on such a coupling existing at all, which is open, and is")
    print("  printed as a profile rather than a prediction:")
    floor = mpm.pi * RH
    print("     r        s(r)          1/s^2        relative to the horizon value")
    ref = 1 / floor ** 2
    for label in ("1.001", "1.5", "2", "3", "5", "10"):
        v = rows[label]
        print("     %-8s %-13s %-12s %s"
              % (label, mpm.nstr(v, 7), mpm.nstr(1 / v ** 2, 6), mpm.nstr((1 / v ** 2) / ref, 5)))
    print("  Two features, and both are Ben's picture rather than mine. It strengthens toward")
    print("  mass rather than away from it, because s shrinks as the horizon is approached.")
    print("  And it SATURATES: s cannot fall below pi r_h on this slice, so the correlation")
    print("  cannot exceed 1/(pi r_h)^2 = %s, a ceiling set by the horizon's own size and by"
          % mpm.nstr(ref, 7))
    print("  nothing else. A grip that tightens with depth and stops at a bound the geometry")
    print("  fixes is the shape he asked for; whether anything couples this way is item 13.")

    print("\n  The gradient is radial and nonzero, and its scale is the hole's mass alone.")
    print("  A geometric fact, not yet a force: making it act needs a term in an effective")
    print("  action that depends on it, and this paper does not have one.")
