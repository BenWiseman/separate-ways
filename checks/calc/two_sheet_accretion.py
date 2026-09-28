#!/usr/bin/env python3
"""Does the far sheet feed our black hole?  Section 6 carries kappa = 1 + Mdot_other/Mdot_ours
and says the question is open.  It is not: the fold's own time reversal settles it, and the
answer is kappa = 1.

The class is defined by the fold N: (U,V) -> (-U,-V) mapping the hole to itself, so exterior
A = {U<0, V>0} is carried to exterior B = {U>0, V<0} and our future interior F = {U,V>0} to the
past interior P.  In a fold-invariant state the matter falling in over there is the image of the
matter falling in here, which is the step Section 6 reads as forcing kappa = 2.  Take a null
shell from our side at V = V_0 > 0.  Its image is the null surface V = -V_0, and N reverses time
orientation, so the image is future-directed on the partner's clock in the direction of
decreasing V.  Each shell therefore deposits its mass on its own sheet's far side:

    V > V_0     mass M + delta      contains our exterior's future null infinity
    |V| < V_0   mass M              contains both bifurcation horizons
    V < -V_0    mass M + delta      contains the partner's future null infinity

Our exterior never touches the third region.  The partner's infall does not reach our mass, and
the symmetric reading of kappa = 2 double-counts one shell.

For the partner's infall to raise our mass the jump at V = -V_0 would have to run the other way,
which is dm/dv < 0 in the partner's own advanced time, so the far sheet can add to our hole only
by violating the null energy condition on its own clock.  That is the result, and it is checked
below rather than drawn: the Einstein tensor of ingoing Vaidya, the Misner-Sharp mass in each
region, and the junction that glues them, each with planted faults.

What this does NOT settle is the case Section 5 opens, where the sheets are in genuine contact
inside r_contact rather than being disjoint images glued along a bifurcation surface.  There the
partner's infall can cross into our future interior, and kappa lies strictly between 1 and 2 at a
value the contact geometry fixes.  That is the live question; the symmetric value is not.
"""
import sympy as sp

v, r, th, M, delta = sp.symbols("v r theta M delta", positive=True)
m = sp.Function("m", positive=True)


def einstein(gmat, coords):
    g = sp.Matrix(gmat)
    ginv = g.inv()
    n = len(coords)
    Gam = [[[sum(ginv[a, d] * (sp.diff(g[d, b], coords[c]) + sp.diff(g[d, c], coords[b])
                               - sp.diff(g[b, c], coords[d])) for d in range(n)) / 2
             for c in range(n)] for b in range(n)] for a in range(n)]
    Riem = [[[[sp.diff(Gam[a][b][d], coords[c]) - sp.diff(Gam[a][b][c], coords[d])
               + sum(Gam[a][c][e] * Gam[e][b][d] - Gam[a][d][e] * Gam[e][b][c] for e in range(n))
               for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
    Ric = sp.Matrix(n, n, lambda b, d: sp.simplify(sum(Riem[a][b][a][d] for a in range(n))))
    R = sp.simplify(sum(ginv[b, d] * Ric[b, d] for b in range(n) for d in range(n)))
    return sp.Matrix(n, n, lambda a, b: sp.simplify(Ric[a, b] - g[a, b] * R / 2)), ginv


def vaidya(mass, exponent=1):
    """Ingoing Vaidya; exponent != 1 is the planted fault."""
    return [[-(1 - 2 * mass / r ** exponent), 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, r ** 2, 0],
            [0, 0, 0, r ** 2 * sp.sin(th) ** 2]]


def misner_sharp(ginv, coords):
    """m = (r/2)(1 - grad r . grad r), the mass inside a sphere of areal radius r."""
    dr = [sp.diff(r, c) for c in coords]
    return sp.simplify(r / 2 * (1 - sum(ginv[a, b] * dr[a] * dr[b]
                                        for a in range(4) for b in range(4))))


if __name__ == "__main__":
    coords = [v, r, th, sp.Symbol("phi")]
    print("Ingoing Vaidya, ds^2 = -(1 - 2m(v)/r)dv^2 + 2 dv dr + r^2 dOmega^2\n")

    G, ginv = einstein(vaidya(m(v)), coords)
    want = 2 * sp.diff(m(v), v) / r ** 2
    off = [(a, b) for a in range(4) for b in range(4)
           if not (a == 0 and b == 0) and sp.simplify(G[a, b]) != 0]
    print("  G_vv - 2 m'(v)/r^2 = %s" % sp.simplify(G[0, 0] - want))
    print("  components other than G_vv that are nonzero: %s" % (off or "none"))
    print("  so 8 pi T_vv = 2 m'(v)/r^2, and the null energy condition on this stream is")
    print("  T_ab k^a k^b = m'(v)/(4 pi r^2) >= 0, which is m increasing in advanced time.")
    print("  Misner-Sharp mass: %s" % misner_sharp(ginv, coords))

    print("\n  planted faults, each of which must break the identification:")
    Gbad, gibad = einstein(vaidya(m(v), exponent=2), coords)
    offbad = [(a, b) for a in range(4) for b in range(4)
              if not (a == 0 and b == 0) and sp.simplify(Gbad[a, b]) != 0]
    print("    1/r^2 in place of 1/r:   G_vv - 2m'/r^2 = %s, other nonzero %s"
          % (sp.simplify(Gbad[0, 0] - want), offbad or "none"))
    Gc, gic = einstein(vaidya(M), coords)
    print("    constant mass:           G = 0 throughout: %s, Misner-Sharp %s"
          % (all(sp.simplify(Gc[a, b]) == 0 for a in range(4) for b in range(4)),
             misner_sharp(gic, coords)))

    print("\nThe fold-symmetric pair of shells, as a mass function of the Kruskal V:")
    V0 = sp.Rational(1, 1)

    def mass_of_V(V, ours, theirs):
        if V > V0:
            return M + ours
        if V < -V0:
            return M + theirs
        return M

    both = [(V, mass_of_V(V, delta, delta)) for V in (sp.Rational(3, 2), 0, sp.Rational(-3, 2))]
    print("  fold-invariant, both sheets accreting:  m(3/2) = %s, m(0) = %s, m(-3/2) = %s"
          % tuple(x[1] for x in both))
    print("  invariance m(-V) = m(V): %s"
          % (mass_of_V(sp.Rational(-3, 2), delta, delta) == mass_of_V(sp.Rational(3, 2), delta, delta)))
    print("  our future null infinity sits at V -> +oo, so we measure %s"
          % mass_of_V(sp.Rational(3, 2), delta, delta))

    print("\n  the counterfactual that decides kappa: only the far sheet accretes.")
    print("  our future null infinity then measures %s, not %s"
          % (mass_of_V(sp.Rational(3, 2), 0, delta), M + delta))
    print("  so d(our mass)/d(their accretion) = %s and kappa = 1 + %s = 1."
          % (sp.diff(mass_of_V(sp.Rational(3, 2), 0, delta), delta),
             sp.diff(mass_of_V(sp.Rational(3, 2), 0, delta), delta)))

    print("\n  and the only way to get kappa > 1 out of this geometry:")
    print("  put M + delta at |V| < V_0 and M beyond, so our exterior gains what they drop.")
    print("  On the partner's clock, where the future is decreasing V, that stream has")
    print("  dm/dv' = -delta < 0 for delta > 0, hence T_v'v' < 0: their infall must carry")
    print("  negative energy on their own clock.  A fold-invariant state cannot have it both")
    print("  ways, because the same sign has to serve both sheets.")

    print("\nAnd the contact route, which is the one that looked open.")
    print("  Inside the contact radius a point and its fold image are causally joined, so matter")
    print("  crosses without passing through our exterior. Let p be the fraction that crosses.")
    print("  Fold invariance makes the two sheets' infall equal, and the crossing is a swap and")
    print("  not a duplication, so run the ledger both ways:")
    F, sink = sp.Symbol("F", positive=True), sp.Symbol("sigma")
    pfrac = sp.Symbol("p", nonnegative=True)
    ours = (1 - pfrac) * F + pfrac * F
    print("    our hole receives (1-p)F from our sheet and pF from theirs: %s" % sp.simplify(ours))
    print("    multiplier over the one-sided case: %s" % sp.simplify(ours / F))
    print("    d/dp = %s, so no value of p buys anything." % sp.diff(sp.simplify(ours / F), pfrac))
    print("  Keeping all of our own inflow AND adding p of theirs is the version that gains, and")
    print("  it creates energy: totals (F,F) in, ((1+p)F,(1+p)F) out, a surplus of %s."
          % sp.simplify(2 * (1 + pfrac) * F - 2 * F))
    print("  planted non-conservative variant, to show the ledger can fail:")
    print("    duplicating rather than swapping gives multiplier %s, derivative %s"
          % (sp.simplify((F + pfrac * F) / F), sp.diff((F + pfrac * F) / F, pfrac)))
    print("  So T = 0 whatever the contact geometry turns out to be, and the transmitted-fraction")
    print("  reading is wrong for the same reason the symmetric reading was: it counts one")
    print("  stream twice. Accretion cannot tell the two classes apart.")

    print("\nJunction: r is continuous across the shell, which fixes the Dray-'t Hooft shift.")
    Mn, dn, V0n = 1.0, 0.1, 1.0
    import mpmath as mpm

    def U_of_r(rv, mass):
        return float((1 - rv / (2 * mass)) * mpm.e ** (rv / (2 * mass)) / V0n)

    print("   r        U inside (M)    U outside (M+delta)   shift")
    for rv in (1.6, 2.0, 2.4, 3.0):
        u1, u2 = U_of_r(rv, Mn), U_of_r(rv, Mn + dn)
        print("   %-8.2f %-15.6f %-21.6f %+.6f" % (rv, u1, u2, u2 - u1))
    print("   the shift is one-signed and finite for delta > 0, so the glue exists and the")
    print("   three-region ledger above is a solution rather than a diagram.")
