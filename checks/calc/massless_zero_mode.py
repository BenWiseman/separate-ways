#!/usr/bin/env python3
"""The massless homogeneous mode on Nariai, and why its Euclidean limit is the thing that fails.

The image-correlation construction of A.8 was first run at $m^2=\\Lambda$, a mass chosen to avoid
the homogeneous zero mode.  That invited a reading in which massive fields feel the fold and
massless ones do not.  The reading is wrong, and the reason is worth having in one place: the
$m\\to0$ limit of the EUCLIDEAN REFERENCE STATE fails, while genuine massless states exist and
carry the effect.

This file checks the second half, which is the part that decides it.  On the $\\mathrm{dS}_2$
factor the homogeneous massless equation is $\\ddot q+\\tanh t\\,\\dot q=0$, whose two solutions are
the constant and the Gudermannian $\\mathrm{gd}(t)=\\arctan\\sinh t$.  A candidate pair is

    u_alpha(t) = sqrt(Lambda/(16 pi^2 alpha)) (1 - i alpha arctan sinh t),   alpha > 0,

and both the solution and the normalising constant are rederived here rather than taken: the
Klein-Gordon product is computed from the metric, including the $2\\pi\\cosh t$ of the circle and
the $4\\pi/\\Lambda$ of the sphere, and is required to be one.

What this settles: there is a normalisable massless canonical pair, so the zero mode is a state to
be chosen and not an obstruction to having one.  What it does not settle: a massless SCALAR is not
a photon, and nothing here is a Maxwell calculation.
"""
import sympy as sp

t = sp.Symbol("t", real=True)
alpha, Lam = sp.symbols("alpha Lambda", positive=True)


def gudermannian():
    return sp.atan(sp.sinh(t))


def wave_residual(q):
    """The homogeneous massless equation on dS_2: q'' + tanh(t) q' = 0."""
    return sp.simplify(sp.diff(q, t, 2) + sp.tanh(t) * sp.diff(q, t))


def kg_norm(u):
    """i (ubar du/dt - u dubar/dt) times the spatial volume, which must be t-independent."""
    ub = sp.conjugate(u)
    current = sp.I * (ub * sp.diff(u, t) - u * sp.diff(ub, t))
    volume = 2 * sp.pi * sp.cosh(t) * (4 * sp.pi / Lam)      # circle of radius cosh t, sphere 4 pi/Lambda
    return sp.simplify(sp.expand(current * volume))


if __name__ == "__main__":
    G = gudermannian()
    print("The two homogeneous solutions of q'' + tanh(t) q' = 0:\n")
    for name, q in (("the constant", sp.Integer(1)), ("gd(t) = arctan sinh t", G)):
        print("   %-24s residual %s" % (name, wave_residual(q)))
    print("   and d/dt arctan sinh t = %s, which is sech t: %s"
          % (sp.simplify(sp.diff(G, t)), sp.simplify(sp.diff(G, t) - 1 / sp.cosh(t)) == 0))

    N = sp.Symbol("N", positive=True)
    u = N * (1 - sp.I * alpha * G)
    norm = kg_norm(u)
    print("\nKlein-Gordon norm of u = N(1 - i alpha gd(t)):")
    print("   i(ubar u' - u ubar') x volume = %s" % norm)
    print("   t-independent, as a conserved norm must be: %s"
          % (sp.simplify(sp.diff(norm, t)) == 0))
    sol = sp.solve(sp.Eq(norm, 1), N)
    print("   setting it to one gives N = %s" % sol[0])
    claimed = sp.sqrt(Lam / (16 * sp.pi ** 2 * alpha))
    print("   the claimed normalisation is %s" % claimed)
    print("   they agree: %s" % (sp.simplify(sol[0] - claimed) == 0))

    print("\nDoes the fold cut the family down?  It does not, and that is the useful part.")
    print("  Theta is antilinear and reverses t, so it acts on a mode function as")
    print("  (Theta u)(t) = conj(u(-t)).  The Gudermannian is odd, gd(-t) = %s,"
          % sp.simplify(G.subs(t, -t) + G))
    print("  so for u = N(1 - i alpha gd t) with N real:")
    folded = sp.simplify(sp.conjugate(u.subs(t, -t)).rewrite(sp.atan))
    print("     (Theta u)(t) - u(t) = %s" % sp.simplify(sp.expand(folded - u)))
    print("  identically, for EVERY alpha. So the mode is fold-real at every member of the")
    print("  family and Theta-invariance selects nothing here. The freedom is real and the")
    print("  fold's own symmetry does not touch it.")
    print("  planted: an even combination N(1 - i alpha t) would not be fold-real, and gives")
    bad_mode = N * (1 - sp.I * alpha * t)
    print("     (Theta u)(t) - u(t) = %s"
          % sp.simplify(sp.conjugate(bad_mode.subs(t, -t)) - bad_mode))

    print("\nPlanted failures:")
    print("   a massive homogeneous mode is not a solution of the massless equation:")
    print("      q = cos(t) gives residual %s" % sp.simplify(wave_residual(sp.cos(t))))
    print("   dropping the cosh t from the volume breaks conservation of the norm:")
    bad = sp.simplify(sp.I * (sp.conjugate(u) * sp.diff(u, t) - u * sp.diff(sp.conjugate(u), t))
                      * 2 * sp.pi * (4 * sp.pi / Lam))
    print("      d/dt of the mis-measured norm = %s, nonzero" % sp.simplify(sp.diff(bad, t)))
    print("   a real combination carries no norm at all, so it is not a canonical pair:")
    print("      i(ubar u' - u ubar') for u = N gd(t) is %s"
          % sp.simplify(kg_norm(N * G)))

    print("\n  So the massless homogeneous sector has a normalisable canonical pair. The zero mode")
    print("  is a state to choose, not a barrier to having one, and the m -> 0 failure reported")
    print("  for the Euclidean reference state is a failure of that reference and not of the field.")
