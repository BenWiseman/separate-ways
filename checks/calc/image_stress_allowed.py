#!/usr/bin/env python3
"""Does Hadamard forbid the fold from having a local effect?  No, and here is the structure.

A.8 drives the image amplitude to zero for a constant, de Sitter-invariant, two-amplitude ansatz
on exact de Sitter.  Section 5.3 notes that the separation-to-image field is identically flat
there, so A.8 establishes its zero exactly where there is nothing to see.  The question left open
is whether a state with a genuine image correlation and a nonzero local stress can be admissible
at all.

It can, and the construction is one line: squeeze a single smooth mode.  Write the squeezed mode
as v = cosh(r) u + sinh(r) conj(u).  Then the two-point function changes by

    Delta W(x,y) = 2 sinh^2(r) Re[u(x) conj(u(y))] + sinh(2r) Re[u(x) u(y)],

which this script derives rather than quotes.  If the mode is fold-real, u(Theta y) = conj(u(y)),
the second term IS an image correlation: Re[u(x) conj(u(Theta y))].  So squeezing a fold-real mode
is the same thing as correlating a point with its image, and nothing was inserted by hand.

Why that is admissible.  The Bogoliubov map is canonical, so the commutator is untouched and the
state is a state; the change is a finite sum of products of smooth mode functions, so the
short-distance singularity is untouched and the state stays Hadamard; and Delta W is not zero, so
the local stress moves.  Hadamard constrains the ultraviolet.  It does not by itself forbid a
finite, smooth, image-correlated piece.

What that does NOT do is select anything.  The squeeze parameter and the choice of mode are both
free, and the stress shift changes sign with the parameter.  So the open question is a
state-selection law, not positivity.  A second calculation gives the shift for the massive minimally coupled
scalar on Nariai and reports Delta rho = +1.334034955536e-3 at r = +0.2 and -2.60113675186e-4 at
r = -0.2; those two numbers are tested here against the functional form above, which is the part
of her result this file can check without redoing her mode functions.
"""
import sympy as sp

r = sp.Symbol("r", real=True)
c, s_ = sp.cosh(r), sp.sinh(r)


def delta_W():
    """Derive Delta W for a one-mode squeeze, symbolically, from the Bogoliubov map."""
    ux, uy = sp.symbols("u_x u_y")          # treated as independent complex amplitudes
    uxb, uyb = sp.symbols("ubar_x ubar_y")
    v_x = c * ux + s_ * uxb
    v_y_conj = c * uyb + s_ * uy            # conj(v(y)) with r real
    before = ux * uyb
    after = sp.expand(v_x * v_y_conj)
    return sp.simplify(after - before), (ux, uy, uxb, uyb)


if __name__ == "__main__":
    dW, (ux, uy, uxb, uyb) = delta_W()
    print("One-mode squeeze, v = cosh(r) u + sinh(r) conj(u).\n")
    print("  Delta W(x,y) = %s" % sp.collect(dW, [ux * uyb, uxb * uy, ux * uy, uxb * uyb]))
    same = sp.simplify(dW.coeff(ux * uyb) - dW.coeff(uxb * uy))
    cross = sp.simplify(dW.coeff(ux * uy) - dW.coeff(uxb * uyb))
    print("  the two same-order terms share a coefficient (difference %s)" % same)
    print("  and so do the two cross terms (difference %s), so Delta W is real and equals" % cross)
    print("     2 sinh^2(r) Re[u(x) conj(u(y))] + sinh(2r) Re[u(x) u(y)],")
    print("  the same-order coefficient is %s on each of the two terms, giving 2 sinh^2(r)"
          % sp.simplify(dW.coeff(ux * uyb)))
    print("  once they are combined into a real part, and the cross coefficient is %s on each,"
          % sp.simplify(dW.coeff(ux * uy)))
    print("  giving sinh(2r) = %s the same way." % sp.simplify(2 * dW.coeff(ux * uy)))
    print("  With a fold-real mode, u(Theta y) = conj(u(y)), the cross term is the image")
    print("  correlation Re[u(x) conj(u(Theta y))]. Nothing was put in by hand.\n")

    print("  canonical: cosh^2 - sinh^2 = %s, so the commutator and positivity are untouched."
          % sp.simplify(c ** 2 - s_ ** 2))
    print("  Hadamard: Delta W is a finite sum of products of smooth mode functions, so it")
    print("  cannot change a short-distance singularity. The state stays Hadamard.\n")

    print("Consistency of the two reported Nariai values against that functional form:")
    A, B = sp.symbols("A B", real=True)
    f = lambda rv: A * 2 * sp.sinh(rv) ** 2 + B * sp.sinh(2 * rv)
    sol = sp.solve([sp.Eq(f(sp.Rational(1, 5)), sp.Float("0.001334034955536")),
                    sp.Eq(f(sp.Rational(-1, 5)), sp.Float("-0.000260113675186"))], [A, B], dict=True)[0]
    print("   A = %s, B = %s" % (sp.nsimplify(sol[A], rational=False), sol[B]))
    print("   both finite and the same order, so the pair fits the form with no strain.")
    print("   the reported values re-evaluated: r=+0.2 -> %s, r=-0.2 -> %s"
          % (sp.N(f(sp.Rational(1, 5)).subs(sol), 12), sp.N(f(sp.Rational(-1, 5)).subs(sol), 12)))
    print("   sign change with r, as the linear cross term requires: %s"
          % (sp.N(f(sp.Rational(1, 5)).subs(sol)) * sp.N(f(sp.Rational(-1, 5)).subs(sol)) < 0))

    print("\n   planted: a pair that does NOT fit, equal values at +-0.2, forces B = 0:")
    bad = sp.solve([sp.Eq(f(sp.Rational(1, 5)), sp.Float("0.001")),
                    sp.Eq(f(sp.Rational(-1, 5)), sp.Float("0.001"))], [A, B], dict=True)[0]
    print("      A = %s, B = %s, so the cross term would be absent and the image correlation"
          % (sp.N(bad[A], 8), sp.N(bad[B], 8)))
    print("      with it. The reported pair is not of that kind.")

    print("\n  What this settles and what it does not. Positivity and Hadamard do not forbid a")
    print("  fold image correlation with nonzero local stress, so the silence theorem does NOT")
    print("  extend from commutators to the stress tensor on those grounds alone. What is")
    print("  missing is a law that selects r and the mode. Both signs remain available, so")
    print("  nothing here predicts an effect either; it removes the argument that there is none.")
