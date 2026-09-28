#!/usr/bin/env python3
"""The propagation problem, solved: what squeeze the fold carries at a horizon.

A.8 leaves the image correlation permitted and its amplitude unfixed, and records that the
amplitude is not a foundational gap but a propagation question: the state at a horizon is whatever
the selected bang state evolves into.  This file does that propagation, and the answer has no free
parameter in it.

The step that makes it easy is one the paper already owns.  Section 2 establishes that the fold's
map is the SQUARE ROOT of the thermal transformation at the primitive period.  A state in
equilibrium at a horizon is KMS at the Hawking temperature, so its cross-sheet correlator is the
direct one shifted by half a period, and that is not a choice but a consequence.  Carry it through
for one mode of frequency omega at inverse temperature beta and the half-period shift turns the
Wightman function into

    W(t - i beta/2) = cos(omega t) / (2 omega sinh(beta omega / 2)),

real and even, which is the fold-real property A.8 needs.  A two-mode state with that cross
correlator is the thermofield double, and the thermofield double is a two-mode squeezed vacuum
with

    tanh r(omega) = exp(-beta omega / 2),       <N> = sinh^2 r = 1/(e^{beta omega} - 1),

so the squeeze parameter at a horizon is fixed by the mode frequency and the surface gravity
alone.  Everything here is checked below rather than quoted, including the identity that makes the
two exponentials of the shifted correlator collapse to one coefficient.

Two caveats travel with it.  The fold is J composed with the antipodal map, not J, so the fold's
state is the thermofield double with an extra transverse parity: an independent check finds
that the partner leaves the marginals unchanged while the cross-sheet term picks up (-1)^l per
multipole, which is why a stress evaluated at one r came out with opposite signs at the pole and
at the equator.  And equilibrium is an assumption about the state at the horizon, not a theorem.
"""
import sympy as sp

w, beta, t = sp.symbols("omega beta t", positive=True)


def thermal_wightman(shift=0):
    """W(t) for one mode at inverse temperature beta, optionally shifted by -i beta/2."""
    n = 1 / (sp.exp(beta * w) - 1)
    tt = t - sp.I * beta * shift
    return sp.simplify(((1 + n) * sp.exp(-sp.I * w * tt) + n * sp.exp(sp.I * w * tt)) / (2 * w))


if __name__ == "__main__":
    print("1. The half-period shift, which the fold IS.\n")
    W0 = thermal_wightman(0)
    Wh = sp.simplify(sp.expand(thermal_wightman(sp.Rational(1, 2))))
    target = sp.cos(w * t) / (2 * w * sp.sinh(beta * w / 2))
    print("   W(t)            = %s" % sp.simplify(W0))
    print("   W(t - i beta/2) = %s" % sp.simplify(Wh))
    print("   claimed closed form cos(wt)/(2 w sinh(beta w/2)); difference = %s"
          % sp.simplify(sp.expand_complex(sp.simplify(Wh - target))))
    print("   the two exponentials collapse because (1+n)e^{-bw/2} and n e^{+bw/2} are equal:")
    n = 1 / (sp.exp(beta * w) - 1)
    print("      (1+n) e^{-bw/2} - n e^{+bw/2} = %s"
          % sp.simplify((1 + n) * sp.exp(-beta * w / 2) - n * sp.exp(beta * w / 2)))
    coeff = sp.simplify(sp.nsimplify(sp.simplify(
        (1 + n) * sp.exp(-beta * w / 2) - 1 / (2 * sp.sinh(beta * w / 2))).rewrite(sp.exp)))
    print("      and each equals 1/(2 sinh(bw/2)): difference %s" % sp.simplify(sp.expand(coeff)))
    print("   So the fold-shifted correlator is real and even, which is the fold-real property.\n")

    print("2. The state with that cross correlator is a two-mode squeeze, and its parameter.\n")
    r = sp.Symbol("r", positive=True)
    tanh_r = sp.exp(-beta * w / 2)
    occ = sp.simplify(sp.sinh(sp.atanh(tanh_r)) ** 2)
    print("   tanh r = e^{-beta omega/2}  =>  <N> = sinh^2 r = %s" % occ)
    print("   which must be the Bose factor 1/(e^{beta omega} - 1); difference = %s"
          % sp.simplify(occ - 1 / (sp.exp(beta * w) - 1)))
    print("   planted: tanh r = e^{-beta omega} instead would give <N> = %s,"
          % sp.simplify(sp.sinh(sp.atanh(sp.exp(-beta * w))) ** 2))
    print("   which is not the Bose factor: difference = %s"
          % sp.simplify(sp.sinh(sp.atanh(sp.exp(-beta * w))) ** 2 - 1 / (sp.exp(beta * w) - 1)))

    print("\n3. Numbers, in units of the horizon temperature.\n")
    import mpmath as mpm
    mpm.mp.dps = 20
    print("   beta*omega   r = arctanh(e^{-bw/2})   <N> = sinh^2 r")
    for bw in ("0.2", "0.5", "1.0", "ln3", "2.0", "5.0", "10.0"):
        bw = mpm.log(3) if bw == "ln3" else mpm.mpf(bw)
        rv = mpm.atanh(mpm.e ** (-bw / 2))
        print("   %-12s %-24s %s" % (mpm.nstr(bw, 5), mpm.nstr(rv, 10),
                                     mpm.nstr(mpm.sinh(rv) ** 2, 10)))
    bang = mpm.asinh(1 / mpm.sqrt(2))
    print("\n   The bang's ceiling was r <= arcsinh(1/sqrt2) = %s, which is <N> = 1/2."
          % mpm.nstr(bang, 10))
    crossing = mpm.log(3)
    print("   A horizon mode reaches <N> = 1/2 at beta*omega = ln 3 = %s, and exceeds it below."
          % mpm.nstr(crossing, 10))
    r_at_ln3 = mpm.atanh(mpm.e ** (-crossing / 2))
    print("   check: at beta*omega = ln 3, <N> = %s and r = %s"
          % (mpm.nstr(mpm.sinh(r_at_ln3) ** 2, 12), mpm.nstr(r_at_ln3, 12)))
    print("   which is the bang's ceiling to %s: the two are the same number, because"
          % mpm.nstr(abs(r_at_ln3 - bang), 3))
    print("   sinh r = 1/sqrt2 gives cosh r = sqrt(3/2) and so tanh r = 1/sqrt3 = e^{-ln3/2}.")
    print("   So the bang's ceiling is not a separate fact: it is the horizon's own value at")
    print("   omega = T ln 3. Softer modes sit above it, harder ones below.")
    print("\n   So the bang's bound does NOT propagate as a ceiling: soft modes are repopulated")
    print("   thermally and sit above it. What replaces the bound is better than a bound, because")
    print("   the horizon squeeze is not bounded but DETERMINED, by beta and omega and nothing")
    print("   else. The free amplitude is gone.")
