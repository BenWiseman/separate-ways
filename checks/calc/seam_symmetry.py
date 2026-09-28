#!/usr/bin/env python3
"""No symmetry of the fold can fix the seam coefficient, and here is why.

The seam is a corner term, $(\\kappa/2)\\int dt\\,(q_1\\dot q_2 - q_2\\dot q_1)$, with $q_1$ and
$q_2$ the field values on the two sheets.  A.11 adopts the transparent member $\\kappa=1$ and says
what deriving it would take.  The first thing anyone tries, and the first thing tried here, is to
demand that the fold itself be a symmetry of the seam, on the grounds that the seam IS the fold
locus and cannot be less symmetric than the map that defines it.

That demand is satisfied at every $\\kappa$, so it fixes nothing, and the reason is one line.  The
fold does two things: it exchanges the sheets, and it reverses time orientation.  The corner term
is antisymmetric in the sheet labels, so the exchange flips its sign.  It is first order in a time
derivative, so reversing time flips its sign again.  The two cancel exactly, for any coefficient.

This is checked below on a concrete pair of profiles rather than argued, with each half of the
fold applied separately as the planted failure: either one alone must flip the sign, and both
together must leave the term alone.  A two-by-two scattering ansatz gives the opposite answer,
that the fold forces zero reflection, and it is the ansatz that is wrong: it tracks the sheet
exchange and drops the exchange of in and out channels that comes with time reversal.

What follows for the paper.  The adoption of $\\kappa=1$ is not a gap in the symmetry analysis
waiting to be filled.  No symmetry of the construction can fix the coefficient, so it has to come
from dynamics or from a completion, which is what A.11 says and is now a statement with a
calculation behind it rather than a report of two failed attempts.
"""
import mpmath as mpm

mpm.mp.dps = 30

PROFILES = [
    ("Gaussian times sine and cosine",
     lambda t: mpm.e ** (-t ** 2) * mpm.sin(3 * t),
     lambda t: mpm.e ** (-t ** 2 / 2) * mpm.cos(2 * t + 1)),
    ("a pair with no parity of its own",
     lambda t: mpm.e ** (-(t - mpm.mpf("0.7")) ** 2) * (1 + t),
     lambda t: mpm.e ** (-(t + mpm.mpf("0.3")) ** 2 / 3) * mpm.sin(t ** 2)),
]


def corner(a, b, span=8):
    """The seam's corner term, up to the coefficient kappa/2 which cancels in every ratio."""
    return mpm.quad(lambda t: a(t) * mpm.diff(b, t) - b(t) * mpm.diff(a, t), [-span, 0, span])


if __name__ == "__main__":
    print("Corner term (q1 q2' - q2 q1') under the fold and under each half of it.\n")
    print("  %-34s %-16s %-16s %s" % ("", "value", "vs original", "verdict"))
    for label, q1, q2 in PROFILES:
        base = corner(q1, q2)
        cases = [
            ("fold: swap and reverse time", lambda t: q2(-t), lambda t: q1(-t), 1),
            ("planted: swap only", q2, q1, -1),
            ("planted: time reversal only", lambda t: q1(-t), lambda t: q2(-t), -1),
        ]
        print("  %s" % label)
        print("  %-34s %-16s" % ("    original", mpm.nstr(base, 12)))
        for name, a, b, want in cases:
            val = corner(a, b)
            ratio = val / base
            ok = abs(ratio - want) < mpm.mpf("1e-20")
            print("  %-34s %-16s %-16s %s"
                  % ("    " + name, mpm.nstr(val, 12), mpm.nstr(ratio, 8),
                     "as required" if ok else "UNEXPECTED"))
        print()

    print("  So the seam action is fold-invariant for every kappa: the sheet exchange flips the")
    print("  sign, the time reversal flips it back, and the coefficient never enters. A symmetry")
    print("  that holds for every member of a family cannot select one.")
