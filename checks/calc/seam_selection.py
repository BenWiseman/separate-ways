#!/usr/bin/env python3
"""What the seam family does with energy and momentum, and the one member that is topological.

A.11 parametrises the seam by the corner term $S_c=(\\kappa/2)\\int dt\\,(q_1\\dot q_2-q_2\\dot q_1)$
and adopts the transparent member.  `seam_symmetry.py` shows no symmetry of the fold can pick a
member out, because the sheet exchange and the time reversal flip the term's sign in turn and
cancel.  This script asks what else the same action says, and the answer selects.

Everything is derived from the action here rather than assumed.  Two half-lines joined at the
origin, $\\phi_L$ on $x<0$ and $\\phi_R$ on $x>0$, with the bulk boundary terms and the corner
term varied together, give the matching conditions, hence the scattering matrix, hence the two
fluxes.  Energy is conserved at every $\\kappa$.  Momentum is not: the seam pushes, and the push
vanishes for arbitrary data at exactly one member of the family.

A defect that absorbs no momentum is one whose stress tensor is continuous across it, which is
what it means for an interface to be topological rather than material.  So the transparent seam is
the unique topological member of the family, and adopting it is adopting that premise rather than
adopting a number.  The premise does not follow from the fold and this script does not claim it
does.
"""
import sympy as sp

t, x, w, k = sp.symbols("t x omega kappa", real=True)
a1, a2 = sp.symbols("a1 a2")
q1d, q2d = sp.symbols("q1dot q2dot", real=True)


def matching():
    """Vary the two bulk actions and the corner term together at x = 0."""
    q1, q2 = sp.symbols("q1 q2")
    d1, d2 = sp.symbols("q1dot q2dot")
    pL, pR = sp.symbols("phiL_x phiR_x")
    # delta S_L gives -pL dq1, delta S_R gives +pR dq2 (opposite outward normals),
    # delta S_c gives kappa (q2dot dq1 - q1dot dq2).
    return sp.solve([sp.Eq(-pL + k * d2, 0), sp.Eq(pR - k * d1, 0)], [pL, pR], dict=True)[0]


def scattering():
    """b = S a for phi_L = a1 e^{iwx} + b1 e^{-iwx}, phi_R = a2 e^{-iwx} + b2 e^{iwx}."""
    b1, b2 = sp.symbols("b1 b2")
    qq1, qq2 = a1 + b1, a2 + b2
    lhsL = sp.I * w * (a1 - b1) - k * (-sp.I * w * qq2)      # phiL_x = kappa q2dot
    lhsR = sp.I * w * (b2 - a2) - k * (-sp.I * w * qq1)      # phiR_x = kappa q1dot
    sol = sp.solve([sp.Eq(lhsL, 0), sp.Eq(lhsR, 0)], [b1, b2], dict=True)[0]
    S = sp.Matrix(2, 2, lambda i, j: sp.simplify(sp.diff(sol[[b1, b2][i]], [a1, a2][j])))
    return sp.simplify(S)


if __name__ == "__main__":
    m = matching()
    print("Matching conditions from the action:")
    print("   phi_L,x(0) = %s      phi_R,x(0) = %s\n" % (m[sp.Symbol("phiL_x")], m[sp.Symbol("phiR_x")]))

    S = scattering()
    r = (1 - k ** 2) / (1 + k ** 2)
    tt = 2 * k / (1 + k ** 2)
    print("Scattering matrix, derived not assumed:")
    sp.pprint(S)
    print("   equals [[r, t], [-t, r]] with r = (1-k^2)/(1+k^2), t = 2k/(1+k^2): %s"
          % sp.simplify(S - sp.Matrix([[r, tt], [-tt, r]])).is_zero_matrix)
    print("   a ROTATION, det = %s, unitary: S^dagger S - 1 = %s"
          % (sp.simplify(S.det()), sp.simplify(S.H * S - sp.eye(2))))
    print("   NOT the symmetric [[r,t],[t,-r]] the naive argument assumes: difference %s\n"
          % sp.simplify(S - sp.Matrix([[r, tt], [tt, -r]])))

    E = sp.Matrix([[0, 1], [1, 0]])
    print("The fold law, exchange of sheets with exchange of in and out channels:")
    print("   E conj(S) E - S^dagger = %s   for every kappa, so no selection here."
          % sp.simplify(E * S.conjugate() * E - S.H))

    print("\nThe two fluxes at the seam, from the same matching conditions:")
    pL, pR = k * q2d, k * q1d
    energy = sp.simplify(-(q1d * pL) - (-(q2d * pR)))            # T_tx continuity
    momentum = sp.simplify(sp.Rational(1, 2) * ((q1d ** 2 + pL ** 2) - (q2d ** 2 + pR ** 2)))
    print("   energy defect   T_tx(0-) - T_tx(0+) = %s" % energy)
    print("   momentum defect T_xx(0-) - T_xx(0+) = %s" % sp.factor(momentum))
    print("   so energy is conserved at every kappa and momentum is not.")

    print("\n   The momentum defect vanishes for ARBITRARY data when its coefficient does:")
    sols = sp.solve(sp.Eq(1 - k ** 2, 0), k)
    print("      1 - kappa^2 = 0  ->  kappa = %s, and kappa = 1 is the transparent member."
          % sols)
    for kv in (sp.Rational(1, 2), 1, 2, 3):
        c = sp.simplify((1 - kv ** 2) / 2)
        print("      kappa = %-4s coefficient of (q1dot^2 - q2dot^2) is %s" % (kv, c))

    print("\n   And what that momentum IS, in the variables the cosmology paper splits fields into.")
    pc, pq = sp.symbols("Phic_dot Phiq_dot", real=True)
    subs = {q1d: pc + pq / 2, q2d: pc - pq / 2}
    rewritten = sp.simplify(sp.expand(momentum.subs(subs)))
    print("      with q1 = Phi_c + Phi_q/2 and q2 = Phi_c - Phi_q/2, the defect becomes")
    print("      %s" % sp.factor(rewritten))
    print("      which is (1 - kappa^2) times the product of the classical and quantum")
    print("      velocities. The seam's force is a classical-quantum cross term: it vanishes")
    print("      when the two sheets agree, when the classical field is static, and")
    print("      identically at the transparent point, and it has no long-range tail at all.")
    print("      planted: with q1 and q2 left alone the same expression is %s, which hides that."
          % sp.factor(momentum))

    print("\n   A defect absorbing no momentum has a continuous stress tensor across it, which is")
    print("   what makes an interface topological rather than material. So the transparent seam")
    print("   is the unique topological member of the family. That premise is not the fold's and")
    print("   is not derived here; what is derived is that it selects, and selects uniquely.")

    print("\nPlanted failures:")
    bad = sp.Matrix([[r, tt], [tt, -r]])
    print("   the symmetric matrix under the same fold law: E conj(S) E - S^dagger = %s"
          % sp.simplify(E * bad.conjugate() * E - bad.H))
    print("   a symmetric corner term, q1 q2dot + q2 q1dot, is a total derivative: d/dt(q1 q2)")
    self_coupled = sp.factor(sp.simplify(sp.Rational(1, 2) *
                             ((q1d ** 2 + (k * q1d) ** 2) - (q2d ** 2 + (k * q2d) ** 2))))
    print("   a self-coupled matching, phi_L,x = kappa q1dot, gives momentum defect %s,"
          % self_coupled)
    print("   whose coefficient is 1 + kappa^2 and never vanishes, so the selection is a")
    print("   property of the CROSS coupling and not of the bookkeeping.")
