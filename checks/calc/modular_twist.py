#!/usr/bin/env python3
"""Is the fold the modular conjugation?  A finite standard pair says not of the obvious state.

The fold is $\\Theta = J \\circ P_\\perp$: an antiunitary composed with the antipodal map.  It is
tempting to read that as Tomita-Takesaki, since an antiunitary carrying an algebra to its
commutant and fixing the vacuum is what a modular conjugation does, and to conclude that the
universe's horizon carries the structure Jacobson turns into the Einstein equation.  Three
properties get checked when people try this: the involution squares to one, it maps the algebra
onto its commutant, and it fixes the state.

Those three do not identify a modular conjugation, and the way they fail is exactly a parity
twist, which is exactly what the fold is.  The defining relation is not any of them.  It is

    S(A Omega) = A-dagger Omega,     S = J Delta^(1/2),

with the polar decomposition mattering: J antiunitary, Delta positive.

The counterexample here came from a second route and is checked independently.  Take the 2x2 matrices with the
Hilbert-Schmidt inner product, M acting by left multiplication, rho = diag(0.8, 0.2), and
Omega = q = sqrt(rho).  Write J(X) = X-dagger, Delta^(1/2)(X) = q X q^(-1), P = diag(1,-1), and
Theta(X) = P X-dagger P.

The repair came with it and it is the part worth having.  The same Theta IS the modular
conjugation of a different purification, Omega' = qP, which has the identical reduced state on
either factor.  So a parity-twisted fold is the modular conjugation of a parity-twisted
thermofield double, indistinguishable from the ordinary one on one sheet alone and differing only
in the cross-sheet correlations.  That is the silence theorem in operator-algebraic dress, and it
is a finite-dimensional statement and not a continuum construction.
"""
import numpy as np

np.set_printoptions(precision=10, suppress=True)
RHO = np.diag([0.8, 0.2])
Q = np.diag(np.sqrt([0.8, 0.2]))
P = np.diag([1.0, -1.0])
BASIS = [np.array(b, dtype=complex).reshape(2, 2)
         for b in ([1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1])]
NAMES = ["E11", "E12", "E21", "E22"]


def hs(X, Y):
    return np.trace(X.conj().T @ Y)


def theta(X, twist=P):
    return twist @ X.conj().T @ twist


def jmod(X):
    return X.conj().T


def delta_half(X):
    return Q @ X @ np.linalg.inv(Q)


def tomita_defect(conj, omega):
    """max_A || conj(Delta^(1/2)(A omega)) - A-dagger omega ||, over the matrix-unit basis."""
    worst, where = 0.0, None
    for name, A in zip(NAMES, BASIS):
        lhs = conj(delta_half(A @ omega))
        rhs = A.conj().T @ omega
        d = np.linalg.norm(lhs - rhs)
        if d > worst:
            worst, where = d, name
    return worst, where


def commutant_image(conj):
    """Does conj(L_A)conj land in the RIGHT multiplications, i.e. the full commutant?"""
    worst = 0.0
    for A in BASIS:
        B = P @ A.conj().T @ P                     # Theta L_A Theta = R_{P A-dagger P}
        for X in BASIS:
            lhs = conj(A @ conj(X))
            worst = max(worst, np.linalg.norm(lhs - X @ B))
    return worst


if __name__ == "__main__":
    print("Standard pair: 2x2 matrices, Hilbert-Schmidt, rho = diag(0.8, 0.2), Omega = sqrt(rho)\n")

    print("The three properties usually offered, for Theta(X) = P X-dagger P:")
    sq = max(np.linalg.norm(theta(theta(X)) - X) for X in BASIS)
    anti = max(abs(hs(theta(X), theta(Y)) - hs(Y, X)) for X in BASIS for Y in BASIS)
    print("   Theta^2 = 1                       residual %.3e" % sq)
    print("   antiunitary, <ThX,ThY> = <Y,X>    residual %.3e" % anti)
    print("   Theta Omega = Omega               residual %.3e" % np.linalg.norm(theta(Q) - Q))
    print("   Theta M Theta lands in M'         residual %.3e" % commutant_image(theta))
    print("   all four pass.\n")

    print("The defining relation, which is the one that decides it:")
    dth, whth = tomita_defect(theta, Q)
    dj, whj = tomita_defect(jmod, Q)
    print("   canonical J    defect %.3e   (worst at %s)" % (dj, whj))
    print("   twisted Theta  defect %.12f   (worst at %s)" % (dth, whth))
    print("   and 2*sqrt(0.8) = %.12f, which is what that number is." % (2 * np.sqrt(0.8)))

    ts = np.zeros((4, 4), dtype=complex)
    for j, X in enumerate(BASIS):
        img = P @ delta_half(X) @ P                 # Theta S = U Delta^(1/2), since J J = 1
        for i, Y in enumerate(BASIS):
            ts[i, j] = hs(Y, img)
    ev = np.sort(np.linalg.eigvals(ts).real)
    print("   no positive polar factor rescues it: Theta S has eigenvalues %s" % ev)
    print("   a positive operator cannot have a negative eigenvalue, so Theta is not a J.\n")

    print("The repair: change the purification, not the conjugation.")
    omega2 = Q @ P
    rho2 = omega2 @ omega2.conj().T
    print("   Omega' = qP has reduced state diag %s, identical to rho" % np.diag(rho2).real)
    d2, wh2 = tomita_defect(theta, omega2)
    print("   Theta's defining-relation defect against Omega' is %.3e (worst at %s)" % (d2, wh2))
    print("   and the canonical J's defect against Omega' is %.6f, so they have swapped roles."
          % tomita_defect(jmod, omega2)[0])
    print("   overlap <Omega, Omega'> = %.4f, so it is a different global state," % hs(Q, omega2).real)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    cross = lambda w: (np.trace(w.conj().T @ sx @ w @ sx) / 2).real   # <Omega| L_sx R_sx |Omega>
    print("   and the cross-leg correlation <L_sx R_sx> moves from %+.2f to %+.2f"
          % (cross(Q), cross(omega2)))
    print("   while the reduced state on either sheet does not move at all.")

    print("\nPlanted failures, so none of the above can pass by matching nothing:")
    for label, twist in (("no twist, P -> 1", np.eye(2)),
                         ("wrong twist, P -> diag(1, 2)", np.diag([1.0, 2.0]))):
        t = lambda X, w=twist: w @ X.conj().T @ w
        print("   %-28s Theta^2 residual %.3e, Tomita defect vs Omega %.4f, vs Omega' %.4f"
              % (label, max(np.linalg.norm(t(t(X)) - X) for X in BASIS),
                 tomita_defect(t, Q)[0], tomita_defect(t, omega2)[0]))
    print("   the untwisted case passes against Omega and fails against Omega', which is the")
    print("   whole content: the conjugation and the purification have to be chosen together.")
