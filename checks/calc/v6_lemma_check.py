#!/usr/bin/env python3
"""Separate Ways v6, section 4.1 and Appendix A: the operator inequality behind the ceiling, on the
four-dimensional Fock space of one fermion pair.

a_+, b_+ annihilate the two members of a pair (momenta p and -p) in the after-region; the
before-region operators are a_- = c a_+ + s b_+^dag, b_- = c b_+ - s a_+^dag, c^2 + s^2 = 1, s^2 = P.
N_+/- = (a^dag a + b^dag b)/2 in each region and Q = (N_+ + N_-)/2. Claims checked:
  1. the transformation is canonical (CAR hold for a_-, b_-);
  2. Q has eigenvalues {n*, 1/2, 1/2, 1 - n*}, n* = (1 - sqrt(1 - P))/2;
  3. every state with <N_+> = <N_-> (which Theta-invariance forces) has <N_+> >= n*, and the
     eigenvector of Q at n* has <N_+> = <N_-> = n* (saturation);
  4. plant: a state with <N_+> != <N_-> can go below n* (so the invariance is doing the work).
"""
import math, sys
import numpy as np
fails = 0
def check(ok, what):
    global fails
    print("   %-70s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)

sm = np.array([[0, 1], [0, 0]], float)       # lowering on one mode, basis (|0>, |1>)
Z = np.diag([1.0, -1.0]); I2 = np.eye(2)
a = np.kron(sm, I2); b = np.kron(Z, sm)      # Jordan-Wigner
dag = lambda A: A.conj().T
acomm = lambda A, B: A @ B + B @ A

rng = np.random.default_rng(7)
for P in (0.05, 0.2, 0.5, 0.8, 0.95):
    c, s = math.sqrt(1 - P), math.sqrt(P)
    am = c*a + s*dag(b); bm = c*b - s*dag(a)
    car = max(np.abs(acomm(am, dag(am)) - np.eye(4)).max(), np.abs(acomm(bm, dag(bm)) - np.eye(4)).max(),
              np.abs(acomm(am, bm)).max(), np.abs(acomm(am, dag(bm))).max())
    Np = (dag(a) @ a + dag(b) @ b)/2; Nm = (dag(am) @ am + dag(bm) @ bm)/2
    Q = (Np + Nm)/2
    formula = (1 - 2*P)*Np + P*np.eye(4) + c*s*(dag(a) @ dag(b) + b @ a)
    check(np.abs(Nm - formula).max() < 1e-14, "P = %.2f: N_- = (1-2P) N_+ + P + cs (a+^dag b+^dag + b+ a+)" % P)
    even = Q[np.ix_([0, 3], [0, 3])]
    check(abs(np.trace(even) - 1) < 1e-14 and abs(np.linalg.det(even) - P/4) < 1e-14, "P = %.2f: even block has trace 1, det P/4" % P)
    ev = np.sort(np.linalg.eigvalsh(Q)); ns = (1 - math.sqrt(1 - P))/2
    print("   P = %.2f: eig(Q) = %s, n* = %.6f" % (P, np.round(ev, 6), ns))
    check(car < 1e-14, "P = %.2f: a_-, b_- satisfy the anticommutation relations" % P)
    check(np.allclose(ev, [ns, 0.5, 0.5, 1 - ns], atol=1e-14), "P = %.2f: eigenvalues are {n*, 1/2, 1/2, 1 - n*}" % P)
    # random states, symmetrised so that <N+> = <N->: minimise <N+> along the constraint
    worst = 1.0
    w, V = np.linalg.eigh(Q); v0 = V[:, 0]
    check(abs(v0 @ Np @ v0 - ns) < 1e-12 and abs(v0 @ Nm @ v0 - ns) < 1e-12, "P = %.2f: Q's lowest eigenvector has <N+> = <N-> = n*" % P)
    for _ in range(4000):
        X = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4)); rho = X @ dag(X); rho /= np.trace(rho).real
        d = np.trace(rho @ (Np - Nm)).real
        # mix with the projector onto v0 shifted along N+ - N- to cancel the imbalance, keeping rho >= 0
        if abs(d) > 1e-12:
            # find a state with the opposite imbalance among eigenvectors of N+ - N-
            we, Ve = np.linalg.eigh(Np - Nm); u = Ve[:, 0] if d > 0 else Ve[:, -1]
            du = (np.conj(u) @ (Np - Nm) @ u).real
            lam = d/(d - du); rho = (1 - lam)*rho + lam*np.outer(u, np.conj(u))
        worst = min(worst, np.trace(rho @ Np).real)
    check(worst >= ns - 1e-12, "P = %.2f: 4000 states with <N+> = <N->: min <N+> = %.6f >= n*" % (P, worst))
    # plant: without the constraint <N+> can be smaller than n* (e.g. the after-region vacuum)
    vac = np.array([1, 0, 0, 0], complex)
    check((vac @ Np @ vac).real < ns, "plant P = %.2f: the after-vacuum, not invariant, sits below n*" % P)
# CPT maps helicity h to -h (Boyle, Finn, Turok, after their Eq. 96), so an invariant state relates the
# after-region of one helicity sector to the before-region of the other. Two sectors, four modes:
# S_+ = N_+^(h) + N_+^(-h) and S_- likewise; invariance gives <S_+> = <S_->, hence <S_+> = <Q_h + Q_-h> >= 2 n*.
print("\n=== two helicity sectors related by CPT ===")
I2b = np.eye(2)
def mode(k, n=4):
    """Jordan-Wigner lowering operator for mode k of n."""
    ops = [Z]*k + [sm] + [I2b]*(n - k - 1)
    out = ops[0]
    for o in ops[1:]: out = np.kron(out, o)
    return out
for P in (0.2, 0.5, 0.8):
    c, s_ = math.sqrt(1 - P), math.sqrt(P)
    ns = (1 - math.sqrt(1 - P))/2
    A1, B1, A2, B2 = (mode(k) for k in range(4))      # sector h: modes 0,1; sector -h: modes 2,3
    def sector(a, b):
        am = c*a + s_*dag(b); bm = c*b - s_*dag(a)
        return (dag(a) @ a + dag(b) @ b)/2, (dag(am) @ am + dag(bm) @ bm)/2
    Np1, Nm1 = sector(A1, B1); Np2, Nm2 = sector(A2, B2)
    Qsum = (Np1 + Nm1 + Np2 + Nm2)/2
    check(abs(np.linalg.eigvalsh(Qsum).min() - 2*ns) < 1e-12, "P = %.1f: lowest eigenvalue of Q_h + Q_-h is 2 n*" % P)
    # the counterexample: after-vacuum in sector h, before-vacuum in sector -h
    w1, V1 = np.linalg.eigh(Nm1 + 0*Np1); w2, V2 = np.linalg.eigh(Nm2)
    vac = np.zeros(16); vac[0] = 1.0                   # |0>: annihilated by A1, B1 (after-vacuum of sector h)
    # before-vacuum of sector -h: ground state of N_-^(-h) inside the after-vacuum of sector h
    M = Nm2 + 10.0*Np1                                 # penalise any sector-h occupation
    w, V = np.linalg.eigh(M); psi = V[:, 0]
    nh = psi @ Np1 @ psi; Sp = psi @ (Np1 + Np2) @ psi; Sm = psi @ (Nm1 + Nm2) @ psi
    print("   P = %.1f: <N+^h> = %.4f (n* = %.4f), <S+> = %.4f, <S-> = %.4f, 2 n* = %.4f" % (P, nh, ns, Sp, Sm, 2*ns))
    check(abs(Sp - Sm) < 1e-12, "P = %.1f: that state has <S+> = <S->, as invariance requires" % P)
    check(nh < ns - 1e-6, "plant P = %.1f: one sector alone falls below n* (the per-sector bound is false)" % P)
    check(Sp >= 2*ns - 1e-12, "P = %.1f: the bound summed over the two helicities holds" % P)
sys.exit(1 if fails else 0)
