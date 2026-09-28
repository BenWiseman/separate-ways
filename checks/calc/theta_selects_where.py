#!/usr/bin/env python3
"""When does fold-invariance pick the state, and when does it pick nothing?

Item 14 asked what selects the amplitude of the image correlation, because the squeeze parameter
is free and the stress shift changes sign with it.  The programme already has a selection, at the
bang: Theta-invariance forces the contact condition and with it a definite occupation.  So the
question is not whether the fold can select.  It is why it selects there and not at a horizon.

The answer is one line and it is checked here rather than argued.  For an antiunitary Theta with
Theta^2 = 1, a Theta-invariant state and a Hermitian A,

    <Theta A Theta> = conj(<A>) = <A>.

If Theta EXCHANGES two things, put A = N_1 and the two occupations are forced equal.  That is the
contact condition, and it is a consequence of antilinearity and hermiticity rather than a
prescription.  If Theta FIXES the thing, as it does for a fold-real mode at a horizon, there is no
second operator to equate and the same identity constrains nothing.

Both halves are exhibited below on explicit finite states, with the failure controls that matter:
a state that is not Theta-invariant must break the first, and the second must remain unselected
across a range of squeezes.
"""
import numpy as np

np.set_printoptions(precision=12, suppress=True)
NMAX = 5                                            # Fock cutoff per mode


def two_mode_ops():
    a = np.diag(np.sqrt(np.arange(1, NMAX)), 1)
    I = np.eye(NMAX)
    return np.kron(a, I), np.kron(I, a)


def swap_matrix():
    S = np.zeros((NMAX * NMAX, NMAX * NMAX))
    for i in range(NMAX):
        for j in range(NMAX):
            S[j * NMAX + i, i * NMAX + j] = 1.0
    return S


def theta_conj(rho, S):
    """Theta = S . K acting on a density matrix: Theta rho Theta = S conj(rho) S."""
    return S @ rho.conj() @ S


def random_state(seed, dim):
    rng = np.random.default_rng(seed)
    M = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
    rho = M @ M.conj().T
    return rho / np.trace(rho).real


def squeezed_vacuum(r):
    """|r> ~ sum_n tanh(r)^n sqrt((2n-1)!!/(2n)!!) |2n>, real amplitudes for real r."""
    amps = np.zeros(NMAX)
    c, t = 1.0 / np.cosh(r), np.tanh(r)
    coef = 1.0
    for n in range(0, NMAX, 2):
        k = n // 2
        if k:
            coef *= np.sqrt((2 * k - 1) / (2 * k))
        amps[n] = np.sqrt(c) * coef * t ** k
    return amps / np.linalg.norm(amps)


if __name__ == "__main__":
    a1, a2 = two_mode_ops()
    N1, N2 = a1.conj().T @ a1, a2.conj().T @ a2
    S = swap_matrix()
    dim = NMAX * NMAX

    print("1. Theta EXCHANGES the two sheets: the contact condition is forced.\n")
    print("   state                         <N1>            <N2>            difference")
    for seed in (1, 7, 23):
        rho = random_state(seed, dim)
        inv = (rho + theta_conj(rho, S)) / 2                 # project onto Theta-invariant
        inv = inv / np.trace(inv).real
        d = abs(np.trace(inv @ N1).real - np.trace(inv @ N2).real)
        print("   Theta-invariant, seed %-3d     %-15.10f %-15.10f %.2e"
              % (seed, np.trace(inv @ N1).real, np.trace(inv @ N2).real, d))
    print("\n   planted: the SAME states before projecting, which must not obey it:")
    for seed in (1, 7, 23):
        rho = random_state(seed, dim)
        d = abs(np.trace(rho @ N1).real - np.trace(rho @ N2).real)
        print("   not invariant, seed %-3d       %-15.10f %-15.10f %.3f"
              % (seed, np.trace(rho @ N1).real, np.trace(rho @ N2).real, d))
    print("\n   The identity behind it: for antiunitary Theta with Theta^2 = 1, a Theta-invariant")
    print("   state and Hermitian A, <Theta A Theta> = conj(<A>) = <A>. Put A = N1 and the two")
    print("   occupations are equal. Antilinearity and hermiticity, nothing else.\n")

    print("2. Theta FIXES a single fold-real mode: nothing is selected.\n")
    print("   A squeezed vacuum with real r has real Fock amplitudes, so Theta = K leaves it")
    print("   invariant at every r, while the occupation runs freely:\n")
    print("   r        Theta-invariance defect     <N> = sinh^2 r")
    for r in (0.0, 0.2, 0.5, 1.0):
        psi = squeezed_vacuum(r)
        defect = np.abs(psi - psi.conj()).max()
        n = float(psi @ (np.diag(np.arange(NMAX))) @ psi)
        print("   %-8.2f %-27.2e %.10f" % (r, defect, n))
    print("\n   planted: give the same state a complex phase per level and Theta-invariance")
    print("   does bite, so the check is not vacuous:")
    for r in (0.2, 0.5):
        psi = squeezed_vacuum(r).astype(complex)
        psi = psi * np.exp(1j * 0.7 * np.arange(NMAX))
        print("   r = %-5.2f complex phases   defect %.4f" % (r, np.abs(psi - psi.conj()).max()))

    print("\n3. And the first case bounds the second, which is what makes this more than a dichotomy.\n")
    print("   At the bang Theta-invariance forces half-and-half occupancy, so the pair-block floor")
    print("   is n_* = (1 - sqrt(1-P))/2 with P the crossing probability. A squeeze carries")
    print("   <N> = sinh^2 r, so the selected occupation IS a selected squeeze:\n")
    import mpmath as mpm
    mpm.mp.dps = 20
    print("      P        n_*            r = arcsinh(sqrt(n_*))")
    for Pv in ("0.05", "0.2", "0.5", "0.9", "1.0"):
        Pv = mpm.mpf(Pv)
        n_ = (1 - mpm.sqrt(1 - Pv)) / 2
        print("      %-8s %-14s %s" % (mpm.nstr(Pv, 3), mpm.nstr(n_, 8),
                                       mpm.nstr(mpm.asinh(mpm.sqrt(n_)), 10)))
    sup = mpm.asinh(1 / mpm.sqrt(2))
    print("\n   n_* is bounded by 1/2 at P = 1, so r is bounded by arcsinh(1/sqrt 2) = %s."
          % mpm.nstr(sup, 12))
    print("   planted: n_* > 1/2 would be needed to exceed it, and n_* = (1-sqrt(1-P))/2 reaches")
    print("      1/2 only at P = 1; at P = 1.2 the expression is complex, %s, so there is no"
          % mpm.nstr((1 - mpm.sqrt(1 - mpm.mpf("1.2"))) / 2, 6))
    print("      real crossing that overshoots.")
    print("\n   Two conditions on reading that as a bound at a horizon, and both are assumptions:")
    print("   that the bang's invariance applies to the mode in question, and that its squeeze is")
    print("   inherited rather than regenerated on the way. Neither is established here. What is")
    print("   established is that the selected quantity and the free quantity are the same")
    print("   quantity, so the freedom was never unbounded in principle.")

    print("\n  So the fold selects where it exchanges and not where it fixes. The bang is the")
    print("  first case and a fold-real mode at a horizon is the second. The horizon state is")
    print("  therefore not a free choice to be justified by a new principle: it is whatever the")
    print("  selected bang state evolves into, which is a propagation problem and not a")
    print("  foundational one.")
