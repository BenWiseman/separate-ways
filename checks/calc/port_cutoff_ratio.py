#!/usr/bin/env python3
"""The one-dimensional port model's C(a), and why it is not a property of the geometry.

A.11 rejects the one-dimensional two-port model as a stand-in for A.6's object, and the
evidence it gives is a table: the ratio C falls steadily as the cutoff is lifted, so it
belongs to the cutoff rather than to the port separation.  The table had no script in the
release until 2026-09-24, found by checking every number of three or more decimals in the
manuscript against CLAIMS.tsv.  The numbers are right; nothing reproduced them.

The quantity, as A.11 writes it:

    C(a) = int dk cos(2 k a) / (2 omega_k)   /   int dk 1 / (2 omega_k),

with omega_k = sqrt(k^2 + m^2) and m = 1, integrated from zero to a cutoff k_max.  The
numerator converges, to K_0(2a)/2 in the limit, while the denominator is arcsinh(k_max)/2
and grows without bound.  So C falls like 1/log(k_max) and has no limit that could be a
property of the ports.  That is the whole of A.11's point, and the table is how it shows.
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import k0

A_SEP = 0.02
CUTOFFS = (20, 100, 400, 1000, 4000, 20000)
QUOTED = (0.958, 0.603, 0.497, 0.441, 0.371, 0.315)


def integrals(a, kmax, m=1.0, power=1.0):
    """Numerator and denominator.  power != 1 deforms omega_k and is only for the plant."""
    w = lambda k: 2.0 * (k * k + m * m) ** (power / 2.0)
    num = quad(lambda k: np.cos(2 * k * a) / w(k), 0, kmax, limit=6000)[0]
    den = quad(lambda k: 1.0 / w(k), 0, kmax, limit=6000)[0]
    return num, den


def ratio(a, kmax, m=1.0, power=1.0):
    num, den = integrals(a, kmax, m, power)
    return num / den


if __name__ == "__main__":
    print(__doc__)
    print("=" * 78)

    print(f"\n[1] the table A.11 quotes, at a = {A_SEP}")
    worst = 0.0
    for kmax, quoted in zip(CUTOFFS, QUOTED):
        got = ratio(A_SEP, kmax)
        worst = max(worst, abs(got - quoted))
        print(f"    k_max {kmax:>6}:  C = {got:.4f}   the paper says {quoted:.3f}")
    print(f"    worst disagreement {worst:.4f}, which is the rounding of 0.4415 to 0.441")
    assert worst < 1e-3, "the table does not reproduce at the precision it is quoted to"

    print("\n[2] the denominator diverges logarithmically, as claimed")
    for kmax in (1e2, 1e4, 1e6):
        _, den = integrals(A_SEP, kmax)
        print(f"    k_max {kmax:.0e}: denominator {den:.6f}, "
              f"arcsinh(k_max)/2 = {np.arcsinh(kmax) / 2:.6f}")
        assert abs(den - np.arcsinh(kmax) / 2) < 1e-6

    print("\n[3] the numerator converges, to K_0(2a)/2")
    for kmax in (1e3, 1e4, 1e5):
        num, _ = integrals(A_SEP, kmax)
        print(f"    k_max {kmax:.0e}: numerator {num:.8f}, "
              f"K_0(2a)/2 = {k0(2 * A_SEP) / 2:.8f}")
    assert abs(integrals(A_SEP, 1e5)[0] - k0(2 * A_SEP) / 2) < 1e-3
    print("    stopped at 1e5 on purpose: beyond it the adaptive quadrature loses the")
    print("    oscillation and returns nonsense, which is the integrator and not the physics")

    print("\n[4] so C falls like 1/log k_max and has no limit")
    prev = None
    for kmax in (1e2, 1e3, 1e4, 1e5, 1e6):
        c = ratio(A_SEP, kmax)
        print(f"    k_max {kmax:.0e}: C = {c:.5f}   C log(k_max) = {c * np.log(kmax):.5f}")
        assert prev is None or c < prev, "C must fall"
        prev = c

    print("\n[5] plant: give the denominator a convergent measure and the fall must stop")
    print("    omega_k^2 in place of omega_k, so both integrals converge")
    ks = (1e2, 1e3, 1e4)      # same range the real case is trusted over, for the same reason
    vals = [ratio(A_SEP, k, power=2.0) for k in ks]
    for k, v in zip(ks, vals):
        print(f"    k_max {k:.0e}: C = {v:.6f}")
    spread = (max(vals) - min(vals)) / max(vals)
    real = [ratio(A_SEP, k) for k in ks]
    print(f"    relative spread {spread:.2e}, against "
          f"{(max(real) - min(real)) / max(real):.2f} for the real measure over the same range")
    assert spread < 1e-2, "the plant still falls, so the check cannot tell the two apart"
