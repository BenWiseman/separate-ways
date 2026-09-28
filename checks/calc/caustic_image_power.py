#!/usr/bin/env python3
"""How strong the image singularity is when the image point sits on a caustic.

Section 5.3 withdrew a power.  The reading there was that the fold censors its own closed
causal curves, because the world function between a point and its fold image vanishes
linearly at the contact sphere, so the image term of the two-point function blows up on it.
Counting that blow-up through the Hadamard parametrix, image term ~ Delta^(1/2)/sigma and a
stress tensor quadratic in derivatives, gives s^-3.  The count needs the Van Vleck factor
finite at the image point, and it is not: the fold's transverse factor is the antipodal map
on a round S^2, every great circle through a point refocuses at its antipode, so the image
point does not merely sit near a caustic, it sits on one.  Delta diverges there and the
parametrix is being expanded about the one place it does not hold.

A parametrix is an approximation.  The exact two-point function is not, and a mode sum needs
no Van Vleck factor at all.  This file computes the singularity in a model that isolates the
one feature at issue, and gets the power the parametrix could not.

MODEL.  ds^2 = -dt^2 + dr^2 + a^2 dOmega^2, a free massless scalar, a = 1.  This is not
Schwarzschild and is not claimed to be.  It is the smallest geometry carrying the feature:
null geodesics leaving a point spread over every great circle of the transverse sphere and
refocus, all of them together, after exactly pi of transverse angle.  The connecting family
at the image point is one-parameter, which is the same degeneracy the fold's contact sphere
carries, and the question of what a two-point function does at such a pair is local to the
caustic and not to the hole.

Reducing on the sphere turns the 4D field into a tower of 2D fields of mass m_l = sqrt(l(l+1))/a:

    W(x,x') = sum_l (2l+1)/(4 pi) P_l(cos gamma) W_2(m_l; Delta t, Delta r).

The fold's transverse factor sends gamma -> pi - gamma, which by P_l(-x) = (-1)^l P_l(x) is
exactly the insertion of (-1)^l that A.8 records on the cross-sheet term.  The two statements
are the same statement, which is the first thing checked below.

TWO CONFIGURATIONS, approached from outside the lightcone by an offset d:

  control  gamma = 0, Delta r = 1, Delta t = 1 - d.  A generic null pair, no caustic.  The
           answer is known: W ~ 1/sigma with sigma linear in d, so the power is -1.  This is
           the case that can fail, and it is here so that the caustic number means something.

  caustic  gamma = pi, Delta r = 0, Delta t = pi - d.  The image pair.  The transverse
           separation is the full pi, so the 2D part is timelike at proper time pi - d, and
           the pair goes null exactly at d = 0.

If the caustic power comes back at -1 the parametrix count survives as written.  If it comes
back steeper, the caustic strengthens the divergence rather than removing it, and Section
5.3's self-censoring reading is back with a different exponent.
"""
import numpy as np
from scipy.special import eval_legendre, kv, hankel2

A = 1.0


def mass(l):
    return np.sqrt(l * (l + 1.0)) / A


def parity_identity_holds():
    """P_l(-x) = (-1)^l P_l(x): the antipodal map IS the (-1)^l of A.8."""
    l = np.arange(0, 40)
    worst = 0.0
    for x in (0.0, 0.3, -0.77, 0.999):
        lhs = eval_legendre(l, -x)
        rhs = (-1.0) ** l * eval_legendre(l, x)
        worst = max(worst, np.max(np.abs(lhs - rhs)))
    return worst


def shell_sum(offset, caustic, eps=None, lmax=None, parity=True):
    """The multipole sum at offset d from the lightcone, in one of the two configurations.

    The caustic sum needs a regulator and there is only one honest choice.  Its terms grow
    as sqrt(l): the antipodal parity cancels the transverse phase, every multipole arrives
    in step at d = 0, and that coherence is the singularity.  A taper in l would be an
    invention and would set the answer.  The Wightman function's own i-epsilon is not an
    invention, it is the definition, and it damps the tower by exp(-m_l eps).  So eps is
    taken well inside the offset and the answer is checked for being free of it.
    """
    if eps is None:
        eps = offset / 40.0
    if lmax is None:
        lmax = int(40.0 / eps)
    l = np.arange(1, lmax + 1, dtype=float)
    m = mass(l)
    if caustic:
        tau = np.pi - offset
        w2 = -0.25j * hankel2(0, m * (tau - 1j * eps))
        angular = (-1.0) ** l if parity else np.ones_like(l)
    else:
        dr, dt = 1.0, 1.0 - offset
        rho = np.sqrt(dr * dr - dt * dt)
        w2 = kv(0, m * rho) / (2.0 * np.pi)
        angular = np.ones_like(l)
    return np.sum((2.0 * l + 1.0) / (4.0 * np.pi) * angular * w2)


def power(caustic, offsets, parity=True):
    vals = np.array([abs(shell_sum(d, caustic, parity=parity)) for d in offsets])
    slope, intercept = np.polyfit(np.log(offsets), np.log(vals), 1)
    resid = np.log(vals) - (slope * np.log(offsets) + intercept)
    return slope, float(np.max(np.abs(resid))), vals


def eps_stability(offset, caustic):
    """The i-epsilon must be small enough to have left the answer."""
    out = [abs(shell_sum(offset, caustic, eps=offset / f)) for f in (20.0, 40.0, 80.0)]
    return (max(out) - min(out)) / max(out), out


if __name__ == "__main__":
    print(__doc__)
    print("=" * 78)

    worst = parity_identity_holds()
    print(f"\n[1] antipodal map == (-1)^l insertion, worst error over 160 cases: {worst:.3e}")
    assert worst < 1e-12

    print("\n[1b] both configurations measure the same variable: sigma vanishes linearly in")
    print("     the offset in each, so the two fitted powers are comparable")
    for s in (1e-3, 1e-4, 1e-5):
        sig_ctrl = 0.5 * (1.0 ** 2 - (1.0 - s) ** 2)
        sig_caus = 0.5 * (np.pi ** 2 - (np.pi - s) ** 2)
        print(f"     offset {s:.0e}:  sigma/s = {sig_ctrl / s:.6f} control, "
              f"{sig_caus / s:.6f} image  (limits 1 and pi = {np.pi:.6f})")
        assert abs(sig_ctrl / s - 1.0) < 1e-3 and abs(sig_caus / s - np.pi) < 1e-3

    offsets = np.array([0.02, 0.014, 0.01, 0.007, 0.005, 0.0035, 0.0025])

    s_ctrl, r_ctrl, _ = power(False, offsets)
    print(f"\n[2] control, a generic null pair with no caustic")
    print(f"    fitted power {s_ctrl:+.4f}   max log-residual {r_ctrl:.4f}")
    print(f"    known answer -1, from W ~ 1/sigma with sigma linear in the offset")
    assert abs(s_ctrl + 1.0) < 0.05, "the control missed its known answer; the method is wrong"

    s_flat, r_flat, _ = power(True, offsets, parity=False)
    print(f"\n[3] the same configuration with the antipodal parity removed by hand")
    print(f"    fitted power {s_flat:+.4f}")
    print(f"    nothing should diverge: without the parity the phases never come into step")
    assert abs(s_flat) < 0.2, "the plant diverged, so the parity is not what is doing it"

    s_caus, r_caus, _ = power(True, offsets)
    print(f"\n[4] caustic, the image pair at the full pi of transverse angle")
    print(f"    fitted power {s_caus:+.4f}   max log-residual {r_caus:.4f}")

    print(f"\n[5] the i-epsilon is not setting either answer")
    for name, flag in (("control", False), ("caustic", True)):
        spread, _ = eps_stability(0.005, flag)
        print(f"    {name}: relative spread over 4x of eps = {spread:.2e}")
        assert spread < 1e-2

    print(f"\n[6] the power sharpens on -3/2 as the window closes on the caustic")
    tight, bywindow = None, []
    for lo, hi in ((0.02, 0.0025), (0.008, 0.001), (0.003, 0.0004)):
        window = np.exp(np.linspace(np.log(lo), np.log(hi), 7))
        s, r, _ = power(True, window)
        print(f"    offsets {lo:g} to {hi:g}:  power {s:+.5f}   max log-residual {r:.4f}")
        bywindow.append(abs(s))
        tight = s
    assert abs(tight + 1.5) < 5e-3, "the caustic power is not -3/2"
    print(f"    so the image term goes as s^-3/2 and not as the parametrix s^-1")

    print(f"\n[7] the coefficient of the divergence is real and positive")
    ratios = []
    for f in (20.0, 40.0, 80.0, 160.0):
        s = shell_sum(0.001, True, eps=0.001 / f)
        ph = np.angle(s) / np.pi
        ratios.append(ph * f)
        print(f"    eps = offset/{f:<5.0f}  arg/pi = {ph:+.7f}   |S| = {abs(s):9.4f}")
    drift = (max(ratios) - min(ratios)) / abs(np.mean(ratios))
    print(f"    the phase halves with eps, so it is the regulator's and not the limit's")
    print(f"    (phase x eps-factor is constant to {drift:.1e}); the limit sits on the")
    print(f"    positive real axis, so the divergence does not oscillate.  That is not")
    print(f"    the sign A.15 asks for, which needs the stress tensor's own components.")
    assert drift < 1e-2

    print(f"\n[8] a stress tensor is quadratic in derivatives, so the two-point powers")
    print(f"    {s_ctrl:+.3f} and {s_caus:+.3f} become {s_ctrl - 2:+.3f} and {s_caus - 2:+.3f}")
    print(f"\n    parametrix count for the image term, Delta taken finite:  s^-3")
    print(f"    mode sum at the caustic, no parametrix and no Delta:      s^{s_caus - 2:+.2f}")
    print(f"    the caustic steepens the divergence by {s_ctrl - s_caus:+.3f} of a power")

    print(f"\n[9] the same exponents as the p in s^-p, which is how the text reads them")
    print(f"    two-point function:  control p = {abs(s_ctrl):.4f}, caustic p = {abs(s_caus):.4f}")
    print(f"    caustic p by window: " + ", ".join(f"{v:.5f}" for v in bywindow))
    print(f"    stress tensor:       control p = {abs(s_ctrl) + 2:.3f}, "
          f"caustic p = {abs(s_caus) + 2:.3f}")
