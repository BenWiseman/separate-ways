#!/usr/bin/env python3
"""Response of the fundamental Schwarzschild ringdown frequency to a reflecting horizon.

The boundary convention, stated because the answer depends on it: at the horizon the
Regge-Wheeler function is written in powers of x = r - 1 (units 2M = 1) as

    psi ~ x^rho + R_h x^(-rho),        rho = -i omega,

with unit leading coefficients on both Frobenius solutions and the principal branch of
log x.  x^rho is the wave transmitted into the horizon and x^(-rho) the reflected wave,
so R_h is the horizon reflectivity in the power-coefficient normalisation.  Measured
instead against unit-amplitude waves in r_* = r + log(r-1), the same reflectivity is
R_wave = e^(2 rho) R_h, because (r-1)^rho = e^(-rho) e^(rho r_*); the two derivatives
therefore differ by e^(-2 rho) and both are reported.  Shifting the additive origin of
r_* moves the second one again, which is why the first is the one to quote.

Method, which shares nothing with a continued fraction.  Let psi be the solution
outgoing at infinity, continued to the horizon and decomposed there as
psi = A_in x^rho(1 + ...) + A_out x^(-rho)(1 + ...).  A horizon of reflectivity R_h
gives the eigenvalue condition S(omega) := A_out/A_in = R_h, so the quasinormal mode is
the root of S and

    d omega / d R_h = 1 / S'(omega_0).

S is built by marching the exact Taylor series of the Regge-Wheeler equation inward from
r = 1 + T e^(i theta) to the matching point and projecting onto the two Frobenius
solutions with a Wronskian.  The ray is rotated off the real axis so that the outgoing
solution is the dominant one in the direction of integration: whatever the starting
asymptotic series gets wrong is then suppressed by e^(-120) on the way in, and the part
of it that is not suppressed is an overall scale, which cancels in A_out/A_in.
"""
from mpmath import mp, mpc, mpf, exp, log, pi

mp.dps = 50

LL = mpf(6)                  # l(l+1), l = 2
BETA = mpf(-3)               # (1 - s^2) 2M with s = 2: coefficient of 1/r^3 in V/f
THETA = mp.pi / 4
T_START = mpf(300)
T_MATCH = mpf("0.5")
NSER = 170
LEAVER = mpc("0.747343368836", "-0.177924631028")     # 2M = 1, l = s = 2, n = 0


def ray(t):
    return 1 + t * exp(1j * THETA)


def operator(omega):
    """A psi'' + B psi' + C psi = 0 with A = r^3(r-1)^2, B = r^2(r-1),
    C = w^2 r^5 - (L r + beta) r (r-1); coefficients ascending in r."""
    A = [mpf(0)] * 3 + [mpf(1), mpf(-2), mpf(1)]
    B = [mpf(0), mpf(0), mpf(-1), mpf(1)]
    C = [mpf(0), BETA, LL - BETA, -LL, mpf(0), omega ** 2]
    return A, B, C


def shift(poly, r0):
    """Re-expand a polynomial about r0 by repeated synthetic division."""
    work, out = list(poly), []
    while work:
        acc = mpc(0)
        for c in reversed(work):
            acc = acc * r0 + c
        out.append(acc)
        rest, carry = [], mpc(0)
        for c in reversed(work[1:]):
            carry = carry * r0 + c
            rest.append(carry)
        work = list(reversed(rest))
    return out


def march(omega, psi, dpsi, r_from, r_to):
    """Analytic continuation along the straight segment by Taylor steps."""
    A, B, C = operator(omega)
    direction = (r_to - r_from) / abs(r_to - r_from)
    r, remaining, worst = r_from, abs(r_to - r_from), mpf(0)
    while remaining > mpf("1e-30"):
        h = min(mpf("0.45") * min(abs(r), abs(r - 1)), 6 / abs(omega), remaining)
        As, Bs, Cs = shift(A, r), shift(B, r), shift(C, r)
        a = [psi, dpsi]
        for m in range(NSER):
            acc = mpc(0)
            for k in range(1, min(len(As), m + 3)):
                acc += As[k] * (m - k + 2) * (m - k + 1) * a[m - k + 2]
            for k in range(0, min(len(Bs), m + 2)):
                acc += Bs[k] * (m - k + 1) * a[m - k + 1]
            for k in range(0, min(len(Cs), m + 1)):
                acc += Cs[k] * a[m - k]
            a.append(-acc / (As[0] * (m + 2) * (m + 1)))
        z = h * direction
        val, der = mpc(0), mpc(0)
        for n in range(len(a) - 1, -1, -1):
            val = val * z + a[n]
            if n:
                der = der * z + n * a[n]
        worst = max(worst, abs(a[-1] * z ** (len(a) - 1)) / max(abs(val), mpf("1e-300")))
        psi, dpsi, r, remaining = val, der, r + z, remaining - h
    return psi, dpsi, worst


def outer_start(omega):
    """The solution outgoing at infinity, e^(i w r)(r-1)^(-rho) sum b_k r^-k, optimally
    truncated.  The constant e^(i w r_start) is dropped; it cancels in A_out/A_in."""
    r, rho = ray(T_START), -1j * omega
    b, size = [mpc(1)], [mpf(1)]
    for m in range(1, 240):
        prev2 = b[m - 2] if m >= 2 else mpc(0)
        b.append(((m * (m - 1) - LL) * b[m - 1] + (-m * m + 2 * m - BETA) * prev2)
                 / (2j * omega * m))
        size.append(abs(b[m]) / abs(r) ** m)
    # b_3 vanishes identically at l = s = 2, so take the minimum of the pairwise
    # maximum: an isolated zero must not be mistaken for the optimal truncation.
    cut = min(range(1, len(size) - 1), key=lambda k: max(size[k], size[k + 1]))
    u, du = mpc(0), mpc(0)
    for k in range(cut, -1, -1):
        u = u / r + b[k]
        du = du / r + (-k * b[k] / r)
    du = du / r
    pref = exp(-rho * log(r - 1))
    return pref * u, pref * ((1j * omega - rho / (r - 1)) * u + du), cut, size[cut]


def frobenius(omega, sign):
    """x^(sign rho)(1 + ...) and its r-derivative at the matching point."""
    sigma = sign * (-1j * omega)
    p2 = {3: mpf(3), 4: mpf(3), 5: mpf(1)}
    p1 = {2: mpf(2), 3: mpf(1)}
    binom5 = [mpf(1), mpf(5), mpf(10), mpf(10), mpf(5), mpf(1)]
    p0 = {k: omega ** 2 * binom5[k] for k in range(6)}
    p0[1] -= LL + BETA
    p0[2] -= 2 * LL + BETA
    p0[3] -= LL
    c = [mpc(1)]
    for m in range(1, NSER):
        acc = mpc(0)
        for k, w in p2.items():
            j = m - k + 2
            if j >= 0:
                acc += w * (j + sigma) * (j + sigma - 1) * c[j]
        for k, w in p1.items():
            j = m - k + 1
            if j >= 0:
                acc += w * (j + sigma) * c[j]
        for k in range(1, 6):
            j = m - k
            if j >= 0:
                acc += p0[k] * c[j]
        c.append(-acc / ((m + sigma) ** 2 + omega ** 2))
    x = T_MATCH * exp(1j * THETA)
    series, dseries = mpc(0), mpc(0)
    for n in range(NSER - 1, -1, -1):
        series = series * x + c[n]
        if n:
            dseries = dseries * x + n * c[n]
    pref = exp(sigma * log(x))
    tail = abs(c[NSER - 1]) * abs(x) ** (NSER - 1) / max(abs(series), mpf("1e-300"))
    return pref * series, pref * (sigma / x * series + dseries), tail


def reflection_ratio(omega):
    psi, dpsi, _, _ = outer_start(omega)
    psi, dpsi, _ = march(omega, psi, dpsi, ray(T_START), ray(T_MATCH))
    pin, dpin, _ = frobenius(omega, +1)
    pout, dpout, _ = frobenius(omega, -1)
    wronskian = pin * dpout - dpin * pout
    return (dpsi * pin - psi * dpin) / (psi * dpout - dpsi * pout)


def secant(f, z0, z1, tol=mpf("1e-32")):
    f0, f1 = f(z0), f(z1)
    for _ in range(40):
        if f1 == f0:
            break
        z0, f0, z1, f1 = z1, f1, z1 - f1 * (z1 - z0) / (f1 - f0), None
        f1 = f(z1)
        if abs(z1 - z0) < tol:
            break
    return z1, f1


def fmt(z, digits=10):
    z = mpc(z)
    return "%s %s %si" % (mp.nstr(z.real, digits), "+" if z.imag >= 0 else "-",
                          mp.nstr(abs(z.imag), digits))


if __name__ == "__main__":
    print("Regge-Wheeler, l = 2, s = 2, units 2M = 1; ray at %s degrees out to T = %s"
          % (mp.nstr(THETA * 180 / pi, 4), mp.nstr(T_START, 4)))
    psi0, dpsi0, cut, cut_size = outer_start(LEAVER)
    _, _, tail = march(LEAVER, psi0, dpsi0, ray(T_START), ray(T_MATCH))
    _, _, ftail = frobenius(LEAVER, +1)
    print("  outer asymptotic series truncated at k = %d, smallest term %s"
          % (cut, mp.nstr(cut_size, 4)))
    print("  worst relative Taylor tail over the march %s; Frobenius relative tail %s"
          % (mp.nstr(tail, 4), mp.nstr(ftail, 4)))

    print("\nQuasinormal frequency as the root of S, with no continued fraction anywhere:")
    root, resid = secant(reflection_ratio, LEAVER, LEAVER * mpf("1.000001"))
    print("  omega_0   = %s   (2M = 1)" % fmt(root, 16))
    print("  M omega_0 = %s" % fmt(root / 2, 14))
    print("  |S(omega_0)| = %s" % mp.nstr(abs(resid), 4))
    print("  distance from Leaver's tabulated 0.747343 - 0.177925i: %s"
          % mp.nstr(abs(root - LEAVER), 4))

    delta = mpf("1e-10")
    d1 = (reflection_ratio(root + delta) - reflection_ratio(root - delta)) / (2 * delta)
    d2 = (reflection_ratio(root + 2 * delta) - reflection_ratio(root - 2 * delta)) / (4 * delta)
    sprime = (4 * d1 - d2) / 3
    dw = 1 / sprime
    rho0 = -1j * root
    print("\nResponse to a reflecting horizon:")
    print("  S'(omega_0)      = %s   (two step sizes differ by %s)"
          % (fmt(sprime, 12), mp.nstr(abs(d1 - d2), 3)))
    print("  M domega/dR_h    = %s" % fmt(dw / 2, 12))
    print("  |M domega/dR_h|  = %s" % mp.nstr(abs(dw / 2), 12))
    print("  e^(-2 rho_0)     = %s" % fmt(exp(-2 * rho0), 12))
    print("  M domega/dR_wave = %s" % fmt(dw * exp(-2 * rho0) / 2, 12))

    pct_re = abs((dw / 2).real) / abs((root / 2).real)
    pct_im = abs((dw / 2).imag) / abs((root / 2).imag)
    print("\n  per one per cent of reflectivity the real frequency moves by %s per cent"
          % mp.nstr(100 * pct_re / 100, 3))
    print("  and the damping rate by %s per cent, so the damping answers %s times harder"
          % (mp.nstr(100 * pct_im / 100, 3), mp.nstr(pct_im / pct_re, 4)))

    print("\nWhat a ringdown measurement would see, which is where Section 4 needs this.")
    print("  tau = -1/Im(omega), so d(tau)/tau = -Im(d omega/dR) R / Im(omega_0):")
    dre, dim = (dw / 2).real, (dw / 2).imag
    w0re, w0im = (root / 2).real, (root / 2).imag
    dtau = -dim / w0im
    dfreq = dre / w0re
    print("    d(Re omega)/Re(omega_0) = %s per unit R" % mp.nstr(dfreq, 8))
    print("    d(tau)/tau              = %s per unit R" % mp.nstr(dtau, 8))
    print("    ratio                   = %s" % mp.nstr(abs(dtau / dfreq), 6))
    print("  The damping response is LINEAR in R and carries a sign, so a reflecting horizon")
    print("  can shorten the ringdown as well as lengthen it. The absorption model in common")
    print("  use, tau -> tau/(1-|R|^2), is quadratic and signless, and it is wrong at first")
    print("  order. Comparing the two at the magnitudes that get quoted:")
    print("     |R|      linear term      quadratic model     ratio")
    for R in (mpf("0.01"), mpf("0.05"), mpf("0.1"), mpf("0.35")):
        lin, quad = abs(dtau) * R, R ** 2 / (1 - R ** 2)
        print("     %-8s %-16s %-19s %s"
              % (mp.nstr(R, 3), mp.nstr(lin, 6), mp.nstr(quad, 6), mp.nstr(lin / quad, 4)))
    need = mpf("0.14") / abs(dtau)
    print("  A damping-time excess of 0.14 needs |R| = %s at linear order, against %s"
          % (mp.nstr(need, 4), mp.nstr(mp.sqrt(mpf("0.14") / mpf("1.14")), 4)))
    print("  from the quadratic model: a factor %s apart, and the linear one fixes the sign."
          % mp.nstr(mp.sqrt(mpf("0.14") / mpf("1.14")) / need, 4))

    print("\nMirror symmetry omega -> -conj(omega) at finite reflectivity:")
    print("  a real reflectivity must keep the pair, an imaginary one must split it")
    for label, refl in (("R = 0.001 ", mpc("0.001", "0")), ("R = 0.01  ", mpc("0.01", "0")),
                        ("R = 0.001i", mpc("0", "0.001")), ("R = 0.01i ", mpc("0", "0.01"))):
        keep = THETA
        THETA = mp.pi / 4
        up, _ = secant(lambda w: reflection_ratio(w) - refl, root, root * mpf("1.000001"))
        THETA = -mp.pi / 4
        mirror_seed = -mp.conj(root)
        dn, _ = secant(lambda w: reflection_ratio(w) - refl, mirror_seed,
                       mirror_seed * mpf("1.000001"))
        THETA = keep
        print("  %s  defect |omega_+ + conj(omega_-)| = %s"
              % (label, mp.nstr(abs(up + mp.conj(dn)), 4)))

    print("\nThe same number off a different contour, since the contour is not physics:")
    keep_t, keep_theta = T_START, THETA
    T_START, THETA = mpf(120), mp.pi / 3
    alt_root, _ = secant(reflection_ratio, LEAVER, LEAVER * mpf("1.000001"))
    a1 = (reflection_ratio(alt_root + delta) - reflection_ratio(alt_root - delta)) / (2 * delta)
    a2 = (reflection_ratio(alt_root + 2 * delta) - reflection_ratio(alt_root - 2 * delta)) / (4 * delta)
    alt = 1 / ((4 * a1 - a2) / 3)
    T_START, THETA = keep_t, keep_theta
    print("  T = 120, ray at 60 degrees: M domega/dR_h = %s" % fmt(alt / 2, 12))
    print("  moves by %s in the frequency, %s in the root itself"
          % (mp.nstr(abs(alt - dw) / 2, 3), mp.nstr(abs(alt_root - root), 3)))

    print("\nPlanted failures, each of which must destroy the root:")
    baseline = abs(resid)
    for label, name, value in (("scalar potential, s = 0 not s = 2", "BETA", mpf(1)),
                               ("wrong multipole, l = 3 not l = 2", "LL", mpf(12)),
                               ("sign of the 1/r^3 term flipped", "BETA", mpf(3))):
        keep = globals()[name]
        globals()[name] = value
        print("  %-36s |S(omega_0)| = %s" % (label, mp.nstr(abs(reflection_ratio(root)), 4)))
        globals()[name] = keep
    print("  %-36s |S(omega_0)| = %s" % ("unperturbed", mp.nstr(baseline, 4)))
