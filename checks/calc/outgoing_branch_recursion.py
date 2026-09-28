#!/usr/bin/env python3
"""The three-term recursion for the outgoing Frobenius branch at a Schwarzschild horizon.

Leaver's series is built on the ingoing branch.  A reflecting horizon needs the other
one, and it is not Leaver's coefficients with a sign moved: the two solutions of a
regular singular point have different recursions.  With 2M = 1, rho = -i omega and
u = (r-1)/r, the outgoing ansatz is fixed before any recursion is written, by requiring
the r-power at infinity to be that of an outgoing wave:

    psi = (r-1)^(-rho) r^0 e^(-rho (r-1)) sum_n a_n u^n,

so A = -rho and B = 0, against Leaver's A = +rho, B = -2 rho.  Under it the recursion is

    alpha_n = (n+1)(n+1-2 rho),  beta_n = -(2n^2+2n+l(l+1)-s^2+1),  gamma_n = n^2-s^2,
    alpha_n a_{n+1} + beta_n a_n + gamma_n a_{n-1} = 0.

This script takes that as given and asks the only question that matters: does the series
it generates solve the Regge-Wheeler equation?  The residual is formed against the
original equation, A psi'' + B psi' + C psi with A = r^3(r-1)^2, B = r^2(r-1),
C = w^2 r^5 - (l(l+1) r + 1 - s^2) r (r-1), and every coefficient is perturbed in turn to
confirm the test can fail.
"""
from mpmath import mp, mpc, mpf, exp, log

mp.dps = 40

LL, SPIN = mpf(6), 2
BETA = mpf(1 - SPIN * SPIN)
NTERM = 400
OMEGA = mpc("0.747343368836", "-0.177924631028")     # 2M = 1, l = s = 2, n = 0


def coefficients(rho, bump=None):
    """a_0 .. a_NTERM from the outgoing recursion; bump plants a fault in one family."""
    def alpha(n):
        v = (n + 1) * (n + 1 - 2 * rho)
        return v * mpf("1.000001") if bump == "alpha" else v

    def beta(n):
        v = -(2 * n * n + 2 * n + LL - SPIN * SPIN + 1)
        return v + mpf("1e-6") if bump == "beta" else v

    def gamma(n):
        v = mpf(n * n - SPIN * SPIN)
        return v + mpf("1e-6") if bump == "gamma" else v

    a = [mpc(1), -beta(0) / alpha(0)]
    for n in range(1, NTERM):
        a.append(-(beta(n) * a[n] + gamma(n) * a[n - 1]) / alpha(n))
    return a


def residual(a, rho, r):
    """Relative size of A psi'' + B psi' + C psi for the series solution at r."""
    omega = 1j * rho
    u = (r - 1) / r
    s0 = s1 = s2 = mpc(0)
    for n in range(len(a) - 1, -1, -1):
        s0 = s0 * u + a[n]
        if n >= 1:
            s1 = s1 * u + n * a[n]
        if n >= 2:
            s2 = s2 * u + n * (n - 1) * a[n]
    tail = abs(a[-1]) * abs(u) ** (len(a) - 1) / max(abs(s0), mpf("1e-300"))
    pref = exp(-rho * log(r - 1) - rho * (r - 1))
    g = -rho / (r - 1) - rho
    gp = rho / (r - 1) ** 2
    du, ddu = 1 / r ** 2, -2 / r ** 3
    psi = pref * s0
    dpsi = pref * (g * s0 + s1 * du)
    ddpsi = pref * ((gp + g * g) * s0 + 2 * g * s1 * du + s2 * du * du + s1 * ddu)
    A = r ** 3 * (r - 1) ** 2
    B = r ** 2 * (r - 1)
    C = omega ** 2 * r ** 5 - (LL * r + BETA) * r * (r - 1)
    scale = max(abs(A * ddpsi), abs(B * dpsi), abs(C * psi))
    return abs(A * ddpsi + B * dpsi + C * psi) / scale, tail


if __name__ == "__main__":
    rho = -1j * OMEGA
    print("Outgoing branch, l = 2, s = 2, 2M = 1, rho = -i omega at the fundamental mode")
    print("  %d terms; residual is |A psi'' + B psi' + C psi| over the largest of the three\n"
          % NTERM)
    good = coefficients(rho)
    for r in (mpf("1.5"), mpf("2.5"), mpf("6"), mpc("3", "2")):
        res, tail = residual(good, rho, r)
        print("  r = %-10s residual %-12s series tail %s"
              % (mp.nstr(r, 6), mp.nstr(res, 4), mp.nstr(tail, 4)))

    print("\n  the same test at a frequency the mode condition does not pick out,")
    print("  because a recursion is a property of the equation and not of its spectrum:")
    for omega in (mpc("0.4", "0"), mpc("0.7", "-0.1")):
        res, _ = residual(coefficients(-1j * omega), -1j * omega, mpf("2.5"))
        print("    omega = %-18s residual %s"
              % (mp.nstr(mpc(omega), 6), mp.nstr(res, 4)))

    print("\n  planted faults, each of which must break it:")
    for bump in ("alpha", "beta", "gamma"):
        res, _ = residual(coefficients(rho, bump), rho, mpf("2.5"))
        print("    %-8s perturbed   residual %s" % (bump, mp.nstr(res, 4)))
    res, _ = residual(good, rho, mpf("2.5"))
    print("    %-8s             residual %s" % ("unperturbed", mp.nstr(res, 4)))

    print("\n  and Leaver's own coefficients under this prefactor, which must also fail:")
    def leaver_a(rho):
        a = [mpc(1)]
        def al(n): return (n + 1) * (n + 1 - 2 * rho)
        def be(n): return -(2 * n * n + 2 * n + LL - SPIN * SPIN + 1) - 4 * rho * (n + 1) + 4 * rho * rho
        def ga(n): return n * n - 4 * rho * n + 4 * rho * rho - SPIN * SPIN
        a.append(-be(0) / al(0))
        for n in range(1, NTERM):
            a.append(-(be(n) * a[n] + ga(n) * a[n - 1]) / al(n))
        return a
    res, _ = residual(leaver_a(rho), rho, mpf("2.5"))
    print("    ingoing coefficients  residual %s" % mp.nstr(res, 4))
