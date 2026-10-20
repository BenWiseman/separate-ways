#!/usr/bin/env python3
"""The inflation companion of Separate Ways v6, section 3.3: where Herring and Boyanovsky's analytic result departs from ours.

Herring and Boyanovsky (PRD 101 (2020) 123522, arXiv:2005.00391) match the Bunch-Davies spinor at
the end of inflation, which for these modes is the massless one (xi, s xi)/sqrt2, to the radiation-era
"out" spinors built from h = D_alpha[sqrt2 e^{i pi/4} z], alpha = -1 - i q^2/2, at z = 0 (their
eqs. (ah), (ahpri)). That is the same calculation as our massless vacuum at the start of the
radiation era. They then write (their eq. bogofin1)
    |B|^2 = (1/2)[1 - N q sqrt(m H_R)(H* h + H h*)],
with N the out-spinor normalisation, and reduce it with Gamma-function identities to the half-angle
formula (1/2)[1 - sqrt(1 - exp(-pi q^2))].

The projection of their own spinors is B = (N s/sqrt2)(H - k h), so |B|^2 = (1/2)[1 - 2 k N^2 Re(H h*)]
(their printed eq. carries N once, a misprint: their final formula is the N^2 version). Evaluated
exactly this is the closed form
    n_vac(x) = (1/2)[1 - sqrt(1 - exp(-pi x^2)) cos Theta(x)],
    Theta(x) = pi/4 + arg Gamma(1/2 + i x^2/4) - arg Gamma(1 + i x^2/4)
(cf. Damski and Zurek 2006 for the half Landau-Zener crossing), which matches the direct mode
integration. Their Appendix C replaces the Gamma-function phases by a Stirling form, which sets
Theta -> 0 and gives the half-angle formula; but Theta(0) = pi/4 and Theta(1) = 0.47, so the
approximation fails where the abundance is made. Also: Beltran-Palau et al.'s I = 0.06305 is the
full integral 0.06419 cut off near x = 31.
Plants: the printed form goes negative; Theta -> 0 must reproduce the half-angle formula exactly.
"""
import math, sys
import mpmath as mp
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from v6_reheating_production import Profile, n_exact
fails = 0
def check(ok, what):
    global fails
    print("   %-70s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)

mp.mp.dps = 30
def hb(q):
    y = q*q/4
    h = mp.sqrt(mp.pi/2)*mp.power(2, -1j*y)/mp.gamma(1 + 1j*y)                         # their (ah), m H_R = 1
    hp = -mp.sqrt(2*mp.pi)*mp.exp(1j*mp.pi/4)*mp.power(2, -1j*y)/mp.gamma(0.5 + 1j*y)    # their (ahpri)
    H = 1j*hp                                                                            # calH = i h' + M h, M -> 0
    N2 = 1/(abs(H)**2 + q*q*abs(h)**2)                                                   # their normalisation
    printed = 0.5*(1 - mp.sqrt(N2)*q*2*(H*h.conjugate()).real)
    corrected = 0.5*(1 - 2*q*N2*(H*h.conjugate()).real)
    direct = N2*abs(H - q*h)**2/2
    return float(printed), float(corrected), float(direct), float(N2)

if __name__ == "__main__":
    prof = Profile(t_rh=1e-4)
    print("   q     as printed      corrected       our n_vac       half-angle      16 q^4 n")
    worst = 0.0; neg = False
    for q in (0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 5.0):
        pr, co, di, n2 = hb(q)
        nv = n_exact(q, prof)
        ns = (1 - math.sqrt(1 - math.exp(-math.pi*q*q)))/2
        print("  %.1f   %+.4e   %.6e   %.6e   %.4e   %.3f" % (q, pr, co, nv, ns, 16*q**4*co))
        worst = max(worst, abs(co/nv - 1)); neg |= pr < 0
        check(abs(co/di - 1) < 1e-20, "q = %.1f: the corrected form equals the direct projection" % q)
        check(abs(n2/(math.exp(-math.pi*q*q/4)/2) - 1) < 1e-12, "q = %.1f: normalisation equals their (Nwidet)" % q)
    print("   largest difference between the corrected projection and the mode integration: %.1e" % worst)
    theta = lambda x: float(mp.pi/4 + mp.arg(mp.gamma(0.5 + 1j*x*x/4)) - mp.arg(mp.gamma(1 + 1j*x*x/4)))
    closed = lambda x, th: 0.5*(1 - math.sqrt(1 - math.exp(-math.pi*x*x))*math.cos(th))
    dev = max(abs(closed(q, theta(q))/hb(q)[1] - 1) for q in (0.3, 0.6, 1.0, 1.5, 2.0, 3.0))
    check(dev < 1e-12, "closed form with Theta equals the corrected projection")
    half = lambda x: (1 - math.sqrt(1 - math.exp(-math.pi*x*x)))/2
    check(max(abs(closed(q, 0.0) - half(q)) for q in (0.3, 1.0, 2.0)) < 1e-15, "plant: Theta -> 0 gives the half-angle formula")
    print("   Theta(0.01) = %.4f (pi/4 = %.4f), Theta(0.5) = %.3f, Theta(1) = %.3f" % (theta(0.01), math.pi/4, theta(0.5), theta(1.0)))
    print("   the half-angle occupation their formula gives has C_M = I_*/pi^1.5 = %.4f" % (0.0127596673634/math.pi**1.5))
    from scipy.integrate import quad
    full = (quad(lambda x: x*x*closed(x, theta(x)), 0, 60, limit=400)[0] + 1/(16*60.0))/math.sqrt(math.pi)
    cut = quad(lambda x: x*x*closed(x, theta(x)), 0, 30.9, limit=400)[0]/math.sqrt(math.pi)
    print("   I for the massless vacuum: %.6f in full; %.6f cut off at x = 30.9 (Beltran-Palau et al. 0.06305)" % (full, cut))
    check(abs(cut - 0.06305) < 5e-5, "the cut-off integral reproduces 0.06305")
    check(worst < 2e-3, "corrected projection reproduces the massless-vacuum occupation")
    check(neg, "plant: the printed form goes negative, so the comparison can fail")
    sys.exit(1 if fails else 0)
