#!/usr/bin/env python3
"""Separate Ways v6, Section 5 and Appendix B: does the plasma make N1 if gravity is classical?

The rate alpha T^8/M_P^4 of v6_grav_rate.py is the standard low-energy treatment of gravity: the
metric's small fluctuations are quantised, and the plasma makes an N1 pair when it hands a graviton
energy k0 > 0. In correlator language the rate weighs the pair-production spectrum of N1 with the
plasma's stress-tensor Wightman function A(k), the spectral weight for the plasma to GIVE energy k0:
for two-body states A carries f1 f2. The reverse weight E(k), for the plasma to TAKE energy k0,
carries (1 +- f1)(1 +- f2) and obeys E = exp(k0/T) A (KMS).

If the metric is classical, three prescriptions are on the market:
  1. mean field (semiclassical gravity): the metric follows <T> only. A homogeneous plasma then
     gives the expanding background and nothing else, so this channel is absent.
  2. stochastic gravity (Hu and Verdaguer): the metric also follows a classical noise whose spectrum
     is the symmetrised correlator, (A + E)/2. A classical field cannot tell giving from taking, so
     N1 pairs are made from E as well. E contains the vacuum term (no f): the rate then grows with
     the momentum cutoff Lambda as Lambda^8/M_P^4, and N1 is made from empty space.
  3. a classical metric that follows only the plasma's thermal, two-particle fluctuations, the
     f1 f2 parts of A and E: (A + E_ff)/2 = A, so the rate equals alpha.
With x = E/T and the s^2 growth of every channel, each two-body weight factorises into moments
  J = Int x^3 f dx (pi^4/15 for bosons, 7 pi^4/120 for fermions) and L = Int_0^(Lambda/T) x^3 dx,
so the symmetrised rate is alpha + (1/2) sum_i N_i c_i (L^2 +- 2 L J_i), with + for bosons and - for
fermions ((1 - f1)(1 - f2) for fermions); each plasma particle on a vacuum leg is cut at energy Lambda in the
plasma's frame, a choice the numbers depend on. N1 is taken massless at these temperatures.
Plants: the KMS identity must hold point by point; with no cutoff (L = 0) prescription 2 must return
alpha; at large Lambda its rate must grow as Lambda^8.
"""
import math, sys
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import numpy as np
import v6_grav_rate as G
fails = 0
def check(ok, what):
    global fails
    print("   %-74s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)

xs, ws = np.polynomial.legendre.leggauss(24)
ang = lambda fn: 2*math.pi*sum(w*fn(x) for x, w in zip(xs, ws))
pref = (32/3)/(1024*math.pi**6)                 # Int dcos12 (1 - cos12)^2 = 8/3, times 4 from s^2 = 4 E1^2 E2^2 (...)^2
JB, JF = math.pi**4/15, 7*math.pi**4/120
# per unit of each plasma species: angular integral of the spin sum into one Dirac N pair (one Majorana N1
# per Dirac-pair event), 1/2 for identical initial bosons, a Weyl fermion half a Dirac spin sum
c = {"scalar": 0.5*ang(lambda x: G.spin_sum("scalar", 1.0, x)),
     "fermion": 0.5*ang(lambda x: G.spin_sum("fermion", 1.0, x)),
     "vector": 0.5*ang(lambda x: G.spin_sum("vector", 1.0, x))}
N = {"scalar": 4, "fermion": 45, "vector": 12}
J = {"scalar": JB, "fermion": JF, "vector": JB}
sign = {"scalar": +1, "fermion": -1, "vector": +1}

def rate(L, J_scale=1.0):
    """symmetrised-noise rate in units T^8/M_P^4; L = Int_0^(Lambda/T) x^3 dx"""
    A = sum(N[k]*c[k]*pref*J[k]**2 for k in N)
    E = sum(N[k]*c[k]*pref*(L*L + 2*sign[k]*L*J[k] + J[k]**2) for k in N)
    return A, 0.5*(A + E)

if __name__ == "__main__":
    alpha, _ = rate(0.0)
    print("=== N1 from the plasma under three treatments of the metric ===")
    print("   quantised fluctuations (standard treatment): alpha = %.4e T^8/M_P^4" % alpha)
    check(abs(alpha/1.2027e-2 - 1) < 1e-3, "reproduces alpha of v6_grav_rate.py")
    # KMS, point by point: (1 +- f1)(1 +- f2) exp(-(x1 + x2)) = f1 f2
    fb = lambda x: 1/math.expm1(x); ff = lambda x: 1/(math.exp(x) + 1)
    kms = max(max(abs((1 + fb(a))*(1 + fb(b))*math.exp(-a - b)/(fb(a)*fb(b)) - 1),
                  abs((1 - ff(a))*(1 - ff(b))*math.exp(-a - b)/(ff(a)*ff(b)) - 1)) for a in (0.3, 1.0, 4.0) for b in (0.5, 2.0, 7.0))
    check(kms < 1e-12, "plant: the weights obey KMS, E = exp(k0/T) A, point by point")
    print("   1. mean field: the metric follows <T> only; a homogeneous plasma makes no N1 this way: rate 0")
    _, th = rate(0.0)
    print("   3. thermal two-particle noise, (A + E_ff)/2: %.4e T^8/M_P^4, %.3f alpha" % (th, th/alpha))
    check(abs(th/alpha - 1) < 1e-12, "thermal two-particle noise gives alpha exactly")
    mb = sum(N[k]*c[k]*pref*36.0 for k in N)        # Maxwell-Boltzmann moments, Int x^3 e^-x dx = 6
    print("      with Maxwell-Boltzmann statistics, a classical gas: %.4e, %.3f alpha" % (mb, mb/alpha))
    print("   2. stochastic gravity, symmetrised noise (A + E)/2, momentum cutoff Lambda:")
    out = []
    for lam in (1.0, 2.0, 3.0, 5.0, 10.0, 30.0):
        L = lam**4/4
        _, s = rate(L)
        out.append((lam, s))
        print("      Lambda = %5.1f T: rate %.3e T^8/M_P^4 = %.3g alpha" % (lam, s, s/alpha))
    check(out[-1][1]/rate((out[-1][0]/2)**4/4)[1] > 200, "at large Lambda the symmetrised rate grows as Lambda^8")
    vac = 0.5*sum(N[k]*c[k]*pref for k in N)*(1/4.0)**2
    print("      the vacuum part alone, at zero temperature: %.3e Lambda^8/M_P^4 with this plasma-frame cutoff; the vacuum"
          % vac)
    print("      noise is Lorentz invariant, so no Lorentz-invariant cutoff keeps that rate finite")
    print("   reading: whether the channel exists is a question about the quantum nature of gravity; a classical")
    print("   metric must say which fluctuations it follows")
    sys.exit(1 if fails else 0)
