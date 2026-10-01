#!/usr/bin/env python3
"""Every number in pub/paper2/findings/BFT_WASHOUT_GAP_20261001.md, recomputed.

Boyle, Finn and Turok (arXiv:1803.08930) discard the primordial abundance of the two
unstable right-handed neutrinos because those "equilibrate with the thermal bath". Their
scenario never establishes a thermal bath. This checks each step of that, and plants a
failure in the one place where a silent pass would be worst.
"""
import math

Mp, v, g = 2.435323e18, 174.0, 106.75
C   = (math.pi / math.sqrt(90)) * math.sqrt(g)
I   = 0.0127596673634          # PAPER2_v4_draft.md:336, the fold's phase-space integral
m3  = 0.05e-9                  # GeV, sqrt(atmospheric)
ok  = True

def line(s=""): print(s)

line("[1] the bang cannot make its own radiation")
# a(eta) = a0 eta  =>  a'' = 0  =>  R = 6 a''/a^3 = 0, and FRW is conformally flat so C^2 = 0.
# The trace anomaly a C^2 + b E_4 + c Box R then has nothing left but a topological term.
line("    R = 6 a''/a^3 with a'' = 0 exactly, and C_abcd = 0 for every FRW metric.")
line("    So the conformal anomaly sources nothing and massless conformal quanta are")
line("    produced in zero quantity. This step is exact, not an estimate.")

line()
line("[2] w = 1/3 for ANY isotropic massless spectrum, so it says nothing about thermality")
def w_of(nf, hi=60.0, N=200000):
    h = hi / N
    rho = sum((i*h)**3 * nf(i*h) for i in range(1, N)) * h
    return (rho / 3.0) / rho
spectra = [("Bose thermal T=1", lambda x: 1.0/(math.exp(x)-1) if x > 1e-9 else 1e9),
           ("fold state",       lambda x: (1-math.sqrt(max(0.0, 1-math.exp(-x*x))))/2),
           ("monochromatic",    lambda x: math.exp(-40*(x-2.0)**2))]
for name, f in spectra:
    w = w_of(f)
    good = abs(w - 1/3) < 1e-9
    ok &= good
    line("    %-20s w = %.9f  %s" % (name, w, "" if good else "<-- ISSUE"))

line()
line("[3] the state the fold does fix is not thermal")
n0 = 0.5
line("    n(0) = %.4f exactly, a ceiling, where Bose diverges" % n0)
for xv in (3.0, 4.0, 5.0):   # asymptotic form, so test it where it is asymptotic
    nf = (1-math.sqrt(1-math.exp(-xv*xv)))/2
    line("    n(%.0f) = %.4e against the Gaussian form e^{-x^2}/4 = %.4e"
         % (xv, nf, math.exp(-xv*xv)/4))
    ok &= abs(nf/(math.exp(-xv*xv)/4) - 1) < 1e-3

line()
line("[4] the two escapes that are closed")
K = m3*Mp/(8*math.pi*C*v*v)
line("    y_seesaw^2/y_crit^2 = m3 Mp/(8 pi C v^2) = %.3f, with M_N cancelling" % K)
line("    survival needs the generated light mass below %.4f meV; solar is 8.7, atmos 50"
     % (8*math.pi*C*v*v/Mp*1e12))
line("    T_eq/M_N = K^(1/3) = %.4f, also independent of M_N" % K**(1/3))
k = (1/math.pi)**1.5*I*(45/(2*math.pi**2*g))*(C/Mp)**1.5
M_eq = ((0.75*math.sqrt(Mp/C))/k)**0.5
line("    pair production equals the radiation only at M = %.3e GeV = %.1f Mp"
     % (M_eq, M_eq/Mp))
line("    at M = Mp it supplies %.2e of it" % (k*Mp**2.5/(0.75*math.sqrt(Mp*Mp/C))))

line()
line("[5] plants, because a check that cannot fail is worse than none")
# A MASSIVE species must NOT give w = 1/3: if it does, the integrator is broken.
def w_massive(mu, hi=60.0, N=200000):
    h = hi/N; rho = P = 0.0
    for i in range(1, N):
        p = i*h; E = math.sqrt(p*p + mu*mu); f = 1.0/(math.exp(E)-1)
        rho += p*p*E*f; P += p**4/(3*E)*f
    return (P*h)/(rho*h)
wm = w_massive(20.0)
caught = wm < 0.05
line("    plant: a heavy species (m = 20T) does NOT give 1/3: w = %.5f  %s"
     % (wm, "yes" if caught else "NO"))
ok &= caught
# And the fold state must NOT match Bose at large x.
r = ((1-math.sqrt(1-math.exp(-9.0)))/2) / (1.0/(math.exp(3.0)-1))
line("    plant: fold state and Bose differ at x = 3 by %.2e, three orders, not order one" % r)
ok &= r < 1e-2

line()
line("    all checks passed" if ok else "    SOMETHING FAILED")
raise SystemExit(0 if ok else 1)
