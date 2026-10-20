#!/usr/bin/env python3
"""Separate Ways v6, section 4: the abundance-to-mass relation, derived and not cited.

A Majorana fermion of mass M made gravitationally in the radiation era, with occupation n(x) in
x = sqrt(pi) k / sqrt(M a1), has comoving number
    n a^3 = 2 Int d^3k/(2 pi)^3 n_k = (M a1/pi)^(3/2) I,   I = (1/pi^2) Int x^2 n(x) dx,
where a1 = a'(eta) in the radiation era. Friedmann in conformal time gives
    a1 = (aT)^2 pi sqrt(g*/15) / mu,   mu = (4 pi G/3)^(-1/2) = sqrt(3/(4 pi)) M_Pl,
and s a^3 = (2 pi^2/45) g* (aT)^3. The (aT)^3 cancels, and the measured rho_DM/s fixes
    M^(5/2) = (rho_DM/s0) (2 pi^2/45) g*^(1/4) 15^(3/4) mu^(3/2) / I.
"""
import math, sys
fails = 0
def check(ok, what):
    global fails
    print("   %-66s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)

M_PL = 1.22089e19                      # GeV, Planck mass (CODATA G); reduced = M_PL/sqrt(8 pi)
MU = math.sqrt(3/(4*math.pi))*M_PL     # (4 pi G/3)^(-1/2)
GEV_PER_INVCM3 = (1.97327e-14)**3      # 1 cm^-3 in GeV^3 (hbar c = 1.97327e-14 GeV cm)
RHO_CRIT_H2 = 1.05368e-5               # GeV cm^-3, critical density / h^2 (PDG)
OMEGA_C_H2, OMEGA_ERR = 0.1200, 0.0012 # Planck 2018 cold dark matter
S0_CM3 = 2891.2                        # cm^-3, today's entropy density (T_CMB = 2.7255 K)
GSTAR = 106.75
I_MIN = 0.0127596673634                # Boyle-Finn-Turok half-angle state

def mass(I, g=GSTAR, omega=OMEGA_C_H2, s0_cm3=S0_CM3):
    rho = omega*RHO_CRIT_H2*GEV_PER_INVCM3
    s0 = s0_cm3*GEV_PER_INVCM3
    M52 = (rho/s0)*(2*math.pi**2/45)*g**0.25*15**0.75*MU**1.5/I
    return M52**0.4

if __name__ == "__main__":
    rho = OMEGA_C_H2*RHO_CRIT_H2*GEV_PER_INVCM3; s0 = S0_CM3*GEV_PER_INVCM3
    print("rho_DM = %.4e GeV^4,  s0 = %.4e GeV^3,  mu = %.4e GeV" % (rho, s0, MU))
    Mc = mass(I_MIN)
    print("\n=== the ceiling ===")
    print("   M_c = %.4e GeV = %.1f PeV,  half-mass line %.1f PeV" % (Mc, Mc/1e6, Mc/2e6))
    dM = 0.4*(OMEGA_ERR/OMEGA_C_H2)*Mc
    print("   width from Omega_c h^2 alone: +- %.1f PeV (mass), +- %.1f PeV (line)" % (dM/1e6, dM/2e6))
    check(abs(Mc/4.916e8 - 1) < 2e-3, "reproduces the 4.916e8 GeV of the release")
    print("\n=== inputs that move the ceiling (all upward or sideways, none widen it) ===")
    s0_alt = 2.3e-38/GEV_PER_INVCM3
    print("   s0 = 2.3e-38 GeV^3 (Boyle et al.): M = %.1f PeV" % (mass(I_MIN, s0_cm3=s0_alt)/1e6))
    for g, lab in [(110.25, "N2, N3 relativistic at production (+2 x 2 x 7/8)")]:
        print("   g* = %.2f, %s: M = %.1f PeV (%+.1f)" % (g, lab, mass(I_MIN, g=g)/1e6, (mass(I_MIN, g=g)-Mc)/1e6))
    print("   an entropy release by a factor F after production: M -> M F^(2/5)... e.g. F = 2: %.1f PeV" % (Mc*2**0.4/1e6))
    # the paper quotes the ceiling for Boyle, Finn and Turok's rounded rho_DM = 9.7e-48 GeV^4, so the
    # numbers that move it are given from that base too
    om_r = 9.7e-48/(RHO_CRIT_H2*GEV_PER_INVCM3)
    Mr = mass(I_MIN, omega=om_r)
    print("\n=== from the rounded rho_DM = 9.7e-48 GeV^4 (Omega_c h^2 = %.5f) ===" % om_r)
    print("   M_* = %.1f PeV; with s0 = 2.3e-38: %.1f PeV; g* = 110.25: %+.1f PeV; F = 2: %.1f PeV"
          % (Mr/1e6, mass(I_MIN, omega=om_r, s0_cm3=s0_alt)/1e6, (mass(I_MIN, g=110.25, omega=om_r) - Mr)/1e6, Mr*2**0.4/1e6))
    check(abs(Mr/4.916e8 - 1) < 1e-3, "rounded input gives the 491.6 PeV the paper quotes")
    print("   a two-body decay to a neutrino puts the line at M_*/2 = %.1f PeV" % (Mr/2e6))
    # plant: a planted change in the abundance must move the mass at the 2/5 power
    check(abs(mass(I_MIN, omega=2*OMEGA_C_H2)/Mc - 2**0.4) < 1e-12, "plant: doubling the abundance raises M by 2^(2/5)")
    check(abs(mass(2*I_MIN)/Mc - 2**-0.4) < 1e-12, "plant: doubling production lowers M by 2^(-2/5)")
    sys.exit(1 if fails else 0)
