#!/usr/bin/env python3
"""Separate Ways v6, section 6: what a neutrino telescope would see of the decay line.

Inputs, each from a cited source (v6 notes, lit_A3):
  NFW halo r_s = 24 kpc, rho_sun = 0.4 GeV/cm^3, R_sun = 8.5 kpc, truncated at 260 kpc.
  Line fraction 0.231 per primary neutrino at M = 5e8 GeV (Bauer, Rodd, Webber 2021, HDMSpectra data
  file), times 1/2 for the branching to nu Z and nu h (the other half, l W, has no line).
  IceCube EHE 12.6 yr sky-averaged all-flavour area: 1.01e3, 1.85e3, 5.13e3 m^2 at 100, 220, 1000 PeV.
  IceCube-Gen2 optical: 5.5e3, 1.07e4, 2.48e4 m^2 at the same energies.
  KM3-230213A at (l, b) = (216.06, -11.13) deg.
  Lifetime scale 1e29 s: air-shower photon limits on hadronic decays give tau >~ 1e29 s for
  1e8 - 1e10 GeV (Ishiwata, Macias, Ando, Arimoto 2020, sec. 3.2, bb channel), and for N1's own
  branching at 440 PeV, 1e29 s is in tension with the EAS-MSU and Auger photon limits (Borah, Das,
  Okada, Sarmah 2025, sec. IV).
"""
import math, sys
import numpy as np
from scipy.integrate import quad
fails = 0
def check(ok, what):
    global fails
    print("   %-66s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)

KPC = 3.0857e21
RS, RHO_SUN, RSUN, RMAX = 24.0, 0.4, 8.5, 260.0
def rho(r):
    x = r/RS; x0 = RSUN/RS
    return RHO_SUN*(x0*(1+x0)**2)/(x*(1+x)**2)
def D(psi):
    """Line-of-sight column (GeV cm^-3 kpc) at angle psi from the Galactic centre."""
    c = math.cos(psi)
    smax = RSUN*c + math.sqrt(max(RMAX**2 - RSUN**2*(1-c*c), 0))
    return quad(lambda s: rho(math.sqrt(RSUN**2 + s*s - 2*RSUN*s*c)), 0, smax, limit=400, points=[RSUN*c] if c > 0 else None)[0]

def rho_burkert(r, rb=9.0):
    """a cored halo with the same local density"""
    f = lambda x: 1/((1 + x)*(1 + x*x))
    return RHO_SUN*f(r/rb)/f(RSUN/rb)
def D_profile(psi, prof):
    c = math.cos(psi)
    smax = RSUN*c + math.sqrt(max(RMAX**2 - RSUN**2*(1-c*c), 0))
    return quad(lambda s: prof(math.sqrt(RSUN**2 + s*s - 2*RSUN*s*c)), 0, smax, limit=400, points=[RSUN*c] if c > 0 else None)[0]

def interp_area(E, table):
    Es = np.log([1e8, 2.2e8, 1e9]); As = np.log(table)
    return math.exp(np.interp(math.log(E), Es, As))*1e4   # cm^2

if __name__ == "__main__":
    psis = np.linspace(1e-4, math.pi, 4001)
    Ds = np.array([D(p) for p in psis])
    Dbar = 0.5*np.trapezoid(Ds*np.sin(psis), psis)      # sky average
    print("sky-averaged column %.2f GeV cm^-3 kpc; toward the centre at 1 deg %.1f" % (Dbar, D(math.radians(1))))
    print("   that is %.2e GeV cm^-2" % (Dbar*KPC))
    l, b = math.radians(216.06), math.radians(-11.13)
    psi_km3 = math.acos(math.cos(l)*math.cos(b))
    print("KM3-230213A: %.1f deg from the Galactic centre; column ratio to the sky mean %.2f" % (math.degrees(psi_km3), D(psi_km3)/Dbar))
    check(D(psi_km3)/Dbar < 1, "the event's direction carries less than the sky-average Galactic column")
    J = 4*math.pi*Dbar*KPC                               # GeV cm^-2 sr, all sky
    # A Z2-breaking Yukawa y N1 L H gives, for M >> m_W (Goldstone equivalence over the four real
    # components of H), BR(l W) = 1/2, BR(nu Z) = 1/4, BR(nu h) = 1/4: half the decays carry a primary
    # neutrino at M/2, and 0.231 of those survive electroweak showering at x = 1 (HDMSpectra).
    BR_NU = 0.5
    frac_line = 0.231*BR_NU
    # the effective areas are summed over flavours (IceCube data release; Gen2 files summed here), so
    # they multiply the flux per flavour: after oscillation over Galactic distances about a third of
    # the line neutrinos are in each flavour
    PER_FLAVOUR = 1.0/3.0
    print("\n=== Galactic line events, sky-integrated (%.3f line neutrinos per decay) ===" % frac_line)
    for M in [3.0e8, 4.4e8, 4.916e8, 6.0e8]:
        E = M/2
        for tau in [1e29, 1e30]:
            flux = PER_FLAVOUR*frac_line*J/(4*math.pi*M*tau)   # cm^-2 s^-1 per flavour (all sky)
            n_ic = flux*interp_area(E, [1.01e3, 1.85e3, 5.13e3])*12.6*3.156e7
            n_g2 = flux*interp_area(E, [5.5e3, 1.07e4, 2.48e4])*10*3.156e7
            print("   M = %.0f PeV, tau = %.2e s: IceCube 12.6 yr %.3f events, Gen2 optical 10 yr %.2f" % (M/1e6, tau, n_ic, n_g2))
    # at the ceiling: N = f D A T/(3 M tau), f line neutrinos per decay, D the sky-averaged column, A the all-flavour
    # sky-averaged effective area at M/2, a third of the line in each flavour
    M = 4.916e8
    print("   at M/2 = %.1f PeV the all-flavour effective areas are %.2e m^2 (IceCube EHE) and %.2e m^2 (Gen2 optical)"
          % (M/2e6, interp_area(M/2, [1.01e3, 1.85e3, 5.13e3])/1e4, interp_area(M/2, [5.5e3, 1.07e4, 2.48e4])/1e4))
    n_formula = frac_line*Dbar*KPC*interp_area(M/2, [1.01e3, 1.85e3, 5.13e3])*12.6*3.156e7/(3*M*1e29)
    n_loop = PER_FLAVOUR*frac_line*J/(4*math.pi*M*1e29)*interp_area(M/2, [1.01e3, 1.85e3, 5.13e3])*12.6*3.156e7
    check(abs(n_formula/n_loop - 1) < 1e-12, "N = f D A T/(3 M tau) reproduces the IceCube count")
    # the lifetime at which Gen2 expects one line event
    n_g2_29 = PER_FLAVOUR*frac_line*J/(4*math.pi*M*1e29)*interp_area(M/2, [5.5e3, 1.07e4, 2.48e4])*10*3.156e7
    print("   at the ceiling, Gen2 expects one line event in ten years for tau = %.1e s" % (1e29*n_g2_29))
    # the halo: a cored Burkert profile with the same local density
    Db = np.array([D_profile(p, rho_burkert) for p in psis])
    Dbar_b = 0.5*np.trapezoid(Db*np.sin(psis), psis)
    print("   a Burkert halo (r_B = 9 kpc) with the same local density has sky-averaged column %.2f, %.2f of the NFW one"
          % (Dbar_b, Dbar_b/Dbar))
    # the Yukawa that allows the decay also gives the lightest neutrino a seesaw mass y^2 v^2/M. With
    # Gamma = y^2 M/(8 pi) for a Majorana N1 summed over channels and tau above 1e29 s:
    HBAR = 6.582e-25                                      # GeV s
    for M in (4.916e8,):
        tau = 1e29/HBAR                                    # GeV^-1
        y2 = 8*math.pi/(M*tau)
        m1 = y2*174.0**2/M*1e9                             # eV
        print("\n   tau > 1e29 s at M = %.1f PeV needs y < %.1e; the lightest neutrino then gets y^2 v^2/M < %.0e eV"
              % (M/1e6, math.sqrt(y2), m1))
    # plant: halving the lifetime doubles the events
    f1 = frac_line*J/(4*math.pi*4.4e8*1e30); f2 = frac_line*J/(4*math.pi*4.4e8*5e29)
    check(abs(f2/f1 - 2) < 1e-12, "plant: halving the lifetime doubles the flux")
    check(abs(Dbar - 7.15) < 0.6, "sky-averaged column near the 7.15 GeV cm^-3 kpc referees reproduced")
    sys.exit(1 if fails else 0)
