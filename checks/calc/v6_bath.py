#!/usr/bin/env python3
"""Separate Ways v6, section 5, and its inflation companion, sections 2 to 4: dark matter made by the plasma
rather than the expansion.

Graviton exchange between Standard Model particles makes N1 pairs at a rate density
R(T) = alpha T^8 / M_P^4 (reduced Planck mass). In radiation domination with constant g*,
dY/dT = -R/(s H T) integrates from the highest temperature T_max to
    Y_grav = alpha (T_max/M_P)^3 / C,   C = 3 (2 pi^2/45) g* sqrt(pi^2 g*/90).
The yield the dark matter needs at mass M is Y_obs = (rho_DM/s0)/M. ALPHA is the fermion
coefficient from the literature (see the v6 notes for the source); every result below scales as
alpha^(-1/3) in temperature and alpha in yield.
"""
import math, sys
fails = 0
def check(ok, what):
    global fails
    print("   %-66s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)

MP = 2.435e18                    # reduced Planck mass, GeV
M_PL = 1.22089e19
GSTAR = 106.75
ALPHA = 1.2e-2                   # Majorana N: Clery et al. 2022 amplitudes, 45 Weyl = 22.5 Dirac counting (1.7e-2 if 45 Dirac)
RHO_OVER_S = 9.7e-48/2.2215e-38  # GeV, the release's inputs
M_C = 4.916e8
C = 3*(2*math.pi**2/45)*GSTAR*math.sqrt(math.pi**2*GSTAR/90)
CH = math.sqrt(math.pi**2*GSTAR/90)

def Y_grav(Tmax, alpha=ALPHA): return alpha*(Tmax/MP)**3/C
def T_star(M, alpha=ALPHA):
    """Highest bath temperature at which graviton exchange alone still makes no more than all the dark matter."""
    return MP*((RHO_OVER_S/M)*C/alpha)**(1/3)

if __name__ == "__main__":
    print("C = %.1f;  Y_obs(M_c) = %.3e" % (C, RHO_OVER_S/M_C))
    Ts = T_star(M_C)
    print("\n=== graviton exchange in the bath ===")
    print("   T* (alone matches the abundance at M_c) = %.3e GeV" % Ts)
    print("   T at H = M_c: %.3e GeV;  T*/T(H=M_c) = %.2f" % (math.sqrt(M_C*MP/CH), Ts/math.sqrt(M_C*MP/CH)))
    for T in [1e15, 1e16, 1e17, MP]:
        print("   T_max = %.2e GeV: graviton exchange alone gives %.2e times the dark matter at M_c" % (T, Y_grav(T)*M_C/RHO_OVER_S))
    # the mass the abundance would then imply if this channel dominated
    for T in [1e15, 1e16, 1e17]:
        print("   if it alone set the abundance at T_max = %.0e GeV: M = %.3g GeV" % (T, RHO_OVER_S/Y_grav(T)))
    # does N1 reach equilibrium with the bath through gravity? The number-changing rate per N1 is
    # Gamma = R/n_eq with n_eq = (3/4)(zeta3/pi^2) 2 T^3 for a Majorana fermion
    print("\n=== equilibration of N1 by graviton exchange: Gamma/H = R/(n_eq H) ===")
    neq = lambda T: 0.75*(1.2020569/math.pi**2)*2*T**3
    for T in [1e16, 1e17, 1e18, MP]:
        print("   T = %.2e: Gamma/H = %.1e" % (T, ALPHA*T**8/MP**4/neq(T)/(CH*T*T/MP)))
    print("   Gamma/H = 1 at T = %.1f M_P" % ((neq(1.0)*CH/ALPHA)**(1/3)))
    print("\n=== inflation dilutes what the pre-inflation bath made ===")
    Y_planck = Y_grav(MP)
    N_need = math.log(Y_planck*M_C/RHO_OVER_S/1e-3)/3
    print("   pre-inflation yield up to T = M_P: %.2e (x%.1e the dark matter)" % (Y_planck, Y_planck*M_C/RHO_OVER_S))
    print("   e-folds to bring it below 1e-3 of the dark matter: N > %.1f" % N_need)
    # if gravity brings N1 to equilibrium near the Planck temperature (Gamma/H ~ 0.26 at 1e18 GeV above),
    # the pre-inflation yield is the equilibrium one for 2 Majorana states, (3/4)(zeta3/pi^2) 2 / ((2 pi^2/45) g*)
    Y_eq = 0.75*(1.2020569/math.pi**2)*2/((2*math.pi**2/45)*GSTAR)
    print("   if N1 reached equilibrium (Y_eq = %.2e): N > %.1f" % (Y_eq, math.log(Y_eq*M_C/RHO_OVER_S/1e-3)/3))
    # Z2-even portal: L = -(c/(2 Lam)) H^+H N^c-bar N. Each of the four real Higgs scalars
    # annihilates through the contact vertex with <sigma v> ~ c^2/(16 pi Lam^2) (Majorana pair,
    # massless limit), so R ~ 4 (zeta(3)/pi^2)^2 T^6 c^2/(16 pi Lam^2). Then R/(s H T) is
    # independent of T and Y = r c^2 T_rh M_P / Lam^2, with r below. Factors of order one from the
    # operator's normalisation are not tracked.
    zeta3 = 1.2020569
    Rcoef = 4*(zeta3/math.pi**2)**2/(16*math.pi)
    r = Rcoef/((2*math.pi**2/45)*GSTAR*CH)
    print("\n=== Z2-even portal ===")
    print("   Y_portal = %.1e c^2 T_rh M_P/Lam^2" % r)
    # without inflation, run to the reduced Planck temperature: the c at which it alone makes the dark matter,
    # and the size the Majorana mass suggests (the operator breaks the same N1-number symmetry as M N1 N1)
    print("   plasma run to M_P (Lam = M_P): c = %.1e makes all the dark matter; c ~ M/M_P = %.1e"
          % (math.sqrt((RHO_OVER_S/M_C)/(r*MP/MP)), M_C/MP))
    for T, M in [(1.87e13, 4.494e8), (5e13, 3.042e8), (1e14, 2.202e8)]:     # perturbative reheating, Table 1 of the inflation companion
        Yobs = RHO_OVER_S/M
        cmax = math.sqrt(0.1*Yobs/(r*T/MP))
        print("   T_rh = %.2e GeV, M = %.0f PeV: c < %.1e (Lam = M_P) keeps it under 10%% of the dark matter" % (T, M/1e6, cmax))
    # inflaton condensate, phi phi -> graviton -> N N, helicity-suppressed (m_N/m_phi)^2.
    # Clery, Mambrini, Olive, Verner, PRD 105 (2022) 075005, Eq. (50), k = 2 (quadratic minimum),
    # Dirac counting (a Majorana N is at most this); rho_end ~ 3 H_inf^2 M_P^2.
    def omega_condensate(mX, m_phi, H_inf, T_rh, g=GSTAR):
        lam = m_phi**2/(2*MP**2)
        rho_end = 3*H_inf**2*MP**2
        rho_rh = math.pi**2*g/30*T_rh**4
        return 0.1*((1/64)/2.4**4)*2*(1e-11/lam)*(1e40/rho_rh)**(-0.25)*(rho_end/1e64)**0.5*(mX/8.3e9)**3
    print("\n=== inflaton condensate (Clery et al. Eq. 50, k = 2), Omega h^2 against the observed 0.12 ===")
    for m_phi, H_inf in [(1e13, 1e12), (3e13, 1e12), (1e13, 1e13), (3e13, 1e13)]:
        for T in [1.87e13, 5e13]:
            if 3*H_inf**2*MP**2 < math.pi**2*GSTAR/30*T**4: continue
            om = omega_condensate(4.2e8, m_phi, H_inf, T)
            print("   m_phi = %.0e, H_inf = %.0e, T_rh = %.2e GeV: Omega h^2 = %.1e (%.3f%% of the dark matter)" % (m_phi, H_inf, T, om, 100*om/0.12))
    # a direct inflaton coupling phi N1 N1, which the Z2 also allows for a singlet inflaton: each decay to
    # N1 pairs gives two N1, and n_phi/s = (3/4) T_rh/m_phi when the inflaton's energy becomes radiation at T_rh
    print("\n=== direct inflaton decay to N1 pairs ===")
    # a Yukawa of the size the Majorana mass suggests, y ~ M/M_P, with L = -(y/2) phi N1c-bar N1, so
    # Gamma(phi -> N1 N1) = y^2 m_phi/(16 pi) and Gamma = H(T_rh); Y_N = 2 B (3/4) T_rh/m_phi = 3 y^2 M_P/(32 pi CH T_rh)
    for T, M in [(1.87e13, 4.494e8), (5e13, 3.042e8)]:
        y = M/MP
        YN = 3*y*y*MP/(32*math.pi*CH*T)
        print("   y = M/M_P = %.1e at T_rh = %.2e GeV: inflaton decays alone make %.0f times the dark matter" % (y, T, YN*M/RHO_OVER_S))
        print("      so y below %.2f M/M_P keeps it under a tenth of the dark matter" % math.sqrt(0.1/(YN*M/RHO_OVER_S)))
    for m_phi, T, M in [(1e13, 1.87e13, 4.494e8), (3e13, 5e13, 3.042e8)]:     # perturbative quadratic masses (Table 1 of the inflation companion)
        Bmax = 0.1*(RHO_OVER_S/M)/(2*0.75*T/m_phi)
        print("   m_phi = %.0e, T_rh = %.2e GeV: branching to N1 pairs below %.1e keeps it under 10%% of the dark matter" % (m_phi, T, Bmax))
    # if inflation also makes the observed perturbations: P_t = 2 H^2/(pi^2 M_P^2) and r = P_t/A_s < 0.036
    # (BICEP/Keck 2021), A_s = 2.1e-9 (Planck 2018)
    H_max = math.pi*MP*math.sqrt(0.036*2.1e-9/2)
    r_b = 2*1e24/(math.pi**2*MP**2*2.1e-9); eps_b = 1e24/(8*math.pi**2*MP**2*2.1e-9)
    print("\n=== at H_inf = 1e12 GeV, if inflation makes the perturbations: r = %.1e, slow-roll epsilon = %.1e ===" % (r_b, eps_b))
    # the range of masses after inflation, H_inf = 1e12 GeV: instant reheating gives the hottest plasma, and
    # the plasma then sets the mass; the heaviest is set by M < H_inf, reached as M ~ T_rh^(-1/2) below 3e12 GeV
    T_inst = (30*3e24*MP**2/(math.pi**2*GSTAR))**0.25
    M_low = RHO_OVER_S/Y_grav(T_inst)
    T_top = 3e12*(1.1116e9/1e12)**2
    print("   instant reheating at H_inf = 1e12 GeV: T_rh = %.1e GeV, the plasma sets M = %.1e GeV" % (T_inst, M_low))
    print("   M reaches H_inf = 1e12 GeV at T_rh = %.1e GeV (M ~ T_rh^-1/2 from 1112 PeV at 3e12 GeV)" % T_top)
    # an inflaton that reheats through a Yukawa coupling to fermions, Gamma = y^2 m_phi/(8 pi) = H(T_rh), gets a
    # thermal mass y T/sqrt(6) from them while they are light; compare with H_inf when the plasma meets the plateau
    y_rh = math.sqrt(8*math.pi*CH*1.87e13**2/MP/1e13)
    print("   reheating to 1.87e13 GeV through a Yukawa (m_phi = 1e13 GeV) needs y = %.3f; at T = %.1e GeV its thermal mass is %.1e GeV, %.0f H_inf"
          % (y_rh, T_inst, y_rh*T_inst/math.sqrt(6), y_rh*T_inst/math.sqrt(6)/1e12))
    # the thermal mass y T/sqrt6 exceeds the Hubble rate of the radiation era before inflation below this T, so the
    # escape needs the fermions heavier than that T on the plateau, y phi above it; their one-loop energy is then at
    # least (y phi)^4/(16 pi^2), against the plateau's 3 H_inf^2 M_P^2
    T_x = y_rh*MP/(math.sqrt(6)*CH)
    V_plat = 3e24*MP**2
    print("   that thermal mass exceeds H below T = %.1e GeV; the escape needs y phi above it, one-loop energy %.1e GeV^4 against the plateau's %.1e"
          % (T_x, T_x**4/(16*math.pi**2), V_plat))
    # a trilinear coupling mu phi H^dag H (4 real scalars, Gamma = mu^2/(8 pi m_phi)) bends the plateau at one loop by
    # about mu^2 phi/(8 pi^2); with phi = M_P and the slope sqrt(2 eps) V/M_P that A_s allows at H_inf = 1e12 GeV
    eps = 1e24/(8*math.pi**2*MP**2*2.1e-9)
    slope = math.sqrt(2*eps)*V_plat/MP
    mu2 = 8*math.pi**2*slope/MP
    T_tri = math.sqrt(mu2/(8*math.pi*1e13)*MP/CH)
    print("   a trilinear coupling that keeps its one-loop slope below the plateau's (phi = M_P, m_phi = 1e13 GeV) reheats to at most %.0e GeV" % T_tri)
    # the inflaton-N1 bound does not depend on m_phi: y_max = 0.05 M/M_P at the reference point, scaling as sqrt(T_rh/M)
    y_ref = math.sqrt(0.1*RHO_OVER_S/4.494e8*32*math.pi*CH*1.87e13/(3*MP))
    print("   y_max = %.1e (T_rh/1.87e13 GeV)^(1/2) (449 PeV/M)^(1/2), whatever m_phi is" % y_ref)
    # the least population that changes production by ten per cent, against the plasma, at N = 58: the energy goes as
    # e^(4N) from its threshold, 13.1 e-folds at H_inf = 1e12, T_rh = 1.87e13 GeV (Table A.1 of the inflation companion)
    print("   at 58 e-folds that population carries %.0e times the plasma's energy at the bang" % math.exp(4*(58 - 13.1)))
    print("\n=== tensor bound on the inflationary Hubble rate (r < 0.036, A_s = 2.1e-9): H_inf < %.1e GeV ===" % H_max)
    # e-folds before the end of inflation at which today's Hubble scale left the horizon:
    # N = ln(a_end H_inf/(a0 H0)), with matter-like reheating a_end/a_rh = (rho_rh/rho_end)^(1/3)
    # and entropy conservation after reheating, a_rh/a0 = (T0/T_rh)(g*s0/g*s)^(1/3)
    H0, T0, GS0 = 1.44e-42, 2.3487e-13, 3.909
    for H_inf, T_rh in [(1e12, 1.87e13), (1e13, 5e13)]:
        rho_end, rho_rh = 3*H_inf**2*MP**2, math.pi**2/30*GSTAR*T_rh**4
        a_end_a0 = (rho_rh/rho_end)**(1/3)*(T0/T_rh)*(GS0/GSTAR)**(1/3)
        print("   H_inf %.0e, T_rh %.2e GeV: today's Hubble scale left the horizon %.1f e-folds before the end of inflation"
              % (H_inf, T_rh, math.log(a_end_a0*H_inf/H0)))
    # Garcia et al. 2017 define T_RH^4 = (40/(pi^2 g)) Gamma^2 M_P^2; here Gamma = H(T_rh) gives
    # T_rh^4 = (90/(pi^2 g)) Gamma^2 M_P^2, so for the same Gamma their T_RH is (40/90)^(1/4) of ours
    print("\n=== reheating-temperature conventions: Garcia et al. T_RH / ours = %.3f for the same Gamma ===" % (40/90)**0.25)
    check(abs(Y_grav(Ts)*M_C/RHO_OVER_S - 1) < 1e-12, "plant: at T* the channel makes exactly the dark matter")
    check(abs(T_star(M_C, 8*ALPHA)/Ts - 0.5) < 1e-12, "plant: eight times the rate halves T*")
    sys.exit(1 if fails else 0)
