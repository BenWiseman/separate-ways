#!/usr/bin/env python3
"""The inflation companion of Separate Ways v6, section 4.1: does a trilinear inflaton-Higgs coupling survive the hot era before inflation?

Section 3 of that paper assumes an inflaton that stays on its plateau through the radiation era before inflation and
still reheats afterwards. v6_bath.py shows that a Yukawa coupling to fermions strong enough to reheat to
1.87e13 GeV gives the inflaton a thermal mass 12 times H_inf when the plasma meets the plateau, and caps a
trilinear coupling mu phi H^dag H at the value whose one-loop slope stays below the plateau's. This script
follows the inflaton through the hot era at that largest mu.

On the plateau the four real Higgs scalars have mass^2 mu phi, and the plasma pushes the inflaton with
    F = d/dphi [4 T^4/(2 pi^2) J_B(mu phi/T^2)] = (mu T^2/pi^2) K(y),  K(y) = int x^2 dx/(E (e^E - 1)),
E = sqrt(x^2 + y^2), y^2 = mu phi/T^2, so F = mu T^2/6 while the Higgs is light. In e-folds N, with
T = T0 e^-N and v = dphi/dN,
    v' = -(3 - eps_H) v - F/H^2,   eps_H = (4 rho_r/3 + H^2 v^2)/(2 M_P^2 H^2),
from T0 = M_P with the field at rest, through the moment the plasma meets the plateau, until the plasma
no longer matters. The Higgs's own thermal mass is left out, which overstates the push.
Plants: (1) a massless Higgs in pure radiation must give the attractor displacement -mu M_P^2/(36 C_H^2 T^2);
(2) the Yukawa thermal mass y T/sqrt(6), y = 0.035, put through the same integrator, must move the field by
more than half of M_P, the failure v6_bath.py describes.
"""
import math, sys
from scipy.integrate import quad, solve_ivp
from v6_bath import MP, GSTAR, CH

fails = 0
def check(ok, what):
    global fails
    print("   %-78s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)

H_INF = 1e12
V_PLAT = 3*H_INF**2*MP**2
EPS = H_INF**2/(8*math.pi**2*MP**2*2.1e-9)              # as in v6_bath.py, A_s = 2.1e-9
SLOPE = math.sqrt(2*EPS)*V_PLAT/MP
MU = math.sqrt(8*math.pi**2*SLOPE/MP)                    # one-loop slope mu^2 phi/(8 pi^2) at phi = M_P equals the plateau's
T_MEET = (30*V_PLAT/(math.pi**2*GSTAR))**0.25
RHO = lambda T: math.pi**2/30*GSTAR*T**4

def K(y):
    if y < 1e-6: return math.pi**2/6
    if y > 600: return 0.0                                # e^-600: the Higgs is gone from the plasma
    def f(x):
        E = math.sqrt(x*x + y*y)
        return x*x*math.exp(-E)/(E*-math.expm1(-E))
    return quad(f, 0, math.inf, limit=200)[0]

def F_tri(phi, T, massless=False):
    return MU*T*T/6 if massless else MU*T*T*K(math.sqrt(max(MU*phi, 0.0))/T)/math.pi**2

def follow(force, T_end, plateau=True, phi0=MP, T0=MP):
    V = V_PLAT if plateau else 0.0
    # the state is the displacement d = phi - phi0, so that the tolerances apply to it and not to phi0
    def rhs(N, s):
        d, v = s
        T = T0*math.exp(-N)
        H2 = (RHO(T) + V)/(3*MP**2 - v*v/2)
        eps = (4*RHO(T)/3 + H2*v*v)/(2*MP**2*H2)
        return [v, -(3 - eps)*v - force(phi0 + d, T)/H2]
    N_end = math.log(T0/T_end)
    return solve_ivp(rhs, [0, N_end], [0.0, 0.0], method="LSODA", rtol=1e-10, atol=1.0, dense_output=True)

if __name__ == "__main__":
    print("plateau V = %.2e GeV^4, slope %.2e GeV^3; largest mu = %.2e GeV; Higgs mass on the plateau sqrt(mu M_P) = %.1e GeV"
          % (V_PLAT, SLOPE, MU, math.sqrt(MU*MP)))
    print("the plasma meets the plateau at T = %.2e GeV" % T_MEET)

    print("\n=== plant 1: massless Higgs, pure radiation, against the attractor ===")
    T1 = 1e15
    sol = follow(lambda p, T: F_tri(p, T, massless=True), T1, plateau=False)
    d_num = sol.y[0][-1]
    d_an = -MU*MP**2/(36*CH**2*T1**2)
    print("   displacement at T = %.0e GeV: integrated %.4e GeV, attractor %.4e GeV" % (T1, d_num, d_an))
    check(abs(d_num/d_an - 1) < 1e-3, "plant 1: the integrator reproduces the radiation-era attractor to 0.1 per cent")

    print("\n=== plant 2: a Yukawa thermal mass y T/sqrt6, y = 0.035, through the same integrator ===")
    y = 0.035
    sol = follow(lambda p, T: (y*T)**2/6*p, T_MEET)
    print("   field when the plasma meets the plateau: %.3f M_P (started at 1)" % (1 + sol.y[0][-1]/MP))
    check(abs(sol.y[0][-1]) > 0.5*MP, "plant 2: a thermal mass above H moves the field by more than M_P/2")

    print("\n=== the trilinear coupling at its largest mu, phi = M_P at the bang ===")
    T_end = T_MEET*1e-3
    sol = follow(lambda p, T: F_tri(p, T), T_end)
    N_meet = math.log(MP/T_MEET)
    d_m, v_m = sol.sol(N_meet)
    H2m = (RHO(T_MEET) + V_PLAT)/(3*MP**2 - v_m**2/2)
    kin = H2m*v_m**2/2
    d_meet = d_m/MP
    d_end = sol.y[0][-1]/MP
    # the push's curvature against H^2 along the way, by finite difference in phi
    curv = 0.0
    for N in [i*sol.t[-1]/400 for i in range(401)]:
        T = MP*math.exp(-N); p = MP + sol.sol(N)[0]
        H2 = (RHO(T) + V_PLAT)/(3*MP**2)
        dF = (F_tri(p*1.001, T) - F_tri(p*0.999, T))/(0.002*p)
        curv = max(curv, abs(dF)/H2)
    print("   when the plasma meets the plateau: displacement %.2e M_P, kinetic energy %.1e of the plateau's" % (d_meet, kin/V_PLAT))
    print("   once the plasma no longer matters (T = %.1e GeV): displacement %.2e M_P (%.1f per cent)" % (T_end, d_end, 100*abs(d_end)))
    print("   largest |dF/dphi|/H^2 on the way: %.1e" % curv)
    check(abs(d_end) < 0.01, "the plasma moves the inflaton by less than 1 per cent of M_P")
    check(kin/V_PLAT < 1e-3, "and leaves it a kinetic energy below a thousandth of the plateau's")
    check(curv < 1e-2, "and gives it no curvature comparable with H^2")
    print("\n%s" % ("all checks pass" if fails == 0 else "%d FAILED" % fails))
    sys.exit(1 if fails else 0)
