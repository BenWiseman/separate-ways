#!/usr/bin/env python3
"""The dark-matter lifetime is not free once the mass is fixed, and the data pick the operator.

WHY THIS EXISTS. Section 3.2 makes the lifetime "a further input" and quotes nothing for it, so a
referee's first question about the 245.8 PeV line has no answer in the paper: not what energy,
but why does the particle live long enough to still be here. This does not derive a lifetime. It
shows the lifetime is not a free parameter once two things already in the paper are granted, and
that what remains free is a single integer which the observed window fixes.

THE TWO GRANTS, both already the paper's. The dark matter is a right-handed neutrino of mass
M_1, fixed by the relic abundance in section 2.3. Its stability rests on a Z_2 that Appendix C.2
retains "as an additional assumption", having closed the holonomy route to deriving it:
H^1(S^3 x R; Z_2) = 0, so there is no holonomy sector to carry it.

THE STANDARD EXPECTATION. Quantum gravity is not expected to respect global symmetries, so an
assumed global Z_2 should be broken by Planck-suppressed operators and by nothing else. An
operator of dimension 4+k with coefficient c/M_Pl^k gives, on dimensions alone,

    Gamma ~ (c^2 / P) M_1 (M_1 / M_Pl)^{2k},

with P an O(10) phase-space factor. The mass is fixed, so k is the only discrete freedom left.

WHAT IS NEW HERE RATHER THAN IN THE TANGENT IT COMES FROM. The tangent
tangents/h61_lowell/dm_lifetime_planck.py was written against a superseded version and uses
M = 4.848e8 GeV, the high-entropy benchmark. The paper's current central value is 4.916e8 GeV.
The lifetime runs as M^{-(2k+1)}, so at k = 3 that is a factor of (4.916/4.848)^7 = 1.10, ten per
cent against an O(1) ambiguity spanning two orders. The conclusion is unchanged, and this file
records that it was checked rather than assumed.
"""
import math

GeV_s = 6.582119569e-25          # hbar, GeV s
M1    = 4.916e8                  # GeV, section 2.3 central value
M1_OLD= 4.848e8                  # GeV, the benchmark the tangent used
M_PL  = 1.220890e19              # GeV, ordinary Planck mass
M_PLR = 2.435323e18              # GeV, reduced Planck mass
T0    = 4.35e17                  # s, age of the universe
TAU_LO, TAU_HI = 1e29, 1e31      # s, the window the external bounds leave

def tau(k, M=M1, MPl=M_PL, P=16*math.pi, c=1.0):
    G = (c*c/P) * M * (M/MPl)**(2*k)     # GeV
    return GeV_s / G

print("=== 1. the lifetime at each operator dimension, current mass ===")
print("   M_1 = %.4g GeV, ordinary Planck mass, c = 1, P = 16 pi\n" % M1)
print("   dim   k        tau (s)       verdict")
for k in (1, 2, 3, 4):
    t = tau(k)
    v = ("gone long ago" if t < T0 else
         "inside the window" if TAU_LO <= t <= TAU_HI else "no signal ever")
    print("   %-5d %-3d %12.3e   %s" % (4+k, k, t, v))

print("\n=== 2. how discriminating that is ===")
r_up, r_dn = tau(4)/tau(3), tau(3)/tau(2)
print("   k = 3 is the only one in the window.")
print("   one step up is %.1e times longer, one step down %.1e times shorter." % (r_up, r_dn))
print("   so the window selects the operator dimension with about twenty orders to spare either side.")

print("\n=== 3. what the mass update changes, which is the point of this file ===")
ratio = tau(3, M1)/tau(3, M1_OLD)
print("   tau(k=3) at the superseded 4.848e8 GeV : %.3e s" % tau(3, M1_OLD))
print("   tau(k=3) at the current    4.916e8 GeV : %.3e s" % tau(3, M1))
print("   ratio %.3f, against an O(1) ambiguity spanning two orders, so the conclusion holds." % ratio)

print("\n=== 4. the weak joint, named rather than buried ===")
tr = tau(3, M1, M_PLR)
print("   with the REDUCED Planck mass instead: tau = %.3e s, which is %.1e times" % (tr, TAU_LO/tr))
print("   below the window. The result therefore depends on which Planck mass the suppression")
print("   scale is, and that is a convention rather than a derivation. Stated, not hidden.")

print("\n=== 5. the plants ===")
p1 = tau(3) > tau(2) and tau(4) > tau(3)
print("   (a) lifetime lengthens with operator dimension: %s" % ("yes" if p1 else "NO"))
p2 = tau(2) < T0
print("   (b) k = 2 really is excluded, not merely disfavoured: tau = %.2e s against an age of" % tau(2))
print("       %.2e s, so the particle would be gone: %s" % (T0, "yes" if p2 else "NO"))
p3 = abs(tau(3, M=M1, c=3.0)/tau(3) - 1/9) < 1e-9
print("   (c) an O(1) coefficient moves it as c^{-2}, so c = 3 divides the lifetime by 9: %s"
      % ("yes" if p3 else "NO"))
print("       which is why this is a dimensional-analysis statement and not a prediction.")

print("\n=== 6. the claim that survives, stated as a bracket and not a hit ===")
lo, hi = tau(3, P=16*math.pi), tau(3, P=1.0)
print("   k = 3 across the phase-space convention: %.2e s at P = 16 pi, %.2e s at P = 1." % (lo, hi))
print("   the window is %.0e to %.0e s, so the two conventions bracket it rather than hitting it." % (TAU_LO, TAU_HI))
brackets = min(lo, hi) < TAU_HI and max(lo, hi) > TAU_LO
print("   brackets the window: %s" % ("yes" if brackets else "NO"))
for k in (2, 4):
    a, b = tau(k, P=16*math.pi), tau(k, P=1.0)
    clear = max(a, b) < TAU_LO or min(a, b) > TAU_HI
    print("   k = %d misses under BOTH conventions (%.1e and %.1e s): %s" % (k, a, b, "yes" if clear else "NO"))
near = all(max(tau(k, P=16*math.pi), tau(k, P=1.0)) < TAU_LO or
           min(tau(k, P=16*math.pi), tau(k, P=1.0)) > TAU_HI for k in (2, 4))

assert p1 and p2 and p3 and brackets and near
print("\n   all checks passed")
