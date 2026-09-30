#!/usr/bin/env python3
"""What the neutrino-mass floor does to structure growth, which is nothing differential.

WHY THIS EXISTS. The fold forces a neutrino-mass sum at or above 58.8 meV and forbids evolving
dark energy. Neutrino mass suppresses the growth of structure, so a cosmology referee will ask
what the floor does to sigma_8 and to the S_8 tension, and neither manuscript mentions sigma_8,
S_8 or growth anywhere. The answer is short and it is a null, and a null that is computed is
worth more here than a silence.

THE ANSWER. Free-streaming neutrinos suppress the small-scale power by about 8 f_nu with
f_nu = Sigma m_nu / (93.14 h^2 Omega_m) the neutrino fraction of the matter budget. At the floor
that is a few per cent. It is not a signature of the fold, because the identical suppression
appears in LCDM at the identical mass: the fold fixes the VALUE of Sigma m_nu, it does not change
what a given Sigma m_nu does. And the usual baseline already carries neutrino mass, since Planck
quotes S_8 for a fit assuming 0.06 eV, so the comparison is 58.8 meV against 60 meV and not
against zero.

THE TRAP, recorded because it is seductive and wrong. In the rejected local-clock branch the
model sits about 3.7 per cent above LCDM in S_8, and a sum near 61.8 meV would cancel exactly
that, which is close enough to 58.8 meV to look like a prediction. It is an artefact of comparing
a model WITH neutrinos against a baseline WITHOUT them. Panel A of
tangents/h61_lowell/s8_neutrino_door.py shows the effect is common, not differential.
"""
SUM_FLOOR = 0.0588        # eV, the exact-stabilisation floor, normal ordering, m_1 = 0
SUM_PLANCK = 0.060        # eV, the value the standard S_8 baseline already assumes
OMEGA_M, H = 0.315, 0.674

def f_nu(sigma_m):  return sigma_m / (93.14 * H * H * OMEGA_M)
def dP_P(sigma_m):  return -8.0 * f_nu(sigma_m)          # small-scale power suppression
def dS8(sigma_m):   return 0.5 * dP_P(sigma_m)           # S_8 goes as sqrt of the power

print("=== 1. the floor's own suppression, against a massless baseline ===")
for s in (SUM_FLOOR, SUM_PLANCK, 0.10, 0.20):
    print("   Sigma m_nu = %.4f eV   f_nu = %.5f   dP/P = %+.2f %%   dS_8 = %+.2f %%"
          % (s, f_nu(s), 100*dP_P(s), 100*dS8(s)))

print("\n=== 2. against the baseline that is actually used, which already has mass ===")
d = dS8(SUM_FLOOR) - dS8(SUM_PLANCK)
print("   the standard S_8 baseline assumes %.3f eV, the fold requires %.4f eV" % (SUM_PLANCK, SUM_FLOOR))
print("   difference in S_8 between them: %+.3f per cent" % (100*d))
print("   so against the baseline anyone actually quotes, the fold moves S_8 by well under a per cent.")

print("\n=== 3. why it is not a signature ===")
print("   the same Sigma m_nu in LCDM gives the same %+.2f per cent." % (100*dS8(SUM_FLOOR)))
print("   the fold fixes the value of Sigma m_nu; it does not change what a given value does,")
print("   so no part of this distinguishes the fold from LCDM at the same neutrino mass.")

print("\n=== 4. the trap, stated structurally because the number is not the point ===")
print("   In the rejected local-clock branch the model sits above LCDM in S_8, and one can solve")
print("   for the neutrino sum that would cancel that excess and land near the floor. The")
print("   comparison is invalid whatever it returns, because it sets a model WITH neutrinos")
print("   against a baseline WITHOUT them, and the same neutrinos lower LCDM by the same amount.")
print("   The tangent that found it, h61_lowell/s8_neutrino_door.py, labels its own number an")
print("   artefact in capitals. Recomputing it here returns a different value anyway, because")
print("   that line divides an S_8 excess by a power suppression, which is a factor of two:")
print("     against the S_8 suppression  : %.4f eV" % (SUM_FLOOR*0.0368/abs(dS8(SUM_FLOOR))))
print("     against the power suppression: %.4f eV  <- the tangent's 0.0618" % (SUM_FLOOR*0.0368/abs(dP_P(SUM_FLOOR))))
print("   Neither is a result. The structural objection is what carries.")

print("\n=== 5. the plants ===")
p1 = dS8(0.0) == 0.0
print("   (a) no neutrino mass, no suppression: %s" % ("yes" if p1 else "NO"))
p2 = dS8(0.2) < dS8(0.1) < 0
print("   (b) more mass suppresses more, and the sign is negative: %s" % ("yes" if p2 else "NO"))
p3 = abs(d) < 0.001
print("   (c) the floor and the standard assumption differ by under 0.1 per cent in S_8,")
print("       which is the whole point of section 2: %s" % ("yes" if p3 else "NO"))
assert p1 and p2 and p3
print("\n   all checks passed")
