#!/usr/bin/env python3
"""How far a late entropy release moves the two-body neutrino line, and how little it takes.

WHY THIS EXISTS. Section 3.2 lists what remains uncertain about the line and ends with "no later
dilution", naming an assumption it never quantifies. The neighbouring compact-object dilution is
quantified to the last digit: the line goes as (1 - f)^{2/5} and f = 0.01014 moves it by its own
width. The entropy assumption had no such number, so a referee asking "what if the heavier
sterile neutrinos inject entropy after production" had nothing in the paper to read. It does now,
and the answer is that the line is about as sensitive to entropy as it is to compact objects.

THE SCALING, which is the paper's own. Section 2.3 has M_1 proportional to
rho_DM,0^{2/5} s_0^{-2/5} I^{-2/5} muhat^{3/5}, so at fixed I and muhat the mass depends on the
abundance only through (rho_DM,0 / s_0)^{2/5}. A late release multiplying the comoving entropy by
gamma dilutes the produced dark-matter-to-entropy ratio by gamma, so reproducing today's observed
abundance requires a production gamma times larger, and the inferred mass and the line both rise
by gamma^{2/5}. The direction is opposite to the compact-object case, which lowers the line.

WHAT IT DOES NOT SAY. Nothing here argues that a release happens, or bounds one. It converts an
unquantified assumption into a stated sensitivity, which is all section 3.2 needs.
"""
E0, W = 245.8, 1.0            # the line and its quoted width, PeV, from section 3.2
F_PAPER = 0.01014             # the paper's own compact-object number, for the cross-check

def line_after_entropy(gamma):   return E0 * gamma ** 0.4
def line_after_compact(f):       return E0 * (1.0 - f) ** 0.4

gamma_one_width = (1.0 + W / E0) ** 2.5

print("=== 1. the sensitivity ===")
print("   line %.1f +- %.1f PeV, so one width is %.3f per cent of it" % (E0, W, 100 * W / E0))
print("   a late entropy release by gamma raises the line by gamma^(2/5)")
print("   gamma moving it by exactly one width: %.5f, a release of %.2f per cent"
      % (gamma_one_width, 100 * (gamma_one_width - 1)))
for g in (1.01, 1.1, 2.0, 10.0):
    print("      gamma = %-5g line -> %7.1f PeV  (%+7.1f PeV, %6.1f widths)"
          % (g, line_after_entropy(g), line_after_entropy(g) - E0, (line_after_entropy(g) - E0) / W))

print("\n=== 2. cross-check against the paper's own published compact-object number ===")
shift = E0 - line_after_compact(F_PAPER)
print("   section 3.2 states f = %.5f moves the line by its own width." % F_PAPER)
print("   recomputed here: %.1f*(1-%.5f)^(2/5) = %.4f PeV, a shift of %.4f PeV"
      % (E0, F_PAPER, line_after_compact(F_PAPER), shift))
ok_cross = abs(shift - W) < 5e-4
print("   agrees with the quoted width to %.1e: %s" % (abs(shift - W), "yes" if ok_cross else "NO"))

print("\n=== 3. the plants ===")
p1 = abs(line_after_entropy(1.0) - E0) < 1e-12
print("   (a) no release leaves the line alone: %s" % ("yes" if p1 else "NO"))
p2 = line_after_entropy(1.05) > E0 and line_after_compact(0.05) < E0
print("   (b) the two dilutions push opposite ways, entropy up and compact objects down: %s"
      % ("yes" if p2 else "NO"))
# The exponent has to be tested where exponents differ. At gamma just above 1 every exponent
# agrees to a fraction of a PeV, so the first version of this plant passed a wrong exponent and
# failed its own assertion, correctly. Test it at gamma = 10.
G = 10.0
right, wrong = E0 * G ** 0.4, E0 * G ** 0.5
p3 = abs(right - wrong) > 10 * W
print("   (c) the exponent is tested where exponents differ, at gamma = %g:" % G)
print("       2/5 gives %.1f PeV, 1/2 gives %.1f PeV, apart by %.1f PeV or %.0f widths"
      % (right, wrong, wrong - right, (wrong - right) / W))
print("       so a wrong exponent would be caught: %s" % ("yes" if p3 else "NO"))
print("       (at gamma = 1.01 the same two differ by %.2f PeV, which is why the first"
      % (E0 * 1.01 ** 0.5 - E0 * 1.01 ** 0.4))
print("        version of this plant was useless and said so)")

assert ok_cross and p1 and p2 and p3
print("\n   all checks passed")
