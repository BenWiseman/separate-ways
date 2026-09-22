# A final error hunt found that the inverted-ordering sum quoted in 3.3 and Figure 3,
# 98.9 meV, does not follow from the oscillation inputs the paper states. Check it, because
# a paper whose whole claim is that its numbers follow from its inputs cannot carry one that
# does not.
#
# 2.5 states ONE input set and uses it everywhere the floor's width is quoted:
#     Delta m^2_21 = (7.53 +- 0.18)e-5 eV^2
#     Delta m^2_31 = (2.510 +- 0.030)e-3 eV^2

d21 <- 7.53e-5;  s21 <- 0.18e-5
d31 <- 2.510e-3; s31 <- 0.030e-3

cat("=== 1. normal ordering, m1 = 0, which is the paper's actual result\n\n")
no <- sqrt(d21) + sqrt(d31)
# propagate: d(sqrt(x)) = dx/(2 sqrt(x))
sno <- sqrt( (s21/(2*sqrt(d21)))^2 + (s31/(2*sqrt(d31)))^2 )
cat(sprintf("   Sum m_nu = sqrt(d21) + sqrt(d31) = %.3f +- %.3f meV\n", no*1000, sno*1000))
cat(sprintf("   paper quotes 58.78 +- 0.32 meV -> agrees to %.4f meV\n", abs(no*1000 - 58.78)))
cat(sprintf("   share of the width from d31: %.1f per cent  (paper says 89)\n",
            100*(s31/(2*sqrt(d31)))^2 / sno^2))
stopifnot(abs(no*1000 - 58.78) < 0.01)

cat("\n=== 2. inverted ordering, m3 = 0, ON THE SAME INPUTS\n\n")
m1 <- sqrt(d31); m2 <- sqrt(d31 + d21)
io <- m1 + m2
cat(sprintf("   m1 = sqrt(|d31|)        = %.3f meV\n", m1*1000))
cat(sprintf("   m2 = sqrt(|d31| + d21)  = %.3f meV\n", m2*1000))
cat(sprintf("   Sum m_nu                = %.2f meV\n\n", io*1000))
cat("   The paper and its neutrino figure USED to carry 98.9 meV, and that was NOT this\n")
cat("   number. Both were corrected to 100.9 meV on 2026-09-21; 3.3 and Figure 7 now agree\n")
cat("   with the line above. This file is the record of that correction.\n")

cat("\n=== 3. where 98.9 came from, since it did not come from the stated inputs\n\n")
cat("   fig5_data.py records the mass triple (49.08, 49.84, 0) meV. Invert it:\n\n")
a <- 49.08e-3; b <- 49.84e-3
cat(sprintf("      implied |d31| = m1^2        = %.4e eV^2   (paper states %.4e)\n", a^2, d31))
cat(sprintf("      implied  d21  = m2^2 - m1^2 = %.4e eV^2   (paper states %.4e)\n", b^2-a^2, d21))
cat("\n   So the solar splitting matches and the atmospheric one does not: the figure used\n")
cat(sprintf("   an atmospheric splitting %.1f per cent below the value the paper states and\n",
            100*(1 - a^2/d31)))
cat("   propagates everywhere else. No source in the repository carries that value, so it\n")
cat("   cannot be defended as a deliberate inverted-ordering global fit; it is drift.\n")

cat("\n=== 4. what to do about it\n\n")
cat("   Two defensible options, and only one of them is available without a new citation.\n\n")
cat("   (a) Quote the inverted-ordering global fit's OWN atmospheric splitting, which does\n")
cat("       differ from the normal-ordering one. That is correct practice, and it needs a\n")
cat("       sourced number the paper does not currently carry. Not available today.\n")
cat("   (b) Compute it from the inputs the paper already states and cites, the same way the\n")
cat(sprintf("       normal-ordering floor is computed. That gives %.1f meV.\n", io*1000))
cat("\n   Take (b). It is internally consistent, it uses no unsourced number, and inverted\n")
cat("   ordering is context in this paper rather than a result: the claim is the NORMAL\n")
cat("   ordering floor, which section 1 reproduces exactly. The cost of (b) is that the\n")
cat("   quoted inverted sum is not the sharpest available; that is worth saying in the text\n")
cat("   and is cheaper than carrying a number with no provenance.\n")

cat("\n=== 5. the neighbouring numbers, checked while here\n\n")
mbb_lo <- 1.5; mbb_hi <- 3.7
cat(sprintf("   m_bb range quoted: %.1f-%.1f meV. Not recomputed here: it depends on the\n", mbb_lo, mbb_hi))
cat("   PMNS phases, which the paper does not state, so it is an input rather than a\n")
cat("   derived quantity and cannot be checked from what is written down.\n")
