# The Gregory-Laflamme scale, and the one conversion Section 5.2 actually uses.
#
# Section 5.2 needs a crossover between the regime where a hole wrapped on a compact
# dimension is a uniform string, where the angular budget is the four-dimensional 2*pi,
# and the regime where it is a localised higher-dimensional hole, where the budget is
# 2*pi/(n+1) and falls short of the bill. What it needs from Gregory-Laflamme is that
# such a crossover exists and sits at r_h/R of order one. It does NOT need the threshold
# to be a phase boundary, and an earlier version of the paragraph wrongly took it for one.
#
# The threshold wavenumber for the neutral five-dimensional uniform string is quoted from
# the literature: k_GL r_h = 0.876 (Gregory-Laflamme 1993; Gubser hep-th/0110193 sec 3.2
# and footnote 6; independently Frolov and Shoom 0903.2893 sec IV). This script does not
# re-solve that eigenvalue problem and says so. What it does check is the conversion the
# manuscript quotes, and that the conversion is not self-fulfilling.

k_gl  <- 0.876                       # quoted, not derived here
L_crit <- 2 * pi / k_gl              # critical circumference in units of r_h
rh_over_R <- 2 * pi / L_crit         # with L = 2 pi R, this is r_h/R at threshold

cat("=== Gregory-Laflamme scale, neutral 5D uniform string ===\n\n")
cat(sprintf("  k_GL r_h              = %.3f   (quoted from the literature, not solved here)\n", k_gl))
cat(sprintf("  critical circumference = %.4f r_h\n", L_crit))
cat(sprintf("  equivalently r_h/R     = %.4f at threshold\n", rh_over_R))
cat(sprintf("  round trip k_GL -> L -> k_GL leaves %.2e\n\n", abs(2 * pi / L_crit - k_gl)))

stopifnot(abs(L_crit - 7.17) < 0.01, abs(rh_over_R - k_gl) < 1e-12)

cat("  planted failures, so the conversion cannot pass by matching nothing:\n")
bad_L <- 2 * pi / 0.5
cat(sprintf("    a threshold of 0.5 would give %.4f r_h, which misses 7.17 by %.2f\n",
            bad_L, abs(bad_L - 7.17)))
bad_conv <- pi / k_gl
cat(sprintf("    dropping the factor of two gives %.4f r_h, out by %.2f\n",
            bad_conv, abs(bad_conv - 7.17)))
cat(sprintf("    and r_h/R computed as 1/L rather than 2 pi/L gives %.4f, not %.3f\n\n",
            1 / L_crit, k_gl))

cat("  What Section 5.2 takes from this: a crossover exists and sits at r_h/R of order\n")
cat("  one. Where the handover finishes is a question about the nonlinear endpoint, which\n")
cat("  Gregory and Laflamme raise as a conjecture and the later evolutions answer with a\n")
cat("  cascade of necks and bulges rather than a clean threshold. Nothing in the angular\n")
cat("  budget depends on that answer.\n")
