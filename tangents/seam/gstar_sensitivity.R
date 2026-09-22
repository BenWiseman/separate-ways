# ==========================================================================================
# An editorial pass asked the right question: is 491.6 PeV a number or a band? It proposed
# moving the radiation history by an effective g_* step and reporting the shift. The paper
# lists g_* = 106.75 as an input "treated as exact" and never gives its exponent, so a reader
# cannot propagate it. Recover the exponent and run the check.
#
# From calc/tangents/greybody/line_energy_error_budget.py, which takes it from BFT's own
# abundance function:
#     M_1 = ( rho_DM0 muhat^1.5 / (S_0 K) )^0.4,   K = (3 I / 2 pi^2) (15/g_star)^0.25
# so K ~ g_*^(-1/4) and M_1 ~ K^(-2/5):
p  <- 1/10
M0 <- 491.6; g0 <- 106.75
cat(sprintf("      M_1 proportional to g_*^(%.2f)\n\n", p))

# The direction that matters. Production is at T ~ 1e13 GeV (2.3), far above every Standard
# Model threshold, so g_* cannot be LOWER there: the QCD and e+e- values are for later epochs
# and are irrelevant to this calculation. What can happen is new content ADDING degrees of
# freedom. Tabulate that.
cat("      g_*      content at T ~ 1e13 GeV              M_1 (PeV)   shift (PeV)  in widths\n")
for (r in list(c(106.75,"Standard Model, adopted"),
               c(112.75,"+ one real scalar and a Weyl pair"),
               c(134.25,"+ a full extra generation"),
               c(228.75,"MSSM"),
               c(427.75,"MSSM plus a second hidden copy"))) {
  g <- as.numeric(r[1]); M <- M0*(g/g0)^p
  cat(sprintf("   %8.2f   %-34s %10.1f %12.1f %10.1f\n", g, r[2], M, M-M0, (M-M0)/2.0))
}
cat("\n  FLATLY, and it answers the question as asked: 491.6 PeV is a ceiling for ONE assumed\n")
cat("  particle content. A supersymmetric spectrum at production moves it to 530.5 PeV, a\n")
cat(sprintf("  shift of %.1f PeV or %.0f times the propagated width of 2.0 PeV.\n",
            M0*(228.75/g0)^p-M0, (M0*(228.75/g0)^p-M0)/2.0))
cat("  2.3 already says the width carries no allowance for the history being different. This\n")
cat("  is that allowance, with a number: the content dominates the budget by an order of\n")
cat("  magnitude over every measured input combined.\n")

cat("\n  WHY THE CEILING SURVIVES BEING WRONG ABOUT IT. The exponent is a tenth, so the\n")
cat("  dependence is very weak:\n\n")
for (f in c(2,4,10)) {
  cat(sprintf("     g_* larger by a factor of %2d  ->  mass larger by %.1f per cent (%.1f PeV)\n",
              f, 100*(f^p-1), M0*(f^p-1)))
}
cat("\n  Quadrupling the number of degrees of freedom at production moves the ceiling by 15\n")
cat("  per cent. That is the sense in which it is a band and not a number, and the band is\n")
cat("  narrow because every input enters at a fractional power.\n")
cat("\n  What this does NOT license: reading 491.6 as a prediction good to its quoted digits.\n")
cat("  The digits record how the arithmetic follows from the inputs, as 2.3 says. The\n")
cat("  physical statement is a ceiling near 490 PeV on Standard Model content at production,\n")
cat("  rising as the tenth power if there is more.\n")
