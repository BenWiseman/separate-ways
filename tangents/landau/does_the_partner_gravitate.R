# ==========================================================================================
# "Do the sheets carry their own energy?" A reader meeting a two-sheeted cosmology will ask
# whether the expansion rate at the bang sees one sheet's content or two, and the paper does
# not answer it in the place the question arises. It matters quantitatively, because
# M_1 ~ g_*^(1/10) (gstar_sensitivity.R), so the answer is worth a number rather than a
# paragraph.
#
# The framework's position turns out to be a precise and slightly odd combination: the
# partner sheet is ESSENTIAL kinematically and ABSENT dynamically. Price both halves.
M0 <- 491.6; p <- 1/10

cat("  (1) KINEMATICALLY ESSENTIAL. The Landau-Zener closed form needs a sweep from\n")
cat("      eta = -infinity to +infinity. One sheet supplies half of one. 2.2 already\n")
cat("      computes what each gives:\n\n")
cat("        full sweep, two sheets : |beta|^2 = exp(-x^2), converges by x ~ 3\n")
cat("        half sweep, one sheet  : the bang-adiabatic gamma^2/16 p^4 tail, whose energy\n")
cat("                                 integral DIVERGES logarithmically\n\n")
cat("      So without the partner there is no finite production integral and hence no mass\n")
cat("      at all. The partner is not decoration; it is what makes the number exist.\n")

cat("\n  (2) DYNAMICALLY ABSENT. 1 states the second sheet is 'a copy related by an\n")
cat("      involution and not a place', and 4.2 shows a one-sided formulation supplied with\n")
cat("      the same state data reproduces the production calculation exactly. So the content\n")
cat("      is one sheet's worth and g_* counts it once.\n\n")
cat("      What it would cost if that were wrong, i.e. if both sheets sourced the same\n")
cat("      Friedmann equation and g_* doubled:\n\n")
cat("        g_*        M_1 (PeV)    line (PeV)   shift (PeV)   in widths of 2.0\n")
for (g in c(1,2)) {
  M <- M0*g^p
  cat(sprintf("   %8.1f %12.1f %12.1f %13.1f %16.1f\n", 106.75*g, M, M/2, M-M0, (M-M0)/2))
}
cat("\n  FLATLY: the ontological choice is worth 35.3 PeV, seventeen and a half times the propagated\n")
cat("  width, so it is not a philosophical preference. A framework in which both sheets\n")
cat("  gravitate predicts a ceiling of 526.9 PeV and a line at 263.4, and this one predicts\n")
cat("  491.6 and 245.8. The two are distinguishable by the same measurement that tests the\n")
cat("  line at all.\n")
cat("\n  WHAT THIS DOES NOT DO: it does not establish the copy reading over the two-place one.\n")
cat("  4.2 is explicit that the fold does not establish a second sheet as the only ontology,\n")
cat("  and this calculation does not change that. What it supplies is the price, so the\n")
cat("  choice stops being free.\n")
cat("\n  The combination is worth stating plainly because it is unusual: the partner is\n")
cat("  required for the crossing to have a finite answer and contributes nothing to the\n")
cat("  expansion rate. Kinematically indispensable, dynamically silent.\n")
