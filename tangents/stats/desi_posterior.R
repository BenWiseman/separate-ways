# Same arithmetic as the KM3NeT check, applied to the paper's most load-bearing observational
# claim. The floor is EXACT at 58.8 meV; the bound is a posterior upper limit, not a hard edge.
# Elbers et al. (DESI DR2): Sigma m_nu < 64.2 meV under LCDM, < 163 meV under w0waCDM, and a
# Feldman-Cousins limit of 53 meV. The paper quotes the survival condition (an 8.4 per cent
# tightening reaches the floor; it read 8.1 before the bound was corrected from 64.0 to 64.2
# on 2026-09-21) but not how much of the CURRENT posterior already sits below us.

floor <- 58.8; lim95 <- 64.2; fc <- 53.0
cat(sprintf("   floor %.1f meV (exact, normal ordering, m1 = 0)\n", floor))
cat(sprintf("   DESI DR2 LCDM 95%% upper limit %.1f meV; Feldman-Cousins %.1f meV\n\n", lim95, fc))

cat("   P(Sigma m_nu > floor) under three posterior shapes all matched to the 95% limit:\n\n")
cat("      shape                              parameter        P(> 58.8 meV)\n")
lam <- -log(0.05)/lim95
p1 <- exp(-lam*floor)
cat(sprintf("      exponential from zero           lambda = %.5f /meV %12.4f\n", lam, p1))
sg <- lim95/qnorm(0.975)
p2 <- 2*(1-pnorm(floor/sg))
cat(sprintf("      half-Gaussian at zero           sigma  = %.2f meV   %12.4f\n", sg, p2))
p3 <- (lim95-floor)/lim95
cat(sprintf("      flat, truncated at the limit    width  = %.1f meV   %12.4f\n", lim95, p3))
cat(sprintf("\n      range across shapes: %.1f to %.1f per cent of the posterior lies above the floor.\n",
    100*min(p1,p2,p3), 100*max(p1,p2,p3)))

cat("\n=== what that means, stated both ways\n\n")
cat(sprintf("   Put positively: the floor is ALLOWED. It sits below the 95 per cent limit by\n"))
cat(sprintf("   %.1f meV, and the data have not excluded it.\n", lim95-floor))
cat(sprintf("   Put honestly: only %.0f to %.0f per cent of the posterior is compatible with a sum\n",
    100*min(p1,p2,p3), 100*max(p1,p2,p3)))
cat("   at or above the floor, so the prediction sits in the upper tail of what DESI allows\n")
cat("   rather than comfortably inside it. Both statements are true and the second is the one\n")
cat("   a referee will compute.\n")

cat("\n=== and the frequentist limit is already past us\n\n")
cat(sprintf("   The Feldman-Cousins limit is %.0f meV, which is %.1f meV BELOW the floor. On that\n",
    fc, floor-fc))
cat("   construction the prediction is already excluded. The paper does not read that as an\n")
cat("   exclusion, and the reason has to be stated precisely rather than waved at: a\n")
cat("   frequentist limit below the oscillation minimum is a statement about EVERY normally\n")
cat("   ordered spectrum, not about this model, since 58.8 meV is the kinematic floor for\n")
cat("   normal ordering whatever the cosmology. A limit that excludes the minimum allowed by\n")
cat("   oscillation data is currently read as a feature of the data.\n")
cat(sprintf("   That defence is available to us only because the floor IS the oscillation\n"))
cat("   minimum. If the model predicted a sum above that minimum, the same limit would bite\n")
cat("   on the model alone and the defence would evaporate.\n")

cat("\n=== the survival condition, restated with the posterior in view\n\n")
need <- floor/lim95
cat(sprintf("   A tightening of the 95 per cent limit by %.1f per cent, from %.1f to %.1f meV,\n",
    100*(1-need), lim95, floor))
cat("   puts the floor at the limit itself. Combined with the tail fraction above, the corner\n")
cat("   is not merely close: the model already occupies the part of the posterior the data\n")
cat("   least prefer, and the next release moves the limit in one direction only.\n")
