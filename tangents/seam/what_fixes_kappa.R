# ==========================================================================================
# The seam coefficient is the paper's largest admitted gap: 4.1 says the algebra has not paid
# for it, and the media pass named it one of three open holes. Unitarity selects the SIGN of
# the branch but not kappa, since r = (1-k^2)/(1+k^2), t = 2k/(1+k^2) gives r^2+t^2 = 1 at
# every kappa. So ask what else in the construction could fix it, and test each candidate
# instead of arguing.
#
# The family is a half-angle: kappa = tan(theta/2) gives r = cos(theta), t = sin(theta).
kap <- c(0, 0.25, 0.5, tan(pi/8), 1, 2, 4, 1e6)
r <- (1-kap^2)/(1+kap^2); t <- 2*kap/(1+kap^2); th <- 2*atan(kap)
cat("      kappa        r        t     r^2+t^2    theta/pi   2*atan(kappa) check\n")
for (i in seq_along(kap))
  cat(sprintf("   %9.4f %8.4f %8.4f %10.6f %10.4f   cos=%.4f sin=%.4f\n",
      kap[i], r[i], t[i], r[i]^2+t[i]^2, th[i]/pi, cos(th[i]), sin(th[i])))
cat("\n  So the seam is a rotation by theta = 2 arctan(kappa). kappa = 1 is theta = pi/2,\n")
cat("  which is r = 0, t = 1: the transparent point the paper ADOPTS rather than derives.\n")

cat("\n  ==========================================================================\n")
cat("  CANDIDATE: does involutivity of the seam map fix kappa? The fold is an involution.\n")
cat("  ==========================================================================\n")
cat("  Two ways a 2x2 orthogonal map can act on the pair of legs.\n\n")
for (nm in c("rotation","reflection")) {
  cat(sprintf("  --- seam acts as a %s\n", nm))
  cat("      kappa    theta/pi    ||S^2 - 1||    involutive?\n")
  for (k in c(0, 0.25, 0.5, 1, 2, 1e6)) {
    thk <- 2*atan(k); c_ <- cos(thk); s_ <- sin(thk)
    S <- if (nm=="rotation") matrix(c(c_,-s_,s_,c_),2,2) else matrix(c(c_,s_,s_,-c_),2,2)
    d <- max(abs(S%*%S - diag(2)))
    cat(sprintf("   %9.4f %10.4f %14.2e    %s\n", k, thk/pi, d, ifelse(d<1e-12,"yes","NO")))
  }
  cat("\n")
}
cat("  FLATLY. If the seam acts as a ROTATION, S^2 = R(2 theta) = 1 only at theta = 0 or pi,\n")
cat("  i.e. kappa = 0 or infinity, the two perfect mirrors. Transparency is then EXCLUDED by\n")
cat("  involutivity. If it acts as a REFLECTION, S^2 = 1 identically and kappa stays free.\n")
cat("\n  So involutivity does not fix kappa, but it is not empty either: it forces a choice\n")
cat("  between two realisations, and the transparent matching the paper adopts commits it to\n")
cat("  the reflection one. That is a stated commitment with a consequence, which is better\n")
cat("  than an unexplained choice, and it is falsifiable inside the formalism: anything that\n")
cat("  showed the seam must act as a rotation would force a perfect mirror and kill\n")
cat("  transparent matching outright.\n")
cat("\n  WHAT STILL DOES NOT FIX KAPPA: unitarity (holds for all), involutivity (holds for all,\n")
cat("  in the reflection realisation), and the half-angle structure itself, which is generic\n")
cat("  to any SU(2) and is NOT evidence that this angle is the bang's angle of 3.1. Do not\n")
cat("  claim that identification; both are two-level problems and the half-angle is what\n")
cat("  two-level problems look like.\n")
