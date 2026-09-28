# Section 4: the generalized Boltzmann reflectivity and its spin dependence.
#
# NOT AN INDEPENDENT CONFIRMATION. The co-rotation error was found first by
# separate_ways/tangents/ringdown/kerr_reflectivity_corotation.R, which computes the same
# |R|(a) from the same fit. This file exists because that one was written against the
# COSMOLOGY paper's claim and prints a verdict about it; what is needed here is a check of
# the four sentences Section 4 actually writes, with the assertions to match. Same arithmetic,
# different claims checked. Do not cite both as if they agreed independently.
#
# Amaro Seoane [13] gives a horizon reflectivity for the black-mirror model equal to the
# generalized Boltzmann factor. Section 4 says one of the bounds placed on it was arithmetic
# and was wrong, and quotes five values to say so. A paper that corrects someone else's
# arithmetic has to be able to show its own, which is what this file is for.
#
# The factor carries the co-rotation term:
#     |R| = exp[ -|omega - m Omega_H| / (2 T_H) ]
# Dropping m Omega_H reverses the spin dependence, because omega and m Omega_H both grow
# with spin and the difference between them shrinks.
#
# Kerr, M = 1, a in units of M:
#     r_+    = 1 + sqrt(1 - a^2)
#     Om_H   = a / (r_+^2 + a^2) = a / (2 r_+)
#     kappa  = (r_+ - r_-) / (2 (r_+^2 + a^2)) = sqrt(1-a^2) / (2 r_+)
#     T_H    = kappa / (2 pi)
# Fundamental l = m = 2 mode from the Berti-Cardoso-Will fit already used in Section 4:
#     M omega = 1.5251 - 1.1568 (1-a)^0.1292

rp    <- function(a) 1 + sqrt(1 - a^2)
OmH   <- function(a) a / (2 * rp(a))
TH    <- function(a) sqrt(1 - a^2) / (4 * pi * rp(a))
Mw    <- function(a) 1.5251 - 1.1568 * (1 - a)^0.1292
Rmag  <- function(a, m = 2) exp(-abs(Mw(a) - m * OmH(a)) / (2 * TH(a)))
Rnaive<- function(a)        exp(-abs(Mw(a))             / (2 * TH(a)))

cat("=== generalized Boltzmann reflectivity, l = m = 2 ===\n\n")
cat("     a      M omega     Om_H        T_H       |R| full   |R| no mOm\n")
for (a in c(0, 0.7, 0.9, 0.95, 0.99)) {
  cat(sprintf("  %5.2f  %9.6f  %9.6f  %9.6f  %10.6f  %10.3e\n",
              a, Mw(a), OmH(a), TH(a), Rmag(a), Rnaive(a)))
}

cat("\n  the five values quoted in Section 4:\n")
q <- c("0"="9.8e-3", "0.7"="0.148", "0.9"="0.444", "0.95"="0.662", "0.99"="0.372")
for (a in c(0, 0.7, 0.9, 0.95, 0.99))
  cat(sprintf("    a = %-5s |R| = %.4f   (paper says %s)\n", a, Rmag(a), q[[as.character(a)]]))

cat(sprintf("\n  at a = 0 the two expressions agree: %.6f vs %.6f\n", Rmag(0), Rnaive(0)))
stopifnot(abs(Rmag(0) - Rnaive(0)) < 1e-12)

cat(sprintf("  energy returned at a = 0.9: |R|^2 = %.4f, about a fifth\n", Rmag(0.9)^2))

# Not monotone: it turns over between 0.95 and 0.99 on this fit.
gr <- seq(0, 0.999, by = 0.0005); pk <- gr[which.max(Rmag(gr))]
cat(sprintf("  maximum of |R| on [0,0.999] at a = %.4f, |R| = %.4f, so the rise is not\n",
            pk, Rmag(pk)))
cat("  monotone to extremality and Section 4 states it only for the interval computed.\n")
stopifnot(pk < 0.99, Rmag(0.95) > Rmag(0.9), Rmag(0.99) < Rmag(0.95))

# The damping-time model, and the temperature such a seam would need.
cat("\n=== the most generous damping model ===\n")
dtau <- 0.14
Rneed <- sqrt(1 - 1 / (1 + dtau))
cat(sprintf("  tau -> tau/(1-|R|^2) reproduces delta-tau-hat_220 = %.2f at |R| = %.4f\n", dtau, Rneed))
for (a in c(0.9, 0.95)) {
  Tneed <- abs(Mw(a) - 2 * OmH(a)) / (2 * log(1 / Rneed))
  cat(sprintf("  a = %.2f: that needs T = %.4f T_H, cooler than Hawking\n", a, Tneed / TH(a)))
}

# A 10 Msun and a 10^9 Msun hole at the same spin reflect identically, because both
# M omega and M T_H depend on a/M alone. Stated in the appendix as a null test.
cat(sprintf("\n  mass independence: omega/T_H at a=0.7 is %.6f for any mass.\n",
            Mw(0.7) / TH(0.7)))
