# F-AO claimed A.13's Hadamard argument cannot see a NON-LINEAR seam, because Hadamard is a
# short-distance condition and the two-point function is "a small-amplitude object". Ben
# asked me to challenge assumptions. That one is wrong, and it is the load-bearing step.
#
# VACUUM FLUCTUATIONS DIVERGE AT SHORT DISTANCE. In 4D, <phi(x)phi(y)> ~ 1/(4 pi^2 s^2), so
# the field amplitude a probe of separation s samples goes as 1/(2 pi s) and grows without
# bound as s -> 0. Short distance is not small amplitude, it is LARGE amplitude. So a seam
# that turns on above some amplitude A0 turns ON for short-distance probes, which is exactly
# where Hadamard regularity looks.

Arms <- function(s) 1/(2*pi*s)                  # vacuum rms amplitude at separation s
rA   <- function(A, A0, p=2) (A/A0)^p/(1+(A/A0)^p)

cat("=== 1. the amplitude a short-distance probe actually samples\n\n")
cat("        separation s      vacuum rms amplitude\n")
for (s in c(1, 1e-2, 1e-4, 1e-8, 1e-12)) cat(sprintf("   %14.0e %24.3e\n", s, Arms(s)))
cat("\n  Divergent. My premise had this backwards.\n")

cat("\n=== 2. so where does a non-linear seam turn on, for any A0?\n\n")
cat("   The crossing separation is s* = 1/(2 pi A0): below it the seam is reflecting.\n\n")
cat("        A0          s* = 1/(2 pi A0)      r at s = s*/10      r at s = s*/100\n")
for (A0 in c(1, 1e3, 1e10, 1e30)) {
  ss <- 1/(2*pi*A0)
  cat(sprintf("   %10.0e %20.3e %20.4f %19.4f\n", A0, ss, rA(Arms(ss/10),A0), rA(Arms(ss/100),A0)))
}
cat("\n  For EVERY finite A0 the reflectivity tends to one as the separation shrinks. A\n")
cat("  non-linear seam is not invisible at short distance; it is MAXIMALLY reflecting\n")
cat("  there, which is worse than the constant case rather than better.\n")

cat("\n=== 3. therefore Hadamard closes the amplitude axis after all\n\n")
cat("  A.13's argument needs the image term to be present as the arguments approach the\n")
cat("  singular locus. With r -> 1 there, the image term arrives at FULL strength, so the\n")
cat("  divergence between spacelike-separated mirror points is not merely present but\n")
cat("  saturated. The conclusion is stronger for the non-linear seam than for the constant\n")
cat("  one, not weaker.\n")
cat("  So F-AO's 'axis the central argument does not cover' does not exist. I opened it by\n")
cat("  conflating SMALL AMPLITUDE with SHORT DISTANCE, and they are opposites in the\n")
cat("  vacuum.\n")

cat("\n=== 4. and this is the scale-dependence, stated properly\n\n")
cat("  The join is scale-dependent, which is what Ben has been saying, but the dependence\n")
cat("  runs the other way from the picture we were both carrying. A seam that reflects at\n")
cat("  large amplitude is transparent at long distance and reflecting at short - not\n")
cat("  transparent in the weak-field exterior and active near a singularity, but\n")
cat("  transparent in the infrared and active in the ultraviolet.\n")
cat("  That is an RG statement and it has the expected sign: a non-linear boundary coupling\n")
cat("  carries a dimensionful scale, so it is IRRELEVANT in the infrared and grows toward\n")
cat("  the ultraviolet. Hadamard regularity is a UV condition, which is precisely why it\n")
cat("  catches it.\n")

cat("\n=== 5. what survives of the pinch picture\n\n")
cat("  The CAUSAL statement is untouched: A.14's F is still the only region where both\n")
cat("  sheets can influence a common event, and that had nothing to do with the seam's\n")
cat("  amplitude response. What does not survive is a seam that is transparent where we\n")
cat("  look and reflecting where we cannot, because the amplitude that would switch it on\n")
cat("  is reached in the vacuum at short distance everywhere, not only near a singularity.\n")
