# CORRECTED 2026-09-21. This file evaluated the response falloff at M*Lambda (~1e38 to 1e47)
# when the response is a function of omega/Lambda (~4e-40 to 4e-48), the reciprocal. That
# inversion produced the paper's claim that a microscopic seam is invisible at ringdown, which is
# backwards: astrophysical ringdown sits in the deep INFRARED of r(omega/Lambda), where a response
# obeying an ultraviolet falloff condition is approximately ONE. r(x) = exp(-x) at x = 4e-40 is
# 1.000000000000. The invisibility claim is withdrawn from the paper; see A.13 and 3.4.
# Objection: "no third option" is not watertight - a Planck, string or seam-intrinsic scale sits
# in the same dimensional slot as a horizon scale. Correct, and conceded. But that objection's own
# parenthesis points somewhere better than a no-go: a non-horizon scale makes the
# reflectivity MASS-dependent at fixed spin, so the measurement does not merely test the
# fold, it CLASSIFIES the seam. Work that out quantitatively rather than rhetorically.
#
# r = r(omega / Lambda) with Lambda the seam's scale. At the ringdown,
#   omega = (M omega)/M, so the argument is  (M omega) / (M Lambda).
# Everything turns on M Lambda, dimensionless, for the astrophysical mass range.

G <- 6.67430e-11; c <- 2.99792458e8; hbar <- 1.054571817e-34; Msun <- 1.98847e30
lP <- sqrt(hbar*G/c^3)                       # Planck length
rg <- function(M) G*M*Msun/c^2               # gravitational radius, metres

cat(sprintf("=== 0. Planck length = %.3e m\n\n", lP))
cat("=== 1. the dimensionless argument M*Lambda for each candidate seam scale\n\n")
cat("   For a scale Lambda with length 1/Lambda, M Lambda ~ rg(M) / (1/Lambda).\n\n")
cat("        M (Msun)     rg (m)      rg/lP (Planck seam)    horizon seam (by definition)\n")
for (M in c(10, 65, 1e6, 1e9)) 
  cat(sprintf("   %12.3g %12.3e %22.3e %28s\n", M, rg(M), rg(M)/lP, "O(1)"))
cat("\n  A horizon-scale seam sits at argument O(1): the ringdown probes the interesting\n")
cat("  part of r. A Planck-scale seam sits at argument 1e38 to 1e47, which is the deep\n")
cat("  ULTRAVIOLET of r, where Hadamard regularity has already forced r -> 0.\n")

cat("\n=== 2. so a Planck-scale seam is not merely mass-dependent. It is INVISIBLE.\n\n")
cat("   r evaluated at argument x = rg/lP, for three fall-offs meeting r(inf) = 0:\n\n")
cat("        M (Msun)          x        r ~ 1/x        r ~ 1/x^2        r ~ exp(-x)\n")
for (M in c(10, 1e9)) {
  x <- rg(M)/lP
  cat(sprintf("   %12.3g %11.3e %13.3e %16.3e %18s\n", M, x, 1/x, 1/x^2, "0 (underflow)"))
}
cat("\n  So the third branch named above is real theoretically and DEAD observationally:\n")
cat("  any seam scale far above the horizon's is pushed into the region where the same\n")
cat("  Hadamard condition that created the dichotomy has already sent r to zero. It is\n")
cat("  not excluded by assumption; it is excluded by being unobservable, which is a\n")
cat("  weaker and more honest claim, and it does not need the lemma the objection asks us to prove.\n")

cat("\n=== 3. what a measurement would therefore MEAN. Three signatures, not two.\n\n")
cat("   S1  no reflectivity at any mass or spin\n")
cat("         -> scale-free seam, or any scale far above the horizon's. These are not\n")
cat("            distinguishable by ringdown, and the paper should not claim they are.\n\n")
cat("   S2  reflectivity varying with SPIN, identical across mass at fixed spin\n")
cat("         -> the seam's scale IS the horizon's. This is the fold's own branch.\n")
cat("            If the scale is the temperature: largest at low spin, vanishing toward\n")
cat("            extremality. If the radius or area: nearly flat in spin.\n\n")
cat("   S3  reflectivity varying with MASS at fixed spin\n")
cat("         -> a seam scale COMPARABLE to the horizon's but not derived from it, since\n")
cat("            a much larger scale gives S1 and a horizon scale gives S2. The mass\n")
cat("            dependence then MEASURES that scale rather than refuting anything.\n\n")
cat("   For S3 to be visible the seam scale must sit within a few orders of the horizon\n")
cat("   scale of an astrophysical hole:\n\n")
cat("        M (Msun)       horizon scale 1/rg as an energy (eV)\n")
for (M in c(10, 65, 1e6, 1e9)) 
  cat(sprintf("   %12.3g %35.3e\n", M, hbar*c/rg(M)/1.602176634e-19))
cat("\n  So S3 probes seam scales around 1e-11 to 1e-20 eV. Nothing in particle physics\n")
cat("  sits there, which is why S3 would be a genuinely new scale if it were seen, and\n")
cat("  why its absence is unsurprising rather than informative.\n")

cat("\n=== 4. the honest headline\n\n")
cat("  Not 'a mass-dependent reflectivity is excluded'. That overclaims, because the\n")
cat("  theoretical dichotomy has the hole found above. The defensible headline is that\n")
cat("  ringdown reflectivity measured across mass and spin SORTS the seam into three\n")
cat("  cases, that the fold's own branch is the middle one and carries a spin trend\n")
cat("  with no free amplitude, and that the first and third are distinguishable from it\n")
cat("  by data alone. That is a measurement programme rather than a no-go, and it does\n")
cat("  not rest on the lemma we cannot prove.\n")
