# If the seam carries the HORIZON'S OWN scale instead of no scale, the UV condition
# A.13 imposes is met without forcing transparency, and the reflectivity becomes a
# function of omega / T_H. Two consequences follow, and the first needs no quasinormal
# frequencies at all.
#
# Kerr, units M = 1:  r+ = 1 + sqrt(1-a^2),  r+^2 + a^2 = 2 r+,
#   surface gravity  kappa_H = sqrt(1-a^2)/(2 r+),   T_H = kappa_H/(2 pi).
# So   M T_H = sqrt(1-a^2) / (4 pi r+)   and
#      omega/T_H = (M omega) * 4 pi (1 + sqrt(1-a^2)) / sqrt(1-a^2).

rp  <- function(a) 1 + sqrt(1-a^2)
MTH <- function(a) sqrt(1-a^2)/(4*pi*rp(a))
S   <- function(a) 4*pi*(1+sqrt(1-a^2))/sqrt(1-a^2)     # omega/T_H per unit M*omega

cat("=== 1. the consequence that needs NO quasinormal data\n\n")
cat("  M*omega and M*T_H are both dimensionless functions of a/M alone, so their ratio\n")
cat("  omega/T_H is a function of SPIN ONLY. A seam whose reflectivity is r(omega/T_H)\n")
cat("  therefore reflects identically for a 10 Msun and a 10^9 Msun hole at the same\n")
cat("  spin. Mass cancels exactly:\n\n")
for (M in c(10, 65, 1e6, 1e9)) cat(sprintf("   M = %9.3g Msun, a = 0.7 : (M omega)/(M T_H) = %.6f x (M omega)\n", M, S(0.7)))
cat("\n  Identical to every digit, by construction, because M never appears in S(a).\n")
cat("  That is a test across eight orders of magnitude in mass with no free parameter,\n")
cat("  and it is available from LIGO-band and LISA-band ringdowns in matching spin bins.\n")

cat("\n=== 2. the spin dependence, exactly\n\n")
cat("        a/M       M T_H      omega/T_H per unit M omega    relative to a=0\n")
for (a in c(0, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 0.998)) 
  cat(sprintf("   %8.3f %11.6f %26.3f %17.2f\n", a, MTH(a), S(a), S(a)/S(0)))
cat("\n  S(0) = 8 pi = 25.13. The band the ringdown probes moves to HIGHER omega/T_H as\n")
cat("  the hole spins up, and diverges at extremality because T_H -> 0 there.\n")

cat("\n=== 3. so how much of that is the quasinormal frequency, and how much is T_H?\n")
cat("   M omega for the l = m = 2 fundamental is NOT computed here. Bound it crudely\n")
cat("   and see whether the conclusion survives the worst case:\n\n")
cat("        a/M      S(a)     omega/T_H if M omega = 0.3   ... if M omega = 1.0\n")
for (a in c(0, 0.5, 0.9, 0.99, 0.998))
  cat(sprintf("   %8.3f %9.2f %26.1f %20.1f\n", a, S(a), 0.3*S(a), 1.0*S(a)))
cat("\n  Across the whole spin range M omega changes by a factor of a few; S(a) changes\n")
cat("  by a factor of %s. The divergence is carried by T_H and not by the mode, so the\n")
cat("  conclusion does not depend on quasinormal values I have not computed.\n")
cat(sprintf("  (S(0.998)/S(0) = %.1f)\n", S(0.998)/S(0)))

cat("\n=== 4. THE PREDICTION, and it runs opposite to the echo literature\n\n")
cat("  Any seam meeting A.13's condition has r(x) -> 0 as x = omega/T_H -> infinity.\n")
cat("  Spinning the hole up pushes the ringdown into larger x. Therefore horizon\n")
cat("  reflectivity effects are LARGEST at LOW spin and VANISH as extremality is\n")
cat("  approached, whatever the function r is:\n\n")
for (rf in list(c("r = exp(-x/10)", function(x) exp(-x/10)), c("r = 1/(1+x/10)", function(x) 1/(1+x/10)),
                c("r = (1+ (x/10)^2)^-1", function(x) 1/(1+(x/10)^2)))) {
  f <- rf[[2]]; cat(sprintf("   %-24s", rf[[1]]))
  for (a in c(0, 0.5, 0.9, 0.99)) cat(sprintf(" a=%.2f:%7.4f", a, f(0.5*S(a))))
  cat("\n")
}
cat("\n  Three unrelated fall-offs, same ordering, monotone decreasing in spin. Echo\n")
cat("  models built on a reflecting wall at fixed proper distance predict the OPPOSITE,\n")
cat("  because the light-crossing delay grows as the horizon is approached and the\n")
cat("  effect becomes easier to see at high spin. The two are distinguishable by the\n")
cat("  SIGN of the trend with spin alone, with no amplitude needed from either side.\n")

cat("\n=== 5. what this does not settle\n\n")
cat("  Whether the seam has the horizon's scale at all. A.13's corner term does not,\n")
cat("  and within it kappa = 1 is forced. This is the one way out of that, and it is a\n")
cat("  different action, so it owes its own matching law. What is shown here is that\n")
cat("  the way out is not free: it commits to a falsifiable spin trend, and to mass\n")
cat("  independence, before any amplitude is chosen.\n")
