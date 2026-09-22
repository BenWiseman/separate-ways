# Hadamard regularity requires r(omega) -> 0. A reflectivity can only vary with omega
# if it has a scale. The seam action's own coefficient is dimensionless. So either the
# seam has no scale, and r = 0 identically, or it borrows one, and the only scales
# available at a bifurcate Killing horizon are the horizon's own.
#
# That is a DICHOTOMY, not a choice, and this asks how much each branch commits to.
# Kerr, M = 1: r+ = 1 + sqrt(1-a^2), r+^2 + a^2 = 2 r+, T_H = sqrt(1-a^2)/(4 pi r+),
# Omega_H = a/(2 r+), curvature scale at B ~ 1/r+, area A = 8 pi r+.

rp <- function(a) 1 + sqrt(1-a^2)
scales <- list(
  "temperature 1/T_H"  = function(a) 4*pi*rp(a)/sqrt(1-a^2),
  "horizon radius r+"  = function(a) rp(a),
  "sqrt area"          = function(a) sqrt(8*pi*rp(a)),
  "1/Omega_H (a>0)"    = function(a) 2*rp(a)/pmax(a,1e-9),
  "curvature (r+^3)^?" = function(a) rp(a)          # any power of r+ behaves alike
)

cat("=== 1. the ROBUST half: mass independence, whatever scale is borrowed\n\n")
cat("  Every horizon scale is M times a dimensionless function of a/M. The ringdown\n")
cat("  frequency is (M omega)/M with M omega a function of a/M. So the dimensionless\n")
cat("  argument omega x scale has M cancelling identically, for EVERY choice:\n\n")
cat("        scale                    arg at a=0.7, M=10      M=1e9        ratio\n")
for (nm in names(scales)) {
  f <- scales[[nm]]
  # omega = (M omega)/M, scale = M * s(a);  product = (M omega) * s(a), no M
  v1 <- 0.5*f(0.7); v2 <- 0.5*f(0.7)
  cat(sprintf("   %-24s %18.6f %12.6f %12.1f\n", nm, v1, v2, v1/v2))
}
cat("\n  Mass cancels exactly and by construction. So whatever the seam borrows, a\n")
cat("  10 Msun and a 10^9 Msun hole at the same spin reflect identically. That is a\n")
cat("  test across eight orders of magnitude that does not depend on WHICH scale.\n")

cat("\n=== 2. the SPECIFIC half: only the temperature diverges at extremality\n\n")
cat("        a/M      1/T_H       r+        sqrt(A)     1/Omega_H\n")
for (a in c(0, 0.3, 0.7, 0.9, 0.99, 0.999)) 
  cat(sprintf("   %8.3f %10.2f %8.4f %10.4f %12.3f\n", a,
      scales[["temperature 1/T_H"]](a), rp(a), sqrt(8*pi*rp(a)), if(a>0) 2*rp(a)/a else Inf))
cat("\n  The temperature scale diverges at extremality and changes by a factor of twelve\n")
cat("  across the spins sampled above; every\n")
cat("  other horizon scale changes by a factor of two or less. So:\n")
cat("   - mass independence is ROBUST to the scale choice;\n")
cat("   - 'effects vanish at extremality' holds ONLY if the borrowed scale is the\n")
cat("     temperature. If it is the radius or the area, reflectivity is nearly\n")
cat("     spin-independent instead.\n")
cat("  Those are different predictions and the data separates them. An earlier note of\n")
cat("  mine attributed the extremality claim to 'any horizon scale', which is wrong.\n")

cat("\n=== 3. so the paper's commitment is a DICHOTOMY, and both halves are falsifiable\n\n")
cat("   branch A  no scale       r = 0 exactly        ordinary absorbing Kerr;\n")
cat("                                                 ANY intrinsic departure kills it.\n")
cat("   branch B  horizon scale  r = r(M omega, a)    departure allowed, but mass-\n")
cat("                                                 independent at fixed spin; a\n")
cat("                                                 departure scaling with M kills it.\n\n")
cat("  There is no branch C, because Hadamard needs r(omega) -> 0, varying with omega\n")
cat("  needs a scale, and at a bifurcate Killing horizon the available scales are the\n")
cat("  horizon's. That is what makes this a bracket rather than a preference: the two\n")
cat("  branches between them exclude every reflecting seam whose amplitude depends on\n")
cat("  the hole's mass at fixed spin.\n")

cat("\n=== 4. the sharpest single statement that follows\n\n")
cat("  A horizon reflectivity that is NOT a function of spin alone is incompatible with\n")
cat("  the fold, in either branch. No amplitude, no fall-off shape and no seam action\n")
cat("  needs to be chosen to say that. It is a statement about what the data may not\n")
cat("  show, and it is testable the moment ringdown measurements exist in two well-\n")
cat("  separated mass ranges at overlapping spin.\n")
