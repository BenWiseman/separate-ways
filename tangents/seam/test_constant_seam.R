# test_constant_seam.R -- does the seam's image weight really multiply as A_I = -t A_0?
#
# The candidate answer is (B), transmission, from a constant corner-term seam
# with frequency-independent r(kappa) = (1-k^2)/(1+k^2), t(kappa) = 2k/(1+k^2).
# An earlier delta-barrier test did NOT reproduce A_I = -t A_0, but a delta
# barrier has t = t(k), so it was not the seam being derived here. This tests that one.
#
# 1D massive scalar, ports at x = -a and x = +a, scattering states
#   psi_L(x) = e^{ikx} + r e^{-ikx}   (x<0),   t e^{ikx}          (x>0)
#   psi_R(x) = t e^{-ikx}             (x<0),   e^{-ikx} + r e^{ikx} (x>0)
# with r, t real, constant, r^2 + t^2 = 1.
#
# Analytically, summing the two channels at the two ports:
#   local  A_0 ~ integral (2 + 2 r cos 2ka) / (2 omega)
#   cross  A_I ~ integral (2 t (cos 2ka + r)) / (2 omega)
# so with C = <cos 2ka> / <1>, both averages weighted by 1/(2 omega),
#   A_I / A_0 = t (C + r) / (1 + r C),
# which equals t exactly when C = 1 and not otherwise. The claim is therefore
# right in the long-wavelength limit and carries a computable correction away
# from it. That is the thing to check.

kmax <- 400; n <- 400000; m <- 1
k <- seq(1e-8, kmax, length.out = n); dk <- k[2]-k[1]; w <- sqrt(k^2 + m^2)

ratio_num <- function(a, r, t) {                  # A_I / A_0 by direct mode sum
  loc   <- sum((2 + 2*r*cos(2*k*a)) / (2*w)) * dk
  cross <- sum((2*t*(cos(2*k*a) + r)) / (2*w)) * dk
  cross / loc
}
Cfun <- function(a) sum(cos(2*k*a)/(2*w))*dk / (sum(1/(2*w))*dk)

cat("=== 1. the closed form A_I/A_0 = t(C+r)/(1+rC), against the mode sum\n\n")
cat("    a     kappa       r        t        C      mode sum   closed form     diff\n")
for (a in c(0.02, 0.1, 0.3, 1.0)) {
  C <- Cfun(a)
  for (kap in c(0.2, 0.5, 1.0, 2.0)) {
    r <- (1-kap^2)/(1+kap^2); t <- 2*kap/(1+kap^2)
    num <- ratio_num(a, r, t); cf <- t*(C+r)/(1+r*C)
    cat(sprintf("%6.2f  %6.2f  %7.4f  %7.4f  %7.4f  %10.6f  %12.6f  %9.2e\n",
                a, kap, r, t, C, num, cf, abs(num-cf)))
  }
}

cat("\n=== 2. does it collapse to A_I/A_0 = t? only as C -> 1\n\n")
cat("     a         C      max over kappa of |A_I/A_0 - t|\n")
for (a in c(1.0, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001)) {
  C <- Cfun(a)
  d <- max(sapply(c(0.1,0.2,0.5,1.0,2.0,5.0), function(kap) {
    r <- (1-kap^2)/(1+kap^2); t <- 2*kap/(1+kap^2)
    abs(ratio_num(a,r,t) - t) }))
  cat(sprintf("%8.3f  %8.5f  %26.6f\n", a, C, d))
}

cat("\n=== 3. the de Sitter case A_I = -A_0 is kappa = 1, and only kappa = 1\n\n")
cat("   kappa       t        r     A_I/A_0 at C=1      even port A_0+A_I\n")
for (kap in c(0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 10.0)) {
  r <- (1-kap^2)/(1+kap^2); t <- 2*kap/(1+kap^2)
  cat(sprintf("%8.2f %8.4f %8.4f %15.6f %20.6f\n", kap, t, r, t, 1-t))
}
cat("\n  even/odd = (1-t)/(1+t) = ((1-kappa)/(1+kappa))^2, zero only at kappa = 1,\n")
cat("  which is the transparent seam r = 0 and reproduces A_I = -A_0.\n")
for (kap in c(0.1,0.5,1,2,10)) {
  t <- 2*kap/(1+kap^2)
  cat(sprintf("    kappa=%5.2f  (1-t)/(1+t) = %.8f   ((1-k)/(1+k))^2 = %.8f\n",
      kap, (1-t)/(1+t), ((1-kap)/(1+kap))^2))
}

cat("\n=== 4. the Mobius form respects positivity of the two-port noise matrix\n\n")
# A_I/A_0 = t (C+r)/(1+rC) is t times a Mobius map of the unit disc to itself,
# so |A_I| <= |A_0| automatically and both eigenvalues A_0 +/- A_I stay >= 0.
worst <- 1
for (a in c(0.001,0.01,0.1,0.3,1,3)) {
  C <- Cfun(a)
  for (kap in c(0.05,0.1,0.25,0.5,1,2,4,10,20)) {
    r <- (1-kap^2)/(1+kap^2); t <- 2*kap/(1+kap^2)
    q <- t*(C+r)/(1+r*C)
    worst <- min(worst, 1-abs(q))
  }
}
cat(sprintf("  min over 54 (a, kappa) of 1 - |A_I/A_0| = %.6f   (>= 0 required)\n", worst))
cat("  so both eigenvalues A_0 +/- A_I stay non-negative with no condition imposed:\n")
cat("  the seam parametrisation cannot produce a negative noise eigenvalue.\n")

cat("\n=== 5. inverting it: a measured ratio constrains kappa to a branch, not a point\n\n")
# t is invariant under kappa -> 1/kappa but r changes sign, so q has no such
# degeneracy. It is, however, NOT monotone: it turns over twice, so a single
# measured ratio leaves up to three kappa consistent with it.
qof <- function(kap, C) { r <- (1-kap^2)/(1+kap^2); tt <- 2*kap/(1+kap^2)
                          tt*(C+r)/(1+r*C) }
gr <- exp(seq(log(1e-3), log(1e3), length.out = 400001))
for (a in c(0.01, 0.1)) {
  C <- Cfun(a); qq <- qof(gr, C)
  turn <- which(diff(sign(diff(qq))) != 0) + 1
  cat(sprintf("  C = %.4f:  q runs %.4f to %.4f, turning at kappa = %s\n",
              C, min(qq), max(qq), paste(sprintf("%.4f", gr[turn]), collapse = ", ")))
  for (kt in c(0.3, 0.7, 1.0, 1.5, 3.0)) {
    tgt <- qof(kt, C)
    roots <- gr[which(abs(diff(sign(qq - tgt))) > 0)]
    cat(sprintf("    kappa_true %5.2f -> q = %+8.5f  consistent kappa: %s\n",
                kt, tgt, paste(sprintf("%.4f", roots), collapse = ", ")))
  }
}
cat("\n  So one even/odd ratio narrows the seam to at most three values and does not\n")
cat("  pick among them. Two port geometries with different C intersect to one:\n")
C1 <- Cfun(0.01); C2 <- Cfun(0.3); kt <- 3.0
r1 <- gr[which(abs(diff(sign(qof(gr,C1) - qof(kt,C1)))) > 0)]
r2 <- gr[which(abs(diff(sign(qof(gr,C2) - qof(kt,C2)))) > 0)]
both <- r1[sapply(r1, function(x) any(abs(r2 - x) < 1e-2*x))]
cat(sprintf("    C=%.4f alone: %s\n    C=%.4f alone: %s\n    both:         %s   (true %.2f)\n",
    C1, paste(sprintf("%.3f", r1), collapse=", "), C2, paste(sprintf("%.3f", r2), collapse=", "),
    paste(sprintf("%.3f", both), collapse=", "), kt))
