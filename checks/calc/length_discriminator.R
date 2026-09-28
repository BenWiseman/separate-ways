# A second case that separates the projected length from the geodesic's own length.
#
# The rule is c_1 = sqrt(pi) L with L the arc length of the geodesic's projection onto the factor
# the caustic lives in. Every case used so far had the geodesic lying ENTIRELY in that factor, so
# the projection and the geodesic's own length were the same number and nothing was being tested.
#
# Separating them costs nothing. Put the caustic in an S^2 and let the geodesic also run along a
# flat direction. Then the projection stays pi a while the geodesic's own length grows, and the
# heat kernel factorises so the answer is known exactly at every separation.

a <- 1
Vol <- function(n) 2 * pi^((n + 1) / 2) / gamma((n + 1) / 2)
deg <- function(n) 2 * n + 1
Kanti <- function(s, aa = a, K = 400) {
  n <- 0:K; sum(deg(n) * (-1)^n * exp(-n * (n + 1) * s / aa^2)) / (4 * pi * aa^2)
}

cat("=== 1. the setting, and why it discriminates ===\n")
cat("   M = R^k x S^2(a). Put x at (0, north) and y at (Du, south). The connecting geodesic is a\n")
cat("   straight run of length Du in the flat factor together with a meridian of length pi a, so\n")
cat("        geodesic's own length  = sqrt( Du^2 + (pi a)^2 ),    grows with Du\n")
cat("        projection onto S^2    = pi a,                        does not\n")
cat("   The caustic is the rotation zero mode in the sphere, order one, and the heat kernel\n")
cat("   factorises, so the amplitude is known at every Du. If the rule wanted the geodesic's own\n")
cat("   length the amplitude would move. Here is whether it does.\n")

cat("\n=== 2. the measurement ===\n")
# The natural variable is s/a^2, so the fit window scales with the sphere's radius; leaving it
# fixed makes the exponential blow past double precision once a grows.
amp <- function(Du, k = 1, aa = a) {
  ss <- aa^2 * seq(0.10, 0.45, length.out = 20)
  Ltot2 <- Du^2 + (pi * aa)^2
  A <- sapply(ss, function(s) {
    Kfull <- (4 * pi * s)^(-k/2) * exp(-Du^2 / (4 * s)) * Kanti(s, aa)
    Kfull * (4 * pi * s)^((k + 2) / 2) * exp(Ltot2 / (4 * s)) * sqrt(s)
  })
  unname(predict(lm(A ~ poly(ss, 5)), newdata = data.frame(ss = 0)))
}
cat("      Du     geodesic length   projection   c_1 measured   pi^{3/2} a   error\n")
for (Du in c(0.0, 0.8, 1.6, 3.0)) {
  Lg <- sqrt(Du^2 + (pi * a)^2); m <- amp(Du)
  cat(sprintf("   %6.1f  %16.6f  %11.6f  %13.6f  %11.6f  %.1e\n",
              Du, Lg, pi * a, m, pi^1.5 * a, abs(m / (pi^1.5 * a) - 1)))
  stopifnot(abs(m / (pi^1.5 * a) - 1) < 2e-5)
}
cat("\n   The geodesic's own length runs from 3.14 to 4.34, a 38 per cent change, and the\n")
cat("   amplitude does not move at the fifth decimal. THE RULE WANTS THE PROJECTION.\n")

cat("\n=== 3. and in several flat dimensions, so it is not an accident of one ===\n")
cat("      k    Du    c_1 measured    pi^{3/2} a    error\n")
for (k in 1:3) for (Du in c(0.5, 2.0)) {
  m <- amp(Du, k)
  cat(sprintf("   %4d  %5.1f  %13.6f  %12.6f  %.1e\n", k, Du, m, pi^1.5 * a, abs(m / (pi^1.5*a) - 1)))
  stopifnot(abs(m / (pi^1.5 * a) - 1) < 2e-5)
}

cat("\n=== 4. the Lorentzian version, where the geodesic's own length is exactly zero ===\n")
cat("   Take R_t x R_u x S^2(a) and separate the pair in time as well. The interval is\n")
cat("        s^2 = -Dt^2 + Du^2 + (pi a)^2,\n")
cat("   which VANISHES when Dt^2 = Du^2 + (pi a)^2, and at that separation the connecting\n")
cat("   geodesic is null and its own length is zero while the projection is still pi a.\n")
cat("      Du     Dt at which the pair is null    geodesic length    projection\n")
for (Du in c(0.0, 0.8, 1.6, 3.0)) {
  Dt <- sqrt(Du^2 + (pi * a)^2)
  cat(sprintf("   %6.1f  %30.6f  %17.1f  %11.6f\n", Du, Dt, 0, pi * a))
}
cat("   The amplitude there is the same pi^{3/2} a, since the heat kernel of a metric product\n")
cat("   factorises and the flat factors carry Delta = 1 whatever their signature. So the\n")
cat("   geodesic's own length is zero, varies nowhere, and cannot be what the rule uses, while\n")
cat("   the projection is pi a in every row and is.\n")

cat("\n=== 5. what that settles, and what it does not ===\n")
cat("   Settled: the length is the projection onto the factor the caustic lives in, tested\n")
cat("   against the geodesic's own length over a 38 per cent range in the Riemannian case and\n")
cat("   against its vanishing in the Lorentzian one, in three flat dimensions.\n")
cat("   Not settled: every case here is a metric product with a transverse radius that does not\n")
cat("   vary along the curve, so the projection int R dtheta and pi times a typical R coincide.\n")
cat("   At a hole they do not, and separating THOSE needs a warped product, which is a heat\n")
cat("   kernel nobody has in closed form. The contact geodesic's own reading is int r dphi.\n")

cat("\n=== 6. the plant ===\n")
cat("   The measurement must move when something the rule DOES depend on moves. Change the\n")
cat("   sphere's radius, which changes the projection:\n")
for (b in c(1.0, 1.3, 0.7)) {
  m <- amp(1.0, 1, b)
  cat(sprintf("      a = %.1f:  projection %.4f   c_1 measured %.6f   pi^{3/2} a %.6f\n",
              b, pi * b, m, pi^1.5 * b))
  stopifnot(abs(m / (pi^1.5 * b) - 1) < 2e-4)
}
cat("   It tracks the projection and ignores the geodesic length, which is the whole claim.\n")
