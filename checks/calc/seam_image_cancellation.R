#!/usr/bin/env Rscript
# Does a constant seam reflectivity put a singularity at spacelike separation? Two passages of
# A.11 gave opposite answers, so this settles which.
#
# THE DISPUTE. A.11's opening kills a regularity route to kappa = 1: the image term
# W_0(x-y) + r W_0(x - ybar) has same-side support, the coincidence x = ybar is unreachable for two
# points on the same side, and across the seam the two incident channels cancel the image
# coefficient, rt + t(-r) = 0. A later passage in the same appendix revives the route and calls it
# load-bearing. Both cannot stand. This computes the coefficient from the mode sum rather than
# quoting either.
#
# WHAT IT FINDS. The opening is right about the ordinary defect and for a reason that is forced
# rather than chosen: unitarity of a symmetric two-port requires the reflectivity seen from the far
# side to be MINUS the near one, and that is exactly what makes the cross-seam image coefficient
# vanish. It is right at every kappa, so no value of kappa is selected.
#
# AND WHAT IT ALSO FINDS, which is why the later passage is not simply a slip. On a Z_2 quotient the
# mode space is projected, the two channels are no longer summed independently, and the cross-seam
# coefficient becomes t^2 - r^2, which is nonzero at kappa = 1 and vanishes at |r| = |t|. So on the
# quotient there IS an image term across the seam at the transparent point, and it is the fold's own
# cross-sheet correlator, the object Section 3 builds the thermal law on. Regularity cannot be used
# against it without discarding the construction, and its singularity is the caustic A.18 computes.
# Either way the route does not select kappa.

TOL <- 1e-12
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

rk <- function(k) (1 - k^2)/(1 + k^2)
tk <- function(k) 2*k/(1 + k^2)

cat("=== 1. unitarity fixes the far-side reflectivity to minus the near-side one ===\n")
cat("   For S = [[r, t],[t, rp]] with real entries, unitarity is r^2 + t^2 = 1 and r t + t rp = 0.\n")
cat("   The second forces rp = -r whenever t is nonzero, so it is not a sign convention.\n\n")
cat("        kappa        r            t         r^2+t^2      rp from unitarity     -r\n")
for (k in c(0.25, 0.5, 1, 2, 4)) {
  r <- rk(k); t <- tk(k)
  rp <- if (abs(t) > 1e-14) -r*t/t else NA
  cat(sprintf("   %10.3f %11.6f %12.6f %12.8f %20.6f %10.6f\n", k, r, t, r^2 + t^2, rp, -r))
  note(abs(r^2 + t^2 - 1) < TOL, "the matching is unitary")
  note(abs(rp + r) < TOL, "unitarity gives rp = -r")
}

cat("\n=== 2. the ordinary defect: both channels summed, and the image term cancels ===\n")
cat("   Channel 1 is incident from the left, channel 2 from the right. For x left and y right the\n")
cat("   products carry an exp(-ik(x+y)) piece with coefficient r t* from the first channel and\n")
cat("   t rp* from the second, so the cross-seam image coefficient is r t + t rp.\n")
cat("   Computed from the mode functions themselves rather than from that sentence: the coefficient\n")
cat("   is read off by projecting the summed product onto exp(-ik(x+y)) over a grid in x and y.\n")
# The coefficient is an algebraic property of the mode functions, so it is extracted by projecting
# onto exp(-ik(x+y)) over a FULL PERIOD in each variable, where the four exponentials appearing are
# orthogonal. Projecting over the physical half-lines instead leaks between them: the first version
# of this did that and returned numbers that matched nothing, which is what the section-4 assertion
# caught.
KK <- 1.7
proj <- function(P, X, Y) mean(P * exp(1i*KK*(X + Y)))
grid <- function(n = 512) {
  L <- 2*pi/KK; u <- (0:(n-1))*L/n
  list(X = outer(u, u, function(a, b) a), Y = outer(u, u, function(a, b) b))
}
crossco <- function(k, projected = FALSE, rpsign = -1) {
  r <- rk(k); t <- tk(k); rp <- rpsign*r
  g <- grid(); X <- g$X; Y <- g$Y
  if (!projected) {
    P <- (exp(1i*KK*X) + r*exp(-1i*KK*X)) * Conj(t*exp(1i*KK*Y)) +
         (t*exp(-1i*KK*X)) * Conj(exp(-1i*KK*Y) + rp*exp(1i*KK*Y))
  } else {
    ex <- (exp(1i*KK*X) + (r + t)*exp(-1i*KK*X))/sqrt(2)
    ey <- ((t + rp)*exp(1i*KK*Y) + exp(-1i*KK*Y))/sqrt(2)
    P  <- ex * Conj(ey)
  }
  proj(P, X, Y)
}
cat("\n        kappa      r t + t rp      coefficient read off the modes\n")
for (k in c(0.25, 0.5, 1, 2, 4)) {
  r <- rk(k); t <- tk(k)
  alg <- r*t + t*(-r); num <- crossco(k)
  cat(sprintf("   %10.3f %14.8f %30s\n", k, alg,
              sprintf("%+.2e%+.2ei", Re(num), Im(num))))
  note(abs(alg) < TOL && Mod(num) < 1e-8, "the cross-seam image coefficient vanishes")
}
cat("   It vanishes at every kappa, so no kappa is selected. The route is empty for the ordinary\n")
cat("   defect, which is what A.11's opening says.\n")
cat("\n   And the plant, because a cancellation that holds for every input proves nothing unless the\n")
cat("   input can be wrong: rp = +r instead, which violates unitarity, must leave 2rt behind.\n")
crossbad <- function(k) crossco(k, rpsign = +1)
for (k in c(0.5, 2)) {
  r <- rk(k); t <- tk(k); b <- crossbad(k)
  cat(sprintf("      kappa = %.2f: rp = +r gives %+.6f against 2rt = %+.6f\n",
              k, Re(b), 2*r*t))
  note(abs(Re(b) - 2*r*t) < 1e-6 && abs(2*r*t) > 1e-3, "the non-unitary choice leaves 2rt, as it must")
}

cat("\n=== 3. the same-side coincidence is out of reach, independently of the coefficient ===\n")
cat("   The reflection term's singular locus is x + y = 0. With both points on the same side that\n")
cat("   sum has a fixed sign and never reaches zero, so the term cannot bite where it exists.\n")
cat("      two points on the left      x + y      reaches 0?\n")
for (p in list(c(-0.2,-0.2), c(-3,-0.001), c(-1e-9,-1e-9))) {
  s <- p[1] + p[2]
  cat(sprintf("   %12.2e %12.2e %10.2e %12s\n", p[1], p[2], s, ifelse(abs(s) < 1e-300, "YES", "no")))
  note(s < 0, "the same-side sum stays negative")
}

cat("\n=== 4. the quotient, where the projection changes the answer ===\n")
cat("   A Z_2 quotient keeps a projected combination rather than summing two independent channels.\n")
cat("   Then the cross-seam coefficient is (r + t)(t - r)/2 = (t^2 - r^2)/2, which does NOT vanish\n")
cat("   at kappa = 1 and vanishes instead where |r| = |t|:\n\n")
cat("        kappa    (t^2 - r^2)/2      coefficient read off the projected modes\n")
for (k in c(0.25, 0.4142136, 1, 2.4142136, 4)) {
  r <- rk(k); t <- tk(k); num <- crossco(k, projected = TRUE)
  cat(sprintf("   %10.4f %14.8f %36s\n", k, (t^2 - r^2)/2, sprintf("%+.6f%+.6fi", Re(num), Im(num))))
  note(Mod(num - (t^2 - r^2)/2) < 1e-8, "the projected coefficient is (t^2 - r^2)/2")
}
note(abs((tk(1)^2 - rk(1)^2)/2 - 0.5) < TOL, "at kappa = 1 the projected coefficient is 1/2, not 0")
cat("\n   So on the quotient the transparent seam carries the LARGEST cross-seam image term, not the\n")
cat("   smallest. That term is the fold's own cross-sheet correlator, which Section 3 needs and\n")
cat("   whose singularity is the caustic A.18 computes. A regularity condition cannot be turned\n")
cat("   against it without discarding the construction it supports.\n")

cat("\n=== 4b. the numbers the manuscripts quote, as magnitudes ===\n")
cat("   The provenance checker compares digits and drops signs, so these are printed unsigned; the\n")
cat("   signs are asserted in sections 2 and 4 and this file stops if either comes out wrong.\n")
cat(sprintf("   projected coefficient at kappa = 1      %.4f\n", abs((tk(1)^2 - rk(1)^2)/2)))
cat(sprintf("   projected coefficient at kappa = 1/4    %.4f\n", abs((tk(0.25)^2 - rk(0.25)^2)/2)))

cat("\n=== 5. the verdict, stated once ===\n")
cat("   Regularity does not select kappa. On the ordinary defect the image coefficient vanishes at\n")
cat("   every kappa by unitarity; on the quotient it is largest exactly at the transparent point and\n")
cat("   is the object the paper wants. What remains behind kappa = 1 is the corner term's uniqueness\n")
cat("   among scale-free seams, together with the physical requirement that an interface reflecting\n")
cat("   with the same amplitude at every frequency is not an interface. That second step is a\n")
cat("   physical expectation and not a theorem, and the manuscripts now say so.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
