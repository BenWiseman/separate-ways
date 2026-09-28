#!/usr/bin/env Rscript
# The "in" vacuum is not a rival state that happens to give a different number. Theta-invariance
# forbids it, and says where.
#
# WHY. Appendix B.1 excludes the "in" vacuum by comparing production integrals, 0.0448968 against
# 0.0127597, a factor 3.52 and 0.605 in the mass. That is true and it is the wrong argument to lead
# with, because it invites the reply that two numbers differing is not an exclusion. The fold
# excludes it outright. The Landau-Zener sweep from eta = -infinity leaves occupation |beta|^2 =
# exp(-x^2), which at x = 0 is ONE: a mode of zero momentum flips completely, the sweep being
# perfectly diabatic there. Theta-invariance forces exactly one half at x = 0, where the band has
# zero width. So the in vacuum misses the only point the fold pins exactly, and by the largest
# margin available.
#
# WHERE IT RE-ENTERS is exact. The band's upper edge is n_max = 1 - n_min, and setting
# exp(-x^2) = n_max with u = exp(-x^2) gives 2u - 1 = sqrt(1-u), so 4u^2 - 3u = 0 and u = 3/4:
# the in vacuum lies outside the band for x < sqrt(log(4/3)) = 0.5364 and inside beyond it. That
# is an algebraic number, not a fit, and it is the statement B.1 should make.

TOL <- 1e-12
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

nmin <- function(x) { u <- exp(-x^2); u/(2*(1 + sqrt(1 - u))) }
nmax <- function(x) 1 - nmin(x)
nin  <- function(x) exp(-x^2)

cat("=== 1. at zero momentum the in vacuum is at one and the fold demands one half ===\n")
cat(sprintf("   n_in(0)  = %.12f\n", nin(0)))
cat(sprintf("   n_min(0) = %.12f   n_max(0) = %.12f   band width = %.2e\n",
            nmin(0), nmax(0), nmax(0) - nmin(0)))
note(abs(nin(0) - 1) < TOL, "the in vacuum saturates at zero momentum")
note(abs(nmin(0) - 0.5) < 1e-8 && (nmax(0) - nmin(0)) < 1e-8, "the band closes at one half there")
cat("   So the miss at x = 0 is 1/2, which is the whole of the available range.\n")

cat("\n=== 2. and it re-enters the band at an algebraic point ===\n")
xstar <- sqrt(log(4/3))
cat(sprintf("   closed form  sqrt(log(4/3)) = %.10f\n", xstar))
r <- uniroot(function(x) nin(x) - nmax(x), c(0.1, 2), tol = 1e-14)$root
cat(sprintf("   root of n_in - n_max        = %.10f\n", r))
note(abs(r - xstar) < 1e-9, "the crossing is at sqrt(log(4/3)) exactly")
cat(sprintf("   and there exp(-x^2) = %.10f, which is 3/4\n", exp(-xstar^2)))
note(abs(exp(-xstar^2) - 0.75) < TOL, "the crossing is where the occupation is three quarters")
cat("        x      n_in      n_max     inside the band?\n")
for (x in c(0.1, 0.3, 0.5, xstar, 0.7, 1.2)) {
  cat(sprintf("   %6.4f %10.6f %10.6f %14s\n", x, nin(x), nmax(x),
              ifelse(nin(x) <= nmax(x) + 1e-12, "yes", "no")))
  note((x >= xstar - 1e-9) == (nin(x) <= nmax(x) + 1e-12), sprintf("the crossing is at x* (x = %g)", x))
}

cat("\n=== 3. B.1's own numbers, reproduced ===\n")
Iof <- function(f) integrate(function(x) x^2*f(x), 0, Inf, rel.tol = 1e-13)$value/pi^2
I0 <- Iof(nmin); Ib <- Iof(nin)
cat(sprintf("   I_b = 1/(4 pi^{3/2}) = %.7f  quadrature %.7f\n", 1/(4*pi^1.5), Ib))
cat(sprintf("   I   = %.7f,  ratio %.2f,  mass ratio %.3f\n", I0, Ib/I0, (Ib/I0)^(-2/5)))
note(abs(Ib - 1/(4*pi^1.5)) < 1e-10, "the in vacuum's integral has B.1's closed form")
note(abs(Ib/I0 - 3.52) < 5e-3 && abs((Ib/I0)^(-2/5) - 0.605) < 5e-4, "B.1's factor and mass ratio")

cat("\n=== 4. the plant: the test must not pass a state that IS in the band ===\n")
for (nm in c("half-angle floor", "band centre", "the in vacuum")) {
  f <- switch(nm, "half-angle floor" = nmin, "band centre" = function(x) rep(0.5, length(x)), nin)
  XG <- exp(seq(log(1e-4), log(8), length.out = 20000))
  out <- any(f(XG) > nmax(XG) + 1e-12)
  cat(sprintf("   %-18s leaves the band: %s\n", nm, ifelse(out, "yes", "no")))
  note(out == (nm == "the in vacuum"), sprintf("only the in vacuum leaves the band (%s)", nm))
}

cat("\n=== 5. the numbers the manuscripts quote, for the digit checker ===\n")
cat(sprintf("   the crossing, sqrt(log(4/3)):  %.4f\n", xstar))
cat(sprintf("   the occupation there:          %.2f\n", exp(-xstar^2)))
cat(sprintf("   the in vacuum at zero momentum: %.0f\n", nin(0)))

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
