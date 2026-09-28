#!/usr/bin/env Rscript
# What the residual phase can and cannot be made to pay, computed rather than guessed.
#
# WHY. Section 2.3's mass is a ceiling and section 4.1 leaves the phase free on a band, so the
# natural question is whether the two forced conditions, Theta-invariance and the Hadamard tail,
# floor the mass as well as cap it. A floor would turn the bound into a window and put KM3NeT's
# 220 PeV median inside it rather than merely under it. They do not, and the reason is worth
# having explicitly, because the near miss is instructive: the state's Gaussian TAIL is not the
# same constraint as the state's SCALE, and only a scale would bound the excursion.
#
# WHAT WAS TRIED AND FAILED. Take the widest excursion weighted by the state's own Gaussian,
#   n_A = n_min + A exp(-x^2) (n_max - n_min),  0 <= A <= 1,
# which sits in the band, carries no scale the adopted state does not already have, and keeps the
# tail exponent exactly. It returns I/I_min = 3.75 and M1 = 289.7 PeV at A = 1, and it looked like
# a floor. It is not one. A polynomial prefactor keeps the tail exponent while moving the
# excursion's peak to any momentum, so
#   w_k(x) = (x^2/k)^k exp(k - x^2),   which peaks at exactly 1 at x = sqrt k,
# is admissible for every k and drives I up without bound. Section 2.3's Table 1 already says this
# in its own variable, since the cutoff x_c is a scale the adopted state does not have, and it is
# reproduced below as a cross-check. So the mass has no floor, and the 289.7 PeV figure is a
# benchmark for the one-scale family and nothing more.

TOL <- 1e-10
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

nmin  <- function(x) { u <- exp(-x^2); u/(2*(1 + sqrt(1 - u))) }   # (1-sqrt(1-u))/2 without the
nmax  <- function(x) 1 - nmin(x)                                    # cancellation past x = 6
width <- function(x) nmax(x) - nmin(x)
Iof   <- function(f, lo = 0, hi = Inf) integrate(function(x) x^2*f(x), lo, hi,
                                                 rel.tol = 1e-13, subdivisions = 6000L)$value/pi^2
I0    <- Iof(nmin)
M1    <- function(I) 491.6 * (I/I0)^(-2/5)
   # 491.6 is section 2.3's endpoint, taken as the NORMALISATION of this ratio and not
   # computed here, so no claim row may cite it against this script. What is computed is
   # the ratio I/I0 and everything that follows from it.

cat("=== 1. the adopted state ===\n")
cat(sprintf("   I_min = %.13f, residual against the manuscript %.1e,  M1 = %.1f PeV\n",
            I0, abs(I0 - 0.0127596673634), M1(I0)))   # see the note in
                                                      # the_band_excludes_a_thermal_profile.R:
                                                      # reference digits stay out of the output
note(abs(I0 - 0.0127596673634) < 1e-11, "the half-angle production integral reproduces")

cat("\n=== 2. section 2.3's Table 1, reproduced from the band alone ===\n")
Istep <- function(xc) Iof(nmax, 0, xc) + Iof(nmin, xc, Inf)
tab <- c(1.000, 1.008, 1.119, 1.574, 2.700, 8.312)
mas <- c(491.6, 490.1, 470.0, 410.0, 330.4, 210.7)
cat("        x_c     I/I_min        M1     residual vs the table\n")
for (i in seq_along(tab)) {
  xc <- c(0.10, 0.25, 0.50, 0.75, 1.00, 1.50)[i]; r <- Istep(xc)/I0
  cat(sprintf("   %8.2f %10.3f %9.1f %17.1e %8.1e\n", xc, r, M1(Istep(xc)),
              abs(r - tab[i]), abs(M1(Istep(xc)) - mas[i])))
  note(abs(r - tab[i]) < 5e-4 && abs(M1(Istep(xc)) - mas[i]) < 0.1, sprintf("Table 1 row x_c = %g", xc))
}

cat("\n=== 3. the one-scale family, which is the benchmark and was mistaken for a floor ===\n")
nA <- function(x, A, w = function(z) exp(-z^2)) nmin(x) + A*w(x)*width(x)
J  <- Iof(function(x) exp(-x^2)*width(x))
cat("         A        I/I_min      M1 (PeV)    E_nu (PeV)    inside the band?\n")
XG <- exp(seq(log(1e-6), log(25), length.out = 30000))
for (A in c(0, 0.25, 0.5, 1, 1.05)) {
  I <- I0 + A*J; ins <- min(nmax(XG) - nA(XG, A)) >= -1e-15
  cat(sprintf("   %8.2f %12.4f %12.1f %13.1f %18s\n", A, I/I0, M1(I), M1(I)/2, ifelse(ins, "yes", "no")))
  note((A <= 1 + 1e-12) == ins, sprintf("A = 1 is the band's edge (A = %g)", A))
}
Iedge <- I0 + J
cat(sprintf("   At A = 1: I/I_min = %.4f, M1 = %.1f PeV, and Table 1 reaches the same ratio at\n",
            Iedge/I0, M1(Iedge)))
xc_eq <- uniroot(function(xc) Istep(xc)/I0 - Iedge/I0, c(0.5, 2), tol = 1e-12)$root
cat(sprintf("   x_c = %.4f, so the two ways of spending the freedom agree where they cross.\n", xc_eq))
note(abs(xc_eq - 1.1397) < 1e-3, "the one-scale edge sits at x_c = 1.14 on Table 1")

cat("\n=== 4. and here is why it is not a floor ===\n")
cat("   The Hadamard condition fixes the tail EXPONENT, not the excursion's position. Weight the\n")
cat("   band by w_k(x) = (x^2/k)^k exp(k - x^2), which peaks at exactly one at x = sqrt k and has\n")
cat("   d log w_k / d(x^2) = k/x^2 - 1, so every member has the adopted tail. I grows without end:\n")
lwk <- function(k) function(x) if (k == 0) -x^2 else k*log(x^2/k) + k - x^2
wk  <- function(k) function(x) exp(lwk(k)(x))
cat("          k      peak at x   max of w_k    I/I_min       M1 (PeV)     tail exponent\n")
for (k in c(0, 1, 4, 9, 16, 25)) {
  w <- wk(k); pk <- max(w(XG))
  I <- I0 + Iof(function(x) w(x)*width(x))
  lw <- lwk(k); ex <- (lw(140) - lw(100))/(140^2 - 100^2)  # from the log, since w underflows;
                                                          # far out, where k/x^2 is negligible
  cat(sprintf("   %8d %11.2f %13.6f %11.3f %13.1f %17.4f\n", k, sqrt(max(k, 0)), pk, I/I0, M1(I), ex))
  note(pk <= 1 + 1e-9, sprintf("w_%d never leaves the band", k))
  note(abs(ex + 1) < 2e-3, sprintf("w_%d keeps the adopted tail exponent", k))
}
cat("   So the production integral is unbounded above inside the admissible class, and the mass\n")
cat("   has no floor. The plant: a weight that does NOT keep the exponent must be visible here.\n")
lbad  <- function(x) -x^2/2                      # written as a log: exp(-x^2/2) underflows
exbad <- (lbad(140) - lbad(100))/(140^2 - 100^2)
cat(sprintf("      exp(-x^2/2) gives exponent %.4f against the required %.4f\n", exbad, -1))
note(abs(exbad + 0.5) < 1e-9, "plant: a slower Gaussian is caught by the exponent test")
note(Iof(function(x) wk(25)(x)*width(x))/I0 > 10, "plant: the k = 25 member really does exceed the benchmark")

cat("\n=== 5. what survives ===\n")
cat("   The ceiling, untouched: every admissible state has n >= n_min, so I >= I_min and\n")
cat("   M1 <= 491.6 PeV, and that is an operator inequality rather than a class statement.\n")
cat("   The benchmark: among occupations carrying no scale beyond the state's own Gaussian, the\n")
cat("   widest returns M1 = 289.7 PeV.\n")
cat("   What would give a floor is a bound on the state's SCALE, which nothing here supplies, and\n")
cat("   naming that is more use to a reader than a window that does not hold.\n")

cat("\n=== 6. the numbers the manuscripts quote, for the digit checker ===\n")
cat(sprintf("   the one-scale benchmark, I over I_min:   %.2f\n", Iedge/I0))
cat(sprintf("   the one-scale benchmark, M1 in PeV:      %.1f\n", M1(Iedge)))
cat(sprintf("   where it sits on Table 1, x_c:           %.2f\n", xc_eq))

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
