#!/usr/bin/env Rscript
# half_period_match.R -- Section 3 and A.10 say the horizon temperature is returned rather
# than assumed: requiring the geometric image kernel to be a thermal half-period shift of
# the Bunch-Davies kernel fixes beta. Two numbers carry that claim, "the match holds to
# 4e-17 at beta = 2 pi/H" and "fails by 0.115 ten per cent away on either side, fifteen
# orders worse", and the second had no script. `tangents/seam/alpha_is_kms.R` checks
# the match AT beta and sweeps beta for the alpha^2 = 1 condition, but never sweeps beta in
# the half-shift comparison itself, which is where 0.115 comes from. Written here so the
# companion cites a file in its own tree.
#
# Setting: de Sitter conformal scalar on the central worldline, H = 1, where the
# Bunch-Davies kernel is W(t) = csch^2(Ht/2) up to a constant and A.6's image kernel is
# -sech^2(Ht/2) up to the same constant. The claim is  W(t - i beta/2) = -sech^2(Ht/2).

H <- 1; beta <- 2*pi/H
W  <- function(t) 1/sinh(H*t/2)^2
img <- function(t) -1/cosh(H*t/2)^2
ts <- c(0.4, 0.9, 1.7, 3.1)
mism <- function(b) max(sapply(ts, function(t) Mod(W(t - 1i*b/2) - img(t))))

cat("=== 1. the match at beta = 2 pi / H, term by term ===\n\n")
cat("        t        W(t - i beta/2)              -sech^2(Ht/2)        |difference|\n")
for (t in ts) { a <- W(t - 1i*beta/2); b <- img(t)
  cat(sprintf("   %7.2f  %10.6f%+10.6fi   %10.6f%+10.6fi   %10.2e\n",
              t, Re(a), Im(a), Re(b), Im(b), Mod(a-b))) }
cat(sprintf("\n   worst over the four times: %.4e, which is the 4e-17 the paper quotes.\n", mism(beta)))
stopifnot(mism(beta) < 1e-16)
cat("   The half shift turns csch^2 into -sech^2 exactly. Nothing thermal was assumed to\n")
cat("   get there: A.6 derives the image kernel from the antipodal map's geometry.\n")

cat("\n=== 2. sweeping beta, which is where the constraint lives ===\n\n")
cat("      beta/(2 pi/H)     worst mismatch      ratio to the value at 2 pi/H\n")
base <- mism(beta)
for (fr in c(0.80, 0.90, 0.95, 0.99, 1.00, 1.01, 1.05, 1.10, 1.20)) {
  v <- mism(fr*beta)
  cat(sprintf("   %14.3f  %17.4e  %27s\n", fr, v,
              if (fr == 1) "-- the match" else sprintf("%.2e", v/base))) }
cat(sprintf("\n   At ten per cent away on either side the mismatch is %.4f, and %.4f/%.2e is\n",
            mism(0.9*beta), mism(0.9*beta), base))
cat(sprintf("   %.2e, so 'fifteen orders worse' is right and symmetric in the two directions:\n",
            mism(0.9*beta)/base))
cat(sprintf("   %.2e below and %.2e above.\n", mism(0.9*beta)/base, mism(1.1*beta)/base))
stopifnot(abs(mism(0.9*beta) - 0.115) < 0.002, mism(0.9*beta)/base > 1e15)

cat("\n=== 3. the odd multiples, which the match does not separate ===\n")
cat("   sinh(z - i(2j+1)pi/2) = (-1)^{j+1} i cosh z, and squaring removes the sign, so the\n")
cat("   match holds at every odd multiple of the fundamental. Checked:\n\n")
cat("      beta/(2 pi/H)      worst mismatch\n")
for (k in c(1, 3, 5, 7)) cat(sprintf("   %14d  %17.2e\n", k, mism(k*beta)))
cat("\n   and fails at the even ones, which is what makes the fundamental the primitive\n")
cat("   period rather than an accident of the search window:\n\n")
for (k in c(2, 4, 6)) cat(sprintf("   %14d  %17.2e\n", k, mism(k*beta)))
for (k in c(1,3,5,7)) stopifnot(mism(k*beta) < 1e-14)
for (k in c(2,4,6)) stopifnot(mism(k*beta) > 1e-3)

cat("\n=== 3b. the free minimisation, for the value A.10 quotes ===\n")
cat("   Treating beta as unknown and minimising the mismatch, with no knowledge of 2 pi/H:\n\n")
o1 <- optimize(function(b) mism(b), c(0.5*beta, 1.5*beta), tol = 1e-10)
o2 <- optimize(function(b) mism(b), c(0.5*beta, 1.5*beta))
cat(sprintf("     tight tolerance:    %.10f\n", o1$minimum))
cat(sprintf("     default tolerance:  %.5f\n", o2$minimum))
cat(sprintf("     2 pi / H:           %.10f\n", beta))
cat("   Both return 2 pi / H. Note what this does NOT reproduce: A.10 also quotes 6.28320\n")
cat("   from a free minimisation, and that figure belongs to a different quantity, the\n")
cat("   alpha^2 = 1 condition minimised in `tangents/seam/alpha_is_kms.R`, which\n")
cat("   returns 6.2832000828. Minimising the half-shift MATCH, as here, lands on\n")
cat("   6.2831853072. Two minimisations of two conditions, both at the same temperature.\n")
stopifnot(abs(o1$minimum - beta) < 1e-6)

cat("\n=== 4. the check has to be able to fail ===\n")
cat("   Replace the image kernel by things it is not and the match at beta must break:\n\n")
for (nm in c("+sech^2", "-sech", "-sech^2 scaled by 1.001")) {
  bad <- switch(nm, "+sech^2" = function(t) 1/cosh(H*t/2)^2,
                    "-sech"   = function(t) -1/cosh(H*t/2),
                    "-sech^2 scaled by 1.001" = function(t) -1.001/cosh(H*t/2)^2)
  v <- max(sapply(ts, function(t) Mod(W(t - 1i*beta/2) - bad(t))))
  cat(sprintf("     %-26s worst mismatch at beta: %.4f   <- fails, as it must\n", nm, v))
  stopifnot(v > 1e-4) }

cat("
=== flatly ===

  Both numbers reproduce. The half-period shift reproduces A.6's image kernel to 4.42e-17
  at beta = 2 pi/H, and ten per cent away in either direction the mismatch is 0.1146, a
  factor 2.59e15 worse, which is the 'fifteen orders' the paper claims and is symmetric
  in the two directions. The match survives every odd multiple of the fundamental and
  fails at every even one, so primitivity is what picks 2 pi/H and the search window is
  not doing the work.

  What this does NOT show, and the paper says so in both places: that the match is a
  derivation of the temperature rather than a consistency check, since csch^2 has period
  2 pi i / H however it was obtained. It shows the half, which is that the fold's map is
  the square root of the thermal transformation.\n")
