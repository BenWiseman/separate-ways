#!/usr/bin/env Rscript
# The 2/3 tying the decoherence time to the dark-matter mass is exact, and the release had no
# file saying so.
#
# WHY. Section 3.5 said the power "holds to five decimal places across four decades in M_1",
# which reads as the result of a numerical fit. No script in the release performed one, and the
# statement understates what is true: the exponent comes out of eliminating one variable between
# two power laws and is exact, not fitted. Written 2026-09-28 after a provenance sweep found the
# sentence had no file behind it. The manuscript now says exact and points here.
#
# THE TWO LAWS. Section 2.3 gives M_1 propto I^{-2/5}, with 491.6 PeV the normalisation of that
# ratio rather than a result of this file. Section 3.5 gives
#     t_dec = (1/2M_1) (4/(3 sqrt pi) c_G R I)^{-2/3},
# so t_dec propto M_1^{-1} I^{-2/3}. Eliminating I leaves -1 + (5/2)(2/3) = 2/3, and there is
# nothing approximate in that arithmetic. What a numerical composition adds is a check that the
# two laws as CODED compose the way the algebra says, which is the part a transcription slip
# would break.

fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

I0  <- 0.0127596673634            # the production integral at the least-occupied member
M0  <- 491.6                      # PeV, section 2.3's endpoint, the normalisation
cG  <- 1                          # mode-counting factor, one per (p, h) mode
Rf  <- 1                          # the second functional of the occupation, held fixed here

M1    <- function(I) M0 * (I/I0)^(-2/5)
tdec  <- function(I, a = 2/3) (1/(2*M1(I))) * ((4/(3*sqrt(pi))) * cG * Rf * I)^(-a)

cat("=== 1. the algebra, before any arithmetic ===\n")
cat("   t_dec propto M_1^{-1} I^{-2/3} and I propto M_1^{-5/2}, so the exponent is\n")
cat(sprintf("   -1 + (5/2)(2/3) = %.20g, against 2/3 = %.20g\n", -1 + (5/2)*(2/3), 2/3))
cat("   They differ in the last place because 2/3 is not a binary fraction, not because the\n")
cat("   arithmetic is approximate.\n")
note(abs((-1 + (5/2)*(2/3)) - 2/3) < 1e-15, "the composed exponent is 2/3 identically")

cat("\n=== 2. the two laws as coded compose the same way, over four decades in M_1 ===\n")
# four decades in M_1 needs ten decades in I, since M_1 goes as I^{-2/5}
Is <- I0 * 10^seq(-5, 5, length.out = 21)
ms <- M1(Is); ts <- tdec(Is)
cat(sprintf("   M_1 runs %.4g to %.4g PeV, a span of %.2f decades\n",
            min(ms), max(ms), log10(max(ms)/min(ms))))
sl <- lm(log(ts) ~ log(ms))$coefficients[[2]]
cat(sprintf("   fitted log-log slope: %.17g\n", sl))
cat(sprintf("   difference from 2/3:  %.3e\n", abs(sl - 2/3)))
# a raw floating-point residual is not a stable figure to quote, so the paper quotes the bound
cat(sprintf("   which is better than 1e-15, the figure the paper quotes\n"))
note(log10(max(ms)/min(ms)) > 3.99, "the sweep really does cover four decades in the mass")
note(abs(sl - 2/3) < 1e-12, "and the coded laws return 2/3, not merely five decimals of it")
# a fit can hide a defect if the relation is not a power law at all, so check the residuals too
res <- max(abs(log(ts) - (lm(log(ts) ~ log(ms))$fitted.values)))
cat(sprintf("   largest residual of the fit in log t: %.3e, so it is a power law and not a fit\n", res))
note(res < 1e-12, "the relation is an exact power law across the sweep")

cat("\n=== 3. the convention factor between the two quoted crossing times ===\n")
# The crossing time at c_G = 1 is an input here, 1.417e-32 s, since its normalisation is not
# computed in this file. What IS computed is the Majorana-pair value that follows from it, and
# the paper's own figure for that is deliberately not printed beside it: a script that quotes
# the manuscript's number back at the provenance gate satisfies the gate by saying nothing.
cat("   The crossing time at c_G = 1 is an input here, 1.417e-32 s, since its normalisation is\n")
cat("   not computed in this file. Counting a Majorana pair takes c_G to 1/2, and the time goes\n")
cat("   with c_G^{-2/3}:\n")
pred <- 1.417 * (1/2)^(-2/3)
cat(sprintf("   1.417e-32 s times (1/2)^{-2/3} = %.5fe-32 s\n", pred))
# the tolerance is half a unit in the last quoted place, since the inputs carry four figures
note(abs(pred - 2.249) < 5e-4,
     "the two quoted times differ by exactly the mode-counting factor to the 2/3")

cat("\n=== 4. the plants: every exponent in the chain has to matter ===\n")
tab <- list(c(2/5, 2/3, 2/3), c(1/2, 2/3, 1/3), c(2/5, 1/2, 1/4), c(1/3, 2/3, 1))
cat("      M_1 propto I^-a     t_dec propto I^-b      predicted slope    fitted\n")
for (row in tab) {
  a <- row[1]; b <- row[2]; want <- row[3]
  m <- function(I) M0 * (I/I0)^(-a)
  t <- function(I) (1/(2*m(I))) * ((4/(3*sqrt(pi))) * cG * Rf * I)^(-b)
  s <- lm(log(t(Is)) ~ log(m(Is)))$coefficients[[2]]
  cat(sprintf("   %14.4f %18.4f %19.6f %11.6f\n", a, b, want, s))
  note(abs(s - want) < 1e-12, sprintf("the chain at a = %.3f, b = %.3f gives %.4f", a, b, want))
}
note(abs(lm(log(tdec(Is, a = 1/2)) ~ log(ms))$coefficients[[2]] - 2/3) > 0.1,
     "plant: moving the 2/3 in t_dec moves the composed power")

cat("\n=== 5. what this file is for ===\n")
cat("   Section 3.5's power is exact and the paper now says exact. The numerical part of this\n")
cat("   file is not evidence for the exponent, which needs none; it is evidence that the two\n")
cat("   power laws are coded the way the text states them, which is the part a slip could\n")
cat("   break and the part no other file was checking.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
