#!/usr/bin/env Rscript
# ###########################################################################
# SUPERSEDED by neff_silence_test.R. Its section 3 treats the events needed as
# resolution/sqrt(N), which is right for random per-event resolution and wrong
# here: the limiting uncertainty on an EeV cascade line position is the ABSOLUTE
# ENERGY SCALE, a common-mode systematic that does not average down. Against a
# 7.18 per cent signal and a 15-25 per cent scale systematic the significance
# saturates below one sigma at any event count. The "160 events at three sigma"
# conclusion below is wrong. The g_*-scaling arithmetic in sections 0 to 2 is
# correct and is why this file is kept.
# ###########################################################################
# partner_sheet_test.R -- the companion's empirical anchor, priced.
#
# The fold's central ontological question is whether the partner sheet is a
# COPY (a relabelling of the same degrees of freedom) or a PLACE (a second
# region with its own stress-energy). The companion is the paper that decides
# what the second sheet IS, so this test belongs here even though the number is
# currently a one-paragraph caveat in section 4.2 of the cosmology paper.
#
# The discriminator: the expansion rate at the bang counts g_* degrees of
# freedom. A copy contributes nothing, so g_* = 106.75. A place doubles it.
# Since M_1 ~ g_*^(1/10), the two readings give different dark-matter ceilings
# and therefore different two-body neutrino line energies.

banner <- function() {
  cat("\n")
  cat("  ############################################################################\n")
cat("  ## SUPERSEDED by neff_silence_test.R. Do not quote the event count below.\n")
  cat("  ############################################################################\n")
cat("  Section 3 treats the events needed as resolution/sqrt(N), which is wrong here: the limiting\n")
cat("  uncertainty on an EeV cascade line is the absolute energy scale, a common-mode systematic that\n")
cat("  does not average down, so the significance saturates below one sigma at any event count. The\n")
cat("  '160 events at three sigma' conclusion below is wrong. Sections 0 to 2 are correct.\n")
  cat("  ############################################################################\n\n")
}
banner()
g_one <- 106.75
M_copy <- 491.6     # PeV, quoted ceiling, one sheet's content
w_M    <- 2.0       # PeV, quoted width on the ceiling
w_line <- 1.0       # PeV, quoted width on the line

cat("=== 0. validation: reproduce the quoted alternative from the scaling ===\n")
M_place <- M_copy * (2*g_one/g_one)^(1/10)
cat(sprintf("  M_1 if the sheet is a copy  : %.1f PeV  (quoted 491.6)\n", M_copy))
cat(sprintf("  M_1 if the sheet is a place : %.2f PeV  (quoted 526.9)\n", M_place))
stopifnot(abs(M_place - 526.9) < 0.1)
# planted failure: a wrong exponent must not reproduce it
M_bad <- M_copy * 2^(1/5)
cat(sprintf("  with the exponent wrong (1/5): %.1f PeV, which misses by %.1f PeV\n",
            M_bad, abs(M_bad - 526.9)))
stopifnot(abs(M_bad - 526.9) > 10)
cat("  scaling reproduces the quoted value and the check can say no\n")

cat("\n=== 1. the separation, in the units each observable is quoted in ===\n")
dM <- M_place - M_copy
line_copy  <- M_copy/2; line_place <- M_place/2
dline <- line_place - line_copy
cat(sprintf("  ceiling : %.1f vs %.1f PeV, gap %.1f PeV = %.1f quoted widths\n",
            M_copy, M_place, dM, dM/w_M))
cat(sprintf("  line    : %.1f vs %.1f PeV, gap %.2f PeV = %.1f quoted widths\n",
            line_copy, line_place, dline, dline/w_line))

cat("\n=== 2. what fractional accuracy on a measured line settles it? ===\n")
cat("  The two hypotheses are separated by a fixed FRACTION of the line energy,\n")
cat("  since both scale together. That fraction is what an experiment must beat.\n\n")
frac <- dline/line_copy
cat(sprintf("  fractional separation: %.4f  (%.2f per cent of the line energy)\n", frac, 100*frac))
cat(sprintf("  so a line measured to better than %.1f per cent distinguishes them at 1 sigma,\n",
            100*frac))
cat(sprintf("  and to better than %.2f per cent at 3 sigma.\n", 100*frac/3))

cat("\n=== 3. how does that compare with the resolution such an instrument has? ===\n")
for (res in c(0.50, 0.30, 0.15, 0.05)) {
  n_needed <- (res/frac)^2            # sigma_mean = res/sqrt(N) must beat frac
  cat(sprintf("  at %2.0f%% per-event energy resolution: %s events to separate at 1 sigma, %s at 3\n",
              100*res,
              format(ceiling(n_needed), big.mark=","),
              format(ceiling(9*n_needed), big.mark=",")))
}

cat("\n=== 4. flatly ===\n")
# interpolated, never typed beside the table: three drifts in one session already
res_ref <- 0.30
n3  <- ceiling(9*(res_ref/frac)^2)
n1  <- ceiling((res_ref/frac)^2)
n_endpoint <- 100                      # the cosmology paper's own population test
cat(sprintf("
  The two readings of what the partner sheet IS are %.1f line-widths apart, which
  sounds decisive and is the wrong way to quote it, because that width is a
  propagated theory error rather than an experimental one. The quantity that matters is
  the fractional separation, %.1f per cent of the line energy, since both
  hypotheses scale the line together. That is what an experiment has to beat.

  At the %.0f per cent per-event energy resolution the cosmology paper assumes in
  this band, telling the two apart needs about %d events at one sigma and %d at
  three. The paper's own endpoint test needs of order %d events. So this is the
  SAME measurement, the same instrument and the same order of exposure, asked a
  different question: %s

  That is the strongest empirical claim this paper has, and it is unusual in kind.
  The question is ontological, whether the second sheet is a copy of our degrees
  of freedom or a place with its own, and the answer is a neutrino energy. It is a
  claim about what a feasible measurement would decide rather than about present
  data, and it has to be written that way.
", dline/w_line, 100*frac, 100*res_ref, n1, n3, n_endpoint,
   if (n3 <= 3*n_endpoint)
     sprintf("%d against %d events is the same regime, not a harder experiment.", n3, n_endpoint)
   else
     sprintf("%d against %d events makes it a materially harder experiment.", n3, n_endpoint)))


banner()
