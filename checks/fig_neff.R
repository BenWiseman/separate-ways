#!/usr/bin/env Rscript
# fig_neff.R -- Section 8: the arithmetic that looks decisive and is not.
#
# Write the radiation density as rho_rad = rho_gamma (1 + 0.2271 N_eff), with
# 0.2271 = (7/8)(4/11)^(4/3). The Standard Model's 1.691 rho_gamma, duplicated by a partner
# sheet contributing at the same events, is an addition of 7.45 in N_eff units. Planck
# measures 2.99 +/- 0.17. The copy reading sits a third of a sigma away; the place reading
# is excluded at forty-four.
#
# The figure exists to show how wide that gap is, and the caption to say why the paper does
# NOT claim the exclusion: the fold puts the two sheets at different events, so the doubling
# never happens and both readings predict 3.044. Arithmetic: checks/calc/neff_silence_test.R

conv  <- (7/8) * (4/11)^(4/3)
sm    <- 3.044
rho_sm<- 1 + conv * sm
place <- (2 * rho_sm - 1) / conv
meas  <- 2.99; sig <- 0.17
stopifnot(abs(conv - 0.2271) < 5e-5, abs(rho_sm - 1.691) < 5e-4,
          abs(place - 10.49) < 5e-3)
cat(sprintf("  conversion %.4f, SM density %.4f rho_gamma, place reading %.2f\n", conv, rho_sm, place))
cat(sprintf("  copy  %.3f -> %.2f sigma from Planck\n", sm, abs(sm - meas)/sig))
cat(sprintf("  place %.2f  -> %.1f sigma from Planck\n", place, abs(place - meas)/sig))

ink <- "grey15"; ok <- "#1f4e79"; no <- "#a8400f"; band <- "#cfd8e3"
draw <- function() {
  par(mar = c(3.9, 1.2, 2.6, 1.2), mgp = c(2.5, 0.7, 0), xpd = FALSE)
  plot(NA, xlim = c(2.35, 11.3), ylim = c(0, 1), axes = FALSE,
       xlab = expression(paste("effective relativistic count  ", italic(N)[eff])), ylab = "")
  rect(2.35, 0, 11.3, 1, col = "#fbfbfa", border = NA)
  axis(1, at = c(3, 4, 5, 6, 7, 8, 9, 10, 11))

  # Planck, and the same band at five and ten sigma so the reader can count the gap.
  for (k in c(10, 5)) rect(meas - k*sig, 0.30, meas + k*sig, 0.70,
                           col = if (k == 10) "#eef2f6" else "#e2e9f0", border = NA)
  rect(meas - sig, 0.30, meas + sig, 0.70, col = band, border = NA)
  segments(meas, 0.30, meas, 0.70, col = "#5b6b7d", lwd = 1.5)
  text(meas, 0.755, expression(paste("Planck: 2.99 ", {} %+-% {}, " 0.17")), col = "#41505f",
       cex = 0.84, adj = 0.5)


  segments(2.35, 0.50, 11.3, 0.50, col = "grey75", lwd = 1)

  points(sm, 0.50, pch = 19, cex = 1.7, col = ok)
  text(sm + 0.30, 0.300, "a copy of our degrees of freedom:", col = ok, cex = 0.84, adj = 0)
  text(sm + 0.30, 0.220, "3.044, a third of a sigma out", col = ok, cex = 0.84, adj = 0)
  text(meas + 10*sig + 0.12, 0.115, "shaded bands: one, five and ten standard deviations",
       col = "#8a97a4", cex = 0.78, adj = 0)

  points(place, 0.50, pch = 19, cex = 1.7, col = no)
  text(place - 0.30, 0.685, "a place carrying its own, contributing at the same events:",
       col = no, cex = 0.84, adj = 1)
  text(place - 0.30, 0.605, "10.49, excluded at forty-four sigma", col = no, cex = 0.84, adj = 1)

  arrows(sm + 0.25, 0.88, place - 0.25, 0.88, code = 3, length = 0.07, col = ink, lwd = 1.2)
  text((sm + place)/2, 0.945, "an addition of 7.45, which is the whole Standard Model counted twice",
       col = ink, cex = 0.84, adj = 0.5)
  mtext("the test that looks decisive, and the reason it does not apply", side = 3,
        line = 1.25, cex = 0.88, col = ink)
}
for (f in c("papers/2_over_the_horizon/fig_companion_neff.pdf", "papers/2_over_the_horizon/fig_companion_neff.png")) {
  if (grepl("pdf$", f)) pdf(f, width = 7.8, height = 3.5) else png(f, width = 1210, height = 545, res = 150)
  draw(); invisible(dev.off()); cat("  wrote", f, "\n")
}
