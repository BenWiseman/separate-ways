# The two splittings of one state, and why they cross.
#
# Section 2 needs the two sheets to be two, and the tempting route is to read the no-boundary
# wavefunction's expanding and contracting WKB branches as the fold's two sheets, since one
# involution alpha exchanges the members of both pairs.  The route does not close, and the reason
# is a partition mismatch that a picture settles faster than a paragraph.
#
# Closed de Sitter in conformal coordinates: ds^2 = (H cos eta)^-2 (-d eta^2 + d Omega_3^2), with
# eta in (-pi/2, pi/2) and chi the polar angle from the observer's worldline.  Then
#
#   the PATCH splitting is by causal position:  the right static patch is chi <= pi/2 - |eta|,
#   the BRANCH splitting is by global time:     the expanding branch is eta > 0,
#
# and the patch diamond is bisected by the branch line.  Every one of the four cells is non-empty,
# so an operator localised in the patch has support on both branches and no element of the patch
# algebra determines the branch label.  The splittings cross; they do not coincide.
#
# The equality of the two halves is exact and not numerical luck: the patch condition
# chi <= pi/2 - |eta| and the conformal metric are both even in eta, so the reflection eta -> -eta
# is a measure-preserving bijection between the halves.  The integral below is a check of the
# arithmetic, not the source of the claim.

half_measure <- function(lo, hi, n = 200001) {
  # conformal 4-volume element on the Einstein cylinder: 4 pi sin^2(chi) d chi d eta
  eta <- seq(lo, hi, length.out = n)
  chimax <- pmax(pi / 2 - abs(eta), 0)
  inner <- 4 * pi * (chimax / 2 - sin(2 * chimax) / 4)   # int_0^chimax sin^2 chi d chi
  w <- c(1, rep(c(4, 2), length.out = n - 2), 1) * ((hi - lo) / (n - 1)) / 3
  sum(w * inner)
}

expanding  <- half_measure(0, pi / 2)
contracting <- half_measure(-pi / 2, 0)
whole      <- half_measure(-pi / 2, pi / 2)

cat("=== the patch straddles both branches in equal measure ===\n\n")
cat(sprintf("  conformal 4-volume, patch on the expanding branch   = %.10f\n", expanding))
cat(sprintf("  conformal 4-volume, patch on the contracting branch = %.10f\n", contracting))
cat(sprintf("  relative difference                                 = %.3e\n",
            abs(expanding - contracting) / whole))
cat(sprintf("  each half as a fraction of the whole patch          = %.10f\n", expanding / whole))
stopifnot(abs(expanding - contracting) / whole < 1e-12,
          abs(expanding / whole - 0.5) < 1e-12)

# A plant: shift the branch line off the throat and the halves must stop matching, or the
# integrator is not measuring what the panel draws.
off <- half_measure(0.3, pi / 2)
cat(sprintf("\n  planted: put the branch line at eta = 0.3 instead and the larger half takes\n"))
cat(sprintf("  %.4f of the patch, so the check is capable of failing.\n", (whole - off) / whole))
stopifnot(abs((whole - off) / whole - 0.5) > 0.1)

ink <- "grey15"; c1 <- "#1f4e79"; c2 <- "#a8400f"; grey <- "#8a97a4"

draw <- function() {
  par(mfrow = c(1, 2), mar = c(4.0, 4.4, 2.2, 1.0), mgp = c(2.5, 0.7, 0), xpd = FALSE)

  # ---- left: the conformal square, the patch diamond, the branch line ----
  plot(NA, xlim = c(0, pi), ylim = c(-pi / 2, pi / 2), axes = FALSE,
       xlab = expression(paste("polar angle  ", chi)),
       ylab = expression(paste("conformal time  ", eta)))
  rect(0, -pi / 2, pi, pi / 2, col = "#fbfbfa", border = NA)
  axis(1, at = c(0, pi / 2, pi), labels = c("0", expression(pi / 2), expression(pi)))
  axis(2, at = c(-pi / 2, 0, pi / 2),
       labels = c(expression(-pi / 2), "0", expression(pi / 2)), las = 1)

  # the two branch regions, shaded by global time
  rect(0, 0, pi, pi / 2, col = "#eef3f8", border = NA)
  rect(0, -pi / 2, pi, 0, col = "#f7efe9", border = NA)

  # right static patch: chi <= pi/2 - |eta|
  ee <- seq(-pi / 2, pi / 2, length.out = 400)
  polygon(c(rep(0, length(ee)), rev(pi / 2 - abs(ee))), c(ee, rev(ee)),
          col = adjustcolor(c1, 0.22), border = c1, lwd = 2.4)

  abline(h = 0, col = c2, lwd = 2.6)
  text(pi - 0.06, pi / 4, "expanding", col = c1, cex = 0.92, adj = 1)
  text(pi - 0.06, -pi / 4, "contracting", col = c2, cex = 0.92, adj = 1)
  text(0.30, 0.62, "right", col = c1, cex = 0.92, adj = 0)
  text(0.30, 0.44, "patch", col = c1, cex = 0.92, adj = 0)
  text(pi / 2 + 0.08, -0.16, expression(paste(chi, " = ", pi / 2, " - |", eta, "|")),
       col = c1, cex = 0.88, adj = 0)
  mtext("one state, two splittings, drawn together", side = 3, line = 0.6, cex = 0.92, col = ink)

  # ---- right: the patch's extent against global time, and the equal halves ----
  plot(NA, xlim = c(-pi / 2, pi / 2), ylim = c(0, 1.02 * pi / 2), axes = FALSE,
       xlab = expression(paste("conformal time  ", eta)),
       ylab = expression(paste("angular extent of the patch  ", chi[max])))
  rect(-pi / 2, 0, pi / 2, 1.02 * pi / 2, col = "#fbfbfa", border = NA)
  axis(1, at = c(-pi / 2, 0, pi / 2),
       labels = c(expression(-pi / 2), "0", expression(pi / 2)))
  axis(2, at = c(0, pi / 4, pi / 2),
       labels = c("0", expression(pi / 4), expression(pi / 2)), las = 1)

  ep <- seq(0, pi / 2, length.out = 300)
  en <- seq(-pi / 2, 0, length.out = 300)
  polygon(c(ep, rev(ep)), c(pi / 2 - abs(ep), rep(0, length(ep))),
          col = adjustcolor(c1, 0.30), border = NA)
  polygon(c(en, rev(en)), c(pi / 2 - abs(en), rep(0, length(en))),
          col = adjustcolor(c2, 0.30), border = NA)
  lines(c(en, ep), pi / 2 - abs(c(en, ep)), col = ink, lwd = 2.6)
  segments(0, 0, 0, pi / 2, col = c2, lwd = 2.6)

  text(-pi / 4, 0.34, sprintf("%.3f", contracting / whole), col = c2, cex = 0.98)
  text(-pi / 4, 0.19, "of the patch", col = c2, cex = 0.82)
  text(pi / 4, 0.34, sprintf("%.3f", expanding / whole), col = c1, cex = 0.98)
  text(pi / 4, 0.19, "of the patch", col = c1, cex = 0.82)
  mtext("the branch line bisects it, exactly", side = 3, line = 0.6, cex = 0.92, col = ink)
}

for (dev in c("pdf", "png")) {
  out <- sprintf("paper/fig_companion_splittings.%s", dev)
  if (dev == "pdf") pdf(out, width = 9.6, height = 4.3, pointsize = 12)
  else png(out, width = 9.6, height = 4.3, units = "in", res = 150, pointsize = 12)
  draw(); dev.off(); cat(sprintf("  wrote %s\n", out))
}
