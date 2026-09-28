#!/usr/bin/env Rscript
# fig_contact_universal.R -- Section 5: the contact region, and which way of saying it survives.
#
# Inside the horizon the causality bound is d(angle) <= dr / sqrt(G) with G = -r^2 f. For
# Reissner-Nordstrom r^2 f = (r - r_+)(r - r_-), so G = (r_+ - r)(r - r_-) is a monic quadratic
# whose roots are the ends of the band, and one leg from r out to the horizon integrates to
#     2 arcsin sqrt( (r_+ - r) / (r_+ - r_-) ) .
# Write x for that depth fraction, 0 at the outer horizon and 1 at the inner. Two legs give
#     turn(x) = 4 arcsin sqrt(x)
# with the charge gone entirely: ONE curve, not a family. The antipodal map bills pi, and
#     4 arcsin sqrt(x) >= pi  <=>  sqrt(x) >= sin(pi/4)  <=>  x >= 1/2 ,
# so the budget crosses the bill at the midpoint of the band at every charge, exactly.
# That is the reaching statement. Confining contact to the inner half needs the maximiser
# argument, which holds for an uncharged hole and is open at charge; see A.15 and 5.1.
#
# Right panel: the same result stated two ways, only one of which survives charge. As a
# fraction of the band the contact region is one half and does not move. As a fraction of the
# horizon AREA it is a quarter at Q = 0 and runs to better than three quarters near extremality,
# so the quarter-area reading is a coincidence of r_- = 0 and not a thermodynamic statement.
# Arithmetic and planted failures: checks/calc/contact_charged.R

turn <- function(x) 4 * asin(sqrt(x))
stopifnot(abs(turn(0.5) - pi) < 1e-12, turn(0.4999) < pi, turn(0.5001) > pi)
s    <- function(Q) sqrt(1 - Q^2)
area <- function(Q) 1 / (1 + s(Q))^2           # (r_c/r_+)^2 with r_c = M
cat("=== the contact condition, normalised to the band ===\n\n")
cat("  turn(x) = 4 arcsin sqrt(x) reaches pi at x = 0.5 exactly, with no charge in it.\n\n")
cat("     Q/M     r_-/r_+    contact as fraction of band    as fraction of horizon area\n")
for (Q in c(0, 0.3, 0.6, 0.9, 0.99)) {
  rp <- 1 + s(Q); rm <- 1 - s(Q)
  cat(sprintf("   %6.3f   %8.5f   %22.7f   %24.5f\n", Q, rm/rp, 0.5, area(Q)))
}
stopifnot(abs(area(0) - 0.25) < 1e-12, area(0.99) > 0.76)
cat("\n  The band fraction does not move. The area fraction runs from 0.250 to 0.768, so the\n")
cat("  quarter-of-the-horizon reading is numerology and charge is what exposes it.\n")

ink <- "grey15"; c1 <- "#1f4e79"; mark <- "#a8400f"; grey <- "#8a97a4"
draw <- function() {
  par(mfrow = c(1, 2), mar = c(3.9, 4.6, 2.2, 1.0), mgp = c(2.5, 0.7, 0), xpd = FALSE)

  x <- seq(0, 1, length.out = 600)
  plot(NA, xlim = c(0, 1), ylim = c(0, 2*pi*1.02), axes = FALSE,
       xlab = "depth across the interior band",
       ylab = "angular path available, two legs")
  rect(0, 0, 1, 2*pi*1.02, col = "#fbfbfa", border = NA)
  axis(1, at = seq(0, 1, 0.25))
  axis(2, at = c(0, pi/2, pi, 3*pi/2, 2*pi),
       labels = c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)), las = 1)
  rect(0.5, 0, 1, 2*pi*1.02, col = "#eef1f5", border = NA)
  segments(0, pi, 1, pi, col = mark, lwd = 1.8)
  # crossed by the dotted midpoint guide at x = 0.5, and it straddles the edge of the shaded
  # right half, so a single-colour mask would show as a box. Accepted instead, in
  # checks/ink_overlaps_accepted.tsv.
  text(0.97, pi - 0.42, "what the antipodal map bills", col = mark, cex = 0.82, adj = 1)
  lines(x, turn(x), col = c1, lwd = 2.8)
  points(0.5, pi, pch = 19, cex = 1.3, col = mark)
  segments(0.5, 0, 0.5, pi, col = mark, lty = 3)
  text(0.53, 0.55, "the sheets touch", col = ink, cex = 0.82, adj = 0)
  text(0.53, 0.28, "in this half", col = ink, cex = 0.82, adj = 0)
  text(0.47, 5.30, "one curve, not a family:", col = c1, cex = 0.82, adj = 1)
  text(0.47, 4.95, expression(paste("turn = 4 arcsin ", sqrt(italic(x)), ", with")), col = c1, cex = 0.82, adj = 1)
  text(0.47, 4.60, "the charge gone entirely", col = c1, cex = 0.82, adj = 1)
  mtext("the budget crosses the bill at the midpoint, at every charge", side = 3, line = 0.7,
        cex = 0.86, col = ink)

  Qs <- seq(0, 0.995, length.out = 500)
  plot(NA, xlim = c(0, 1), ylim = c(0.14, 0.88), axes = FALSE,
       xlab = expression(paste("charge  ", italic(Q), "/", italic(M))),
       ylab = "size of the contact region")
  rect(0, 0.14, 1, 0.88, col = "#fbfbfa", border = NA)
  axis(1, at = seq(0, 1, 0.25)); axis(2, at = seq(0.2, 0.8, 0.2), las = 1)
  lines(Qs, rep(0.5, length(Qs)), col = c1, lwd = 2.8)
  lines(Qs, area(Qs), col = grey, lwd = 2.6, lty = 2)
  text(0.03, 0.615, "as a fraction of the interior band:", col = c1, cex = 0.82, adj = 0)
  text(0.03, 0.572, "one half, flat at every charge", col = c1, cex = 0.82, adj = 0)
  text(0.03, 0.215, "as a fraction of the horizon area: a quarter only", col = grey, cex = 0.82, adj = 0)
  text(0.03, 0.172, "at zero charge, and nothing special beyond it", col = grey, cex = 0.82, adj = 0)
  points(0, 0.25, pch = 19, cex = 1.1, col = grey); points(0, 0.5, pch = 19, cex = 1.1, col = c1)
  mtext("and the quarter-area coincidence does not survive it", side = 3, line = 0.7,
        cex = 0.86, col = ink)
}
for (f in c("paper/fig_companion_contact.pdf", "paper/fig_companion_contact.png")) {
  if (grepl("pdf$", f)) pdf(f, width = 8.6, height = 4.2) else png(f, width = 1330, height = 650, res = 150)
  draw(); invisible(dev.off()); cat("  wrote", f, "\n")
}
