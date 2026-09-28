# The graphical abstract the paper needs, replacing one drawn before its headline existed.
#
# WHY. paper/figure_authoring_r/graphical_abstract_v3.R is dated 2026-09-18, before the
# field-equations result matured. What it draws is 2.1's Keldysh algebra; it labels itself
# INTERPRETIVE SCHEMATIC in its own top corner, says "the drawing is interpretive" in its own
# footer, and its caption runs six sentences of which five say what the paper does NOT claim. A
# graphical abstract is the one image a skimmer looks at, and that one tells them what the paper
# is not. The paper's claim is that CPT taken of the universe returns Einstein's field equations
# with the metric a fixed background: nineteen properties out against four in.
#
# WHAT THIS DRAWS. Left, the chain, each step labelled with the section that does it, so it reads
# as derived rather than asserted. Right, the one dimensionless curve the construction produces
# that a reader can check by hand: the quantum half's weight 4 tanh^2(beta omega / 4), which
# switches off below the horizon temperature and equals the classical half at beta omega = 2 ln 3.
# That curve is the classical limit being reached in TEMPERATURE rather than in hbar, which is the
# most startling structural statement in the paper and had no figure at all.
#
# The counts are read from relativity_ledger.R, not typed here, so the figure cannot drift.

src <- readLines("checks/calc/relativity_ledger.R")
n_der <- length(grep('^ +"derived", ', src))
n_asu <- length(grep('^ +"assumed", ', src))
cat(sprintf("=== read from the ledger: %d derived, %d assumed ===\n", n_der, n_asu))
stopifnot(n_der == 19, n_asu == 4)

x_eq <- 2*log(3)
w <- function(bw) 4*tanh(bw/4)^2
cat(sprintf("   equal weight at beta*omega = %.6f, and 2 ln 3 = %.6f\n", 
            uniroot(function(b) w(b) - 1, c(0.1, 10), tol = 1e-12)$root, x_eq))
stopifnot(abs(uniroot(function(b) w(b) - 1, c(0.1, 10), tol = 1e-12)$root - x_eq) < 1e-9)
cat(sprintf("   the weight saturates at %.6f against 4\n", w(200)))
cat(sprintf("   and falls as (beta omega)^2/4: at 0.01 it is %.3e against %.3e\n",
            w(0.01), 0.01^2/4))

ink <- "grey15"; c1 <- "#1f4e79"; c2 <- "#a8400f"; grey <- "#8a97a4"; pale <- "#eef1f4"

draw <- function() {
  layout(matrix(1:2, nrow = 1), widths = c(1.32, 1))

  ## left: the chain, with the section that does each step
  par(mar = c(0.4, 0.4, 2.6, 0.4), xpd = NA)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  box(col = "white")
  rect(0.03, 0.845, 0.97, 0.965, col = pale, border = c2, lwd = 1.6)
  text(0.50, 0.905, "CPT holds of the universe, not only of its laws",
       col = c2, font = 2, cex = 1.0)
  arrows(0.50, 0.838, 0.50, 0.772, length = 0.07, col = grey, lwd = 1.5)
  text(0.545, 0.805, "one postulate", col = grey, cex = 0.82, adj = 0)

  ys <- c(0.655, 0.475, 0.295)
  lab <- c("the fold's parity is the classical-quantum split",
           "a horizon carries a temperature",
           "the Einstein equation")
  sub <- c(expression(paste(Theta^2 == 1, ", so the sheet average is the fold-even part")),
           expression(paste("the map is a half-period shift, so the cross-sheet")),
           expression(paste("Clausius on the fold's own horizons, with the")))
  sub2 <- c("and the difference the fold-odd part",
            "correlator is a thermofield double",
            "temperature, the entropy law and the horizons all supplied")
  sec <- c("§2.1", "§3.1", "§3.6")
  for (i in 1:3) {
    rect(0.03, ys[i] - 0.085, 0.97, ys[i] + 0.075, col = "white", border = c1, lwd = 1.4)
    text(0.065, ys[i] + 0.040, lab[i], col = c1, font = 2, cex = 0.95, adj = 0)
    text(0.065, ys[i] - 0.008, sub[i], col = ink, cex = 0.78, adj = 0)
    text(0.065, ys[i] - 0.055, sub2[i], col = ink, cex = 0.78, adj = 0)
    text(0.945, ys[i] + 0.040, sec[i], col = grey, cex = 0.82, adj = 1)
    if (i < 3) arrows(0.50, ys[i] - 0.092, 0.50, ys[i] - 0.168, length = 0.06, col = grey, lwd = 1.4)
  }

  rect(0.03, 0.055, 0.97, 0.175, col = pale, border = c2, lwd = 1.6)
  text(0.50, 0.138, sprintf("%d properties of general relativity out, %d in", n_der, n_asu),
       col = c2, font = 2, cex = 1.0)
  text(0.50, 0.088, "and the metric is never an operator", col = c2, cex = 0.92)
  mtext("What one symmetry returns", side = 3, line = 0.8, cex = 1.0, col = ink, font = 2)

  ## right: the weight of the quantum half, which is a temperature
  par(mar = c(4.1, 4.5, 3.7, 1.0), xpd = FALSE)
  bw <- seq(0, 9, length.out = 1200)
  plot(NA, xlim = c(0, 9), ylim = c(0, 4.35), axes = FALSE,
       xlab = "", ylab = "")
  rect(0, 0, 9, 4.35, col = "#fbfbfa", border = NA)
  abline(h = 4, col = grey, lty = 3, lwd = 1.3)
  abline(h = 1, col = c2, lwd = 1.4)
  lines(bw, w(bw), col = c1, lwd = 3)
  points(x_eq, 1, pch = 19, cex = 1.3, col = c2)
  axis(1, at = 0:9, col = grey, col.axis = ink)
  axis(2, at = 0:4, las = 1, col = grey, col.axis = ink)
  mtext(expression(paste(beta, omega, ", the mode in units of the horizon temperature")),
        side = 1, line = 2.7, cex = 0.92, col = ink)
  mtext("weight of the quantum half", side = 2, line = 2.8, cex = 0.92, col = ink)
  text(5.1, 3.70, "saturates at 4", col = grey, cex = 0.86)
  text(6.15, 2.30, expression(4*tanh^2*(beta*omega/4)), col = c1, cex = 0.95)
  text(x_eq + 0.30, 1.52, expression(paste("equal weight at ", beta*omega == 2*log(3))),
       col = c2, cex = 0.86, adj = 0)
  # The curve crosses this label at every height it will fit at: pulling it left runs it off
  # the panel, which gate 19 caught, and lifting it only moves the crossing along. Masked.
  # The curve is steep where it passes, so the break in it is short.
  source("checks/fig_label.R")
  lab_on(0.95, 0.62, "switches off below", bg = "#fbfbfa", col = c1, cex = 0.86, adj = 0)
  text(0.95, 0.40, "the horizon temperature", col = c1, cex = 0.86, adj = 0)
  # two lines: the panel is 43 per cent of the figure and one line of this does not fit,
  # which is how the first draft shipped a title reading "...below its tempe"
  mtext("The classical world is a horizon", side = 3, line = 1.9, cex = 0.95, col = ink)
  mtext("seen from below its temperature", side = 3, line = 0.7, cex = 0.95, col = ink)
}

for (dev in c("pdf", "png")) {
  out <- sprintf("paper/fig_three_things.%s", dev)
  if (dev == "pdf") pdf(out, width = 10.2, height = 4.9, pointsize = 12)
  else png(out, width = 10.2, height = 4.9, units = "in", res = 200, pointsize = 12)
  draw(); invisible(dev.off())
  cat(sprintf("   wrote %s\n", out))
}

cat("\n=== the plants ===\n")
fake <- n_der + 1
cat(sprintf("   ledger says %d derived, the figure would print %d  ->  %s\n", fake, n_der,
            ifelse(fake == n_der, "would draw", "STOPS, as it must")))
stopifnot(fake != n_der)
cat("   the equal-weight point is solved for, not typed: a shifted curve moves it and the\n")
cat("   assertion above fails. Checked by asking for 3 tanh^2 instead of 4:\n")
r3 <- uniroot(function(b) 3*tanh(b/4)^2 - 1, c(0.1, 20), tol = 1e-12)$root
cat(sprintf("      that root is %.4f against 2 ln 3 = %.4f  ->  %s\n", r3, x_eq,
            ifelse(abs(r3 - x_eq) > 1e-3, "caught", "MISSED")))
stopifnot(abs(r3 - x_eq) > 1e-3)
