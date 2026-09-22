# FIGURE: "The dark budget is closed, so a competitor subtracts."
#
# Section 3.4's argument is one a reader keeps getting backwards: a primordial black hole
# component does not sit ALONGSIDE the relic, it eats the budget the abundance assigns to it
# and drags the ceiling down. Draw the subtraction, then show that everything observation
# currently allows sits inside the quoted width -- the actual conclusion, which is a
# robustness result and not an exclusion.
#
# Labels are NUMBERED into a key rather than placed next to their markers. Four leader lines
# in a ten-decade panel cannot be laid out by estimating text width in data units; that has
# collided twice before. A key cannot collide.
#
# Base R only. Authored at 6.5 in, the width the manuscript places it at.

INK <- "#1b1b1b"; REL <- "#0b6e9e"; PBH <- "#c2461f"; GREY <- "#7e7e7e"
BAND <- rgb(0.11, 0.43, 0.62, 0.18)

M0 <- 491.6; SIG <- 2.0
ceil <- function(f) M0 * (1 - f)^(2/5)
f_for <- function(dM) uniroot(function(f) M0 - ceil(f) - dM, c(1e-12, 0.99))$root
F_WIDTH <- f_for(SIG)

# --- validate against the numbers 3.4 quotes, so the figure cannot drift from the text
stopifnot(abs(ceil(0.1) - 471.3) < 0.1,
          abs(ceil(0.5) - 372.6) < 0.1,
          abs(ceil(0.9) - 195.7) < 0.1,
          abs(F_WIDTH - 0.0101) < 5e-4)

# Green & Kavanagh 2007.10722 in the seed window, via pbh_ceiling.R; seeds via seed_margins.R
MARK <- list(
  list(f = 3.0e-10, col = REL,  key = "the red dots' seeds, as observed"),
  list(f = 3.0e-9,  col = GREY, key = "CMB accretion limit (tightest)"),
  list(f = 1.0e-4,  col = GREY, key = "dwarf heating limit (weakest)"),
  list(f = F_WIDTH, col = PBH,  key = "ceiling moves by one width"))

draw <- function() {
  layout(matrix(c(1, 2), 1, 2), widths = c(1, 1.46))

  # ---------------- (a) the budget itself ----------------
  par(mar = c(0.5, 0.5, 2.1, 0.5), family = "sans", xpd = NA)
  plot(NA, xlim = c(0, 1.42), ylim = c(-0.30, 2.62), axes = FALSE, xlab = "", ylab = "")
  mtext("(a)  one budget, not two", side = 3, line = 0.5, adj = 0, cex = 0.98, font = 2, col = INK)

  ys <- c(1.84, 0.98, 0.12); fs <- c(0, 0.1, 0.5); BW <- 0.86; H <- 0.36
  for (i in 1:3) {
    f <- fs[i]; y <- ys[i]
    rect(0, y, BW * (1 - f), y + H, col = REL, border = NA)
    if (f > 0) rect(BW * (1 - f), y, BW, y + H, col = PBH, border = NA)
    rect(0, y, BW, y + H, col = NA, border = INK, lwd = 1.1)
    text(0, y + H + 0.14,
         if (f == 0) "no black holes" else sprintf("%g%% black holes", 100 * f),
         col = INK, cex = 0.90, adj = 0)
    text(BW + 0.06, y + H/2, sprintf("%.1f", ceil(f)), col = INK, cex = 1.00, adj = 0)
  }
  text(BW + 0.06, ys[1] + H + 0.14, "ceiling, PeV", col = INK, cex = 0.86, adj = 0, font = 2)
  text(0, ys[3] - 0.26, "The relic keeps only what is left over.", col = REL, cex = 0.90, adj = 0)

  # ---------------- (b) what observation allows ----------------
  par(mar = c(3.6, 4.2, 2.1, 1.1))
  xs <- 10^seq(-10, log10(0.90), length.out = 900)
  plot(NA, xlim = c(-10, 0), ylim = c(175, 520), axes = FALSE, xlab = "", ylab = "")
  mtext("(b)  and all of it sits inside the width",
        side = 3, line = 0.5, adj = 0, cex = 0.98, font = 2, col = INK)

  rect(-10, M0 - SIG, 0, M0 + SIG, col = BAND, border = NA)
  lines(log10(xs), ceil(xs), col = INK, lwd = 2.6)

  axis(1, at = seq(-10, 0, 2), labels = parse(text = sprintf("10^%d", seq(-10, 0, 2))),
       cex.axis = 0.92, col = GREY, col.axis = INK, tck = -0.02, mgp = c(2, 0.6, 0))
  axis(2, at = seq(200, 500, 100), cex.axis = 0.92, las = 1,
       col = GREY, col.axis = INK, tck = -0.02, mgp = c(2, 0.7, 0))
  mtext("fraction of the dark matter in black holes", side = 1, line = 2.3, cex = 0.94, col = INK)
  mtext("ceiling on the relic mass  (PeV)", side = 2, line = 2.8, cex = 0.94, col = INK)

  # numbered markers ON the curve; nothing else goes near the curve
  for (i in seq_along(MARK)) {
    m <- MARK[[i]]; lx <- log10(m$f)
    points(lx, ceil(m$f), pch = 21, bg = "white", col = m$col, lwd = 2.0, cex = 1.55)
    text(lx, ceil(m$f), i, col = m$col, cex = 0.76, font = 2)
  }
  text(-6.0, M0 + 14, "the quoted width", col = REL, cex = 0.88, adj = 0.5)

  # key, in the empty lower-left the flat curve leaves free
  for (i in seq_along(MARK)) {
    yk <- 342 - (i - 1) * 33
    points(-9.72, yk, pch = 21, bg = "white", col = MARK[[i]]$col, lwd = 2.0, cex = 1.55)
    text(-9.72, yk, i, col = MARK[[i]]$col, cex = 0.76, font = 2)
    text(-9.35, yk, MARK[[i]]$key, col = INK, cex = 0.88, adj = 0)
  }
}

dir.create("pub/paper2/figs", showWarnings = FALSE, recursive = TRUE)
pdf("pub/paper2/figs/fig_budget.pdf", width = 6.5, height = 4.2); draw(); dev.off()
png("pub/paper2/figs/fig_budget.png", width = 1950, height = 1260, res = 300); draw(); dev.off()
cat(sprintf("wrote fig_budget.{pdf,png}   one width at f = %.5f\n", F_WIDTH))
