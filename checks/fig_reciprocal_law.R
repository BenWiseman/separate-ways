#!/usr/bin/env Rscript
# fig_reciprocal_law.R -- Section 3's result made visible.
#
# A.10 fixes the two kernels of the conformal scalar, with r the radial variable of the
# static patch, dt the boost separation and gamma the angle between the two points:
#     Z_alpha = -(1-r^2) cosh(dt) - r^2 cos gamma
#     Z_J     = -(1-r^2) cosh(dt) + r^2 cos gamma
# and W ~ 1/(1-Z) for the conformal scalar, so the ratio the reciprocal law is about is
#     R(gamma) = G_alpha/G_J = (1 - Z_J)/(1 - Z_alpha) = (c - b cos gamma)/(c + b cos gamma)
# with c = 1 + (1-r^2) cosh(dt) > 0 and b = r^2 > 0.
#
# The whole of Section 3 is then one substitution. Sending gamma -> pi - gamma sends
# cos gamma -> -cos gamma, which swaps numerator and denominator:
#     log R(pi - gamma) = -log R(gamma),
# exactly, for every r and every boost. log R is an ODD function about ninety degrees.
# Expand it in Legendre polynomials of cos gamma: P_l(-x) = (-1)^l P_l(x), so every
# EVEN multipole integrates against an odd integrand and vanishes identically. Not to
# some order, and not for some state: identically.
#
# The right panel is the point of the figure. The even coefficients are not small, they
# are zero, and what the plot shows at 1e-17 is the arithmetic floor of the quadrature.

logR <- function(gam, r, dt) {
  cc <- 1 + (1 - r^2) * cosh(dt); b <- r^2
  log((cc - b * cos(gam)) / (cc + b * cos(gam)))
}

# Gauss-Legendre nodes and weights, Newton from the Chebyshev guess. Written out rather
# than taken from a package so the figure has no dependency beyond base R.
gauss_legendre <- function(n) {
  x <- cos(pi * (seq_len(n) - 0.25) / (n + 0.5))
  for (it in 1:100) {
    p0 <- rep(1, n); p1 <- x
    for (k in 2:n) { p2 <- ((2*k-1) * x * p1 - (k-1) * p0) / k; p0 <- p1; p1 <- p2 }
    dp <- n * (x * p1 - p0) / (x^2 - 1)
    dx <- -p1 / dp; x <- x + dx
    if (max(abs(dx)) < 1e-15) break
  }
  p0 <- rep(1, n); p1 <- x
  for (k in 2:n) { p2 <- ((2*k-1) * x * p1 - (k-1) * p0) / k; p0 <- p1; p1 <- p2 }
  dp <- n * (x * p1 - p0) / (x^2 - 1)
  list(x = x, w = 2 / ((1 - x^2) * dp^2))
}
Pl <- function(l, x) { if (l == 0) return(rep(1, length(x))); if (l == 1) return(x)
  p0 <- rep(1, length(x)); p1 <- x
  for (k in 2:l) { p2 <- ((2*k-1) * x * p1 - (k-1) * p0) / k; p0 <- p1; p1 <- p2 }
  p1 }

LMAX <- 12
gl <- gauss_legendre(400)
coefs <- function(r, dt) {
  f <- logR(acos(gl$x), r, dt)
  sapply(0:LMAX, function(l) (2*l + 1)/2 * sum(gl$w * f * Pl(l, gl$x)))
}

cases <- list(list(r = 0.60, dt = 0.0, col = "#1f4e79", lab = "r = 0.6, no boost"),
              list(r = 0.90, dt = 1.2, col = "#2e7d5b", lab = "r = 0.9, boosted"))

cat("=== the reciprocal law, and the multipoles it kills ===\n\n")
for (cs in cases) {
  g  <- seq(1e-6, pi - 1e-6, length.out = 4001)
  asym <- max(abs(logR(pi - g, cs$r, cs$dt) + logR(g, cs$r, cs$dt)))
  a  <- coefs(cs$r, cs$dt)
  ev <- max(abs(a[seq(1, LMAX + 1, by = 2)]))     # l = 0, 2, 4, ...
  od <- max(abs(a[seq(2, LMAX + 1, by = 2)]))     # l = 1, 3, 5, ...
  cat(sprintf("  %-22s  |log R(pi-g) + log R(g)| <= %.1e\n", cs$lab, asym))
  cat(sprintf("  %-22s  R(pi/2) - 1 = %.1e\n", "", exp(logR(pi/2, cs$r, cs$dt)) - 1))
  cat(sprintf("  %-22s  largest EVEN multipole %.2e, largest ODD %.4f\n\n", "", ev, od))
  stopifnot(asym < 1e-14, ev < 1e-13, od > 0.1)
  cs$a <- a
}

# The check must be able to fail: a kernel that is not antipodally reflected has even
# multipoles of the same size as its odd ones.
bad <- function(gam) log(1.6 - 0.4*cos(gam) + 0.3*cos(gam)^2)
f <- bad(acos(gl$x))
abad <- sapply(0:LMAX, function(l) (2*l+1)/2 * sum(gl$w * f * Pl(l, gl$x)))
cat(sprintf("  control, a kernel with no antipodal symmetry: largest even %.3f, odd %.3f\n",
            max(abs(abad[seq(1, LMAX+1, 2)])), max(abs(abad[seq(2, LMAX+1, 2)]))))
stopifnot(max(abs(abad[seq(1, LMAX+1, 2)])) > 0.05)
cat("  so the floor on the left is the law, not the quadrature.\n")

ink <- "grey15"; mark <- "#a8400f"
draw <- function() {
  par(mfrow = c(1, 2), mar = c(3.9, 4.0, 2.0, 0.8), mgp = c(2.4, 0.7, 0), xpd = NA)

  g <- seq(0.5, 179.5, by = 0.25)
  plot(NA, xlim = c(0, 180), ylim = c(-1.45, 1.45), axes = FALSE,
       xlab = expression(paste("angle ", gamma, " between the two points, degrees")),
       ylab = expression(paste("log ", italic(R), "(", gamma, ")")))
  rect(0, -1.45, 180, 1.45, col = "#fbfbfa", border = NA)
  axis(1, at = seq(0, 180, 45)); axis(2, at = seq(-1, 1, 0.5), las = 1)
  # segments, not abline: xpd = NA is needed for the axis titles and would let an
  # abline run clean across the device and into the other panel.
  segments(0, 0, 180, 0, col = "grey70", lwd = 1)
  segments(90, -1.45, 90, 1.45, col = "grey70", lty = 2)
  for (cs in cases) lines(g, logR(g * pi/180, cs$r, cs$dt), col = cs$col, lwd = 2.4)
  points(90, 0, pch = 19, cex = 1.2, col = mark)
  text(5, 1.20, expression(paste(italic(R), " = 1 exactly at ninety degrees,")), col = mark, cex = 0.82, adj = 0)
  # the dashed ninety-degree guide ran through this; masked, the guide being a thin one
  source("checks/fig_label.R")
  lab_on(5, 1.01, "at every radius and every boost", bg = "#fbfbfa",
         col = mark, cex = 0.82, adj = 0)
  text(174, -0.78, cases[[1]]$lab, col = cases[[1]]$col, cex = 0.82, adj = 1)
  text(174, -0.99, cases[[2]]$lab, col = cases[[2]]$col, cex = 0.82, adj = 1)
  mtext("odd about ninety degrees", side = 3, line = 0.4, cex = 0.86, col = ink)

  a1 <- coefs(cases[[1]]$r, cases[[1]]$dt); a2 <- coefs(cases[[2]]$r, cases[[2]]$dt)
  flo <- 1e-18
  h <- function(v) pmax(abs(v), flo)
  plot(NA, xlim = c(-0.5, LMAX + 0.5), ylim = c(log10(flo), 1.1), axes = FALSE,
       xlab = expression(paste("multipole ", italic(l))),
       ylab = expression(paste("|", italic(a[l]), "|")))
  rect(-0.5, log10(flo), LMAX + 0.5, 1.1, col = "#fbfbfa", border = NA)
  axis(1, at = seq(0, LMAX, 2))
  axis(2, at = seq(-18, 0, 3), labels = parse(text = sprintf("10^%d", seq(-18, 0, 3))), las = 1)
  for (l in 0:LMAX) {
    for (j in 1:2) {
      v <- if (j == 1) a1[l+1] else a2[l+1]
      xo <- if (j == 1) -0.16 else 0.16
      rect(l + xo - 0.14, log10(flo), l + xo + 0.14, log10(h(v)),
           col = cases[[j]]$col, border = NA)
    }
  }
  segments(-0.5, log10(2e-16), LMAX + 0.5, log10(2e-16), col = mark, lty = 3, lwd = 1.2)
  # The floor label sits at a height every odd bar passes through, so it was printed over by
  # the bars at l = 5, 7 and 9 in the built PDF. Gate 19 measures labels against the panel and
  # not against what else is drawn, so it passed. A panel-coloured patch behind the text is the
  # fix, and it is drawn from the string's own measured extent rather than a guessed box.
  lab <- "double-precision floor"; lx <- 3.6; ly <- log10(2e-16) + 1.1
  rect(lx - 0.08, ly - 0.55*strheight(lab, cex = 0.80),
       lx + strwidth(lab, cex = 0.80) + 0.08, ly + 0.70*strheight(lab, cex = 0.80),
       col = "#fbfbfa", border = NA)
  text(lx, ly, lab, col = mark, cex = 0.80, adj = 0)
  text(0, -16.4, "even", col = ink, cex = 0.82, adj = 0.5)
  text(1, 0.55, "odd", col = ink, cex = 0.82, adj = 0.5)

  mtext("so every even multipole is zero", side = 3, line = 0.4, cex = 0.86, col = ink)
}
for (f in c("papers/2_over_the_horizon/fig_companion_reciprocal.pdf", "papers/2_over_the_horizon/fig_companion_reciprocal.png")) {
  if (grepl("pdf$", f)) pdf(f, width = 8.4, height = 3.9) else png(f, width = 1300, height = 620, res = 150)
  draw(); invisible(dev.off()); cat("  wrote", f, "\n")
}
