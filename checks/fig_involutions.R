#!/usr/bin/env Rscript
# fig_involutions.R -- A.15's classification, drawn. Three candidates, one survivor.
#
# Kruskal's metric depends on U and V only through the product UV, so a linear isometry is
# either a diagonal (U,V) -> (lam U, V/lam) or a swap (U,V) -> (lam V, U/lam). Squaring kills
# all but three up to a boost:
#     N  : (U,V) -> (-U,-V)      S : (U,V) -> (V,U)      SN : (U,V) -> (-V,-U)
# In Kruskal coordinates T = (U+V)/2 and X = (V-U)/2 these are, respectively, the point
# reflection through the origin, the mirror in the vertical axis and the mirror in the
# horizontal one. That is the whole classification, and it is why the fold's map is forced.
#
# What separates them:
#   S  preserves time orientation, so CPT excludes it. It is the RP^3 geon.
#   SN sends each exterior to ITSELF, so it offers no mirror sheet at all, and
#      whatever contact meant it would mean it where we can look.
#   N  reverses time orientation, carries our exterior to the other one, and is rigid: no
#      boost composes with it to give another involution.
# Table and numerics: checks/calc/fold_map_classification.R

maps <- list(
  list(k = "N",  f = function(T, X) c(-T, -X), fix = "point",
       ttl = expression(paste(italic(N), ":  (", italic(U), ",", italic(V), ") ", symbol("\256"),
                              " (-", italic(U), ", -", italic(V), ")")),
       v1 = "reverses time orientation and carries",
       v2 = "our exterior to the other one: the fold", good = TRUE),
  list(k = "S",  f = function(T, X) c(T, -X), fix = "vert",
       ttl = expression(paste(italic(S), ":  (", italic(U), ",", italic(V), ") ", symbol("\256"),
                              " (", italic(V), ", ", italic(U), ")")),
       v1 = "preserves time orientation, so CPT",
       v2 = "excludes it: this is the RP3 geon", good = FALSE),
  list(k = "SN", f = function(T, X) c(-T, X), fix = "horiz",
       ttl = expression(paste(italic(SN), ":  (", italic(U), ",", italic(V), ") ", symbol("\256"),
                              " (-", italic(V), ", -", italic(U), ")")),
       v1 = "sends each exterior to itself, so there",
       v2 = "is no mirror sheet for it to make", good = FALSE))

# The three are isometries and involutions; a fourth map is neither. Checked, not asserted.
UV <- function(T, X) (T - X) * (T + X) * -1   # U = T - X, V = T + X, so UV = T^2 - X^2
set.seed(1); Ts <- runif(500, -2, 2); Xs <- runif(500, -2, 2)
for (m in maps) {
  im <- mapply(function(t, x) m$f(t, x), Ts, Xs)
  stopifnot(max(abs(UV(im[1,], im[2,]) - UV(Ts, Xs))) < 1e-12)          # isometry
  im2 <- mapply(function(t, x) m$f(t, x), im[1,], im[2,])
  stopifnot(max(abs(im2[1,] - Ts)) < 1e-12, max(abs(im2[2,] - Xs)) < 1e-12)  # involution
}
bad <- function(T, X) c(-T, 2*X)
stopifnot(max(abs(UV(-Ts, 2*Xs) - UV(Ts, Xs))) > 1)                      # and one that is not
cat("=== all three are isometries and involutions of the (U,V) plane ===\n")
cat("   N, S and SN each preserve UV to 1e-12 over 500 random points and square to the\n")
cat("   identity; the control map (T,X) -> (-T,2X) fails the first test, so it can fail.\n\n")
p <- c(0.30, 1.20)
for (m in maps) {
  q <- m$f(p[1], p[2])
  reg <- function(v) if (abs(v[2]) > abs(v[1])) (if (v[2] > 0) "R (our exterior)" else "L (the mirror)") else
                     (if (v[1] > 0) "F (interior)" else "P (interior)")
  cat(sprintf("   %-3s sends %s to %s\n", m$k, reg(p), reg(q)))
}

ink <- "grey15"; hor <- "#9fb3c8"; sing <- "#8a8f95"; c1 <- "#1f4e79"; mark <- "#a8400f"
ok <- "#2e7d5b"
panel <- function(m) {
  par(mar = c(3.4, 1.0, 2.6, 1.0), xpd = FALSE)
  plot(NA, xlim = c(-2.05, 2.05), ylim = c(-2.05, 2.05), axes = FALSE, xlab = "", ylab = "", asp = 1)
  xs <- seq(-2.05, 2.05, length.out = 400)
  polygon(c(-2.05, xs, 2.05), c(2.05, sqrt(1 + xs^2), 2.05), col = "#f0eeec", border = NA)
  polygon(c(-2.05, xs, 2.05), c(-2.05, -sqrt(1 + xs^2), -2.05), col = "#f0eeec", border = NA)
  lines(xs, sqrt(1 + xs^2), col = sing, lwd = 1.8); lines(xs, -sqrt(1 + xs^2), col = sing, lwd = 1.8)
  segments(-2.05, -2.05, 2.05, 2.05, col = hor, lwd = 1.4)
  segments(-2.05, 2.05, 2.05, -2.05, col = hor, lwd = 1.4)
  # both sat on the horizontal axis line; lifted off it
  text(1.62, 0.22, "R", col = ink, cex = 0.92); text(-1.62, 0.22, "L", col = ink, cex = 0.92)
  text(0, 1.30, "F", col = ink, cex = 0.92);  text(0, -1.30, "P", col = ink, cex = 0.92)

  if (m$fix == "point") points(0, 0, pch = 4, cex = 1.3, col = ok, lwd = 2)
  if (m$fix == "vert")  segments(0, -1, 0, 1, col = ok, lwd = 2.4)
  if (m$fix == "horiz") segments(-2.05, 0, 2.05, 0, col = ok, lwd = 2.4)

  q <- m$f(p[1], p[2])
  arrows(p[2], p[1], q[2], q[1], length = 0.08, col = mark, lwd = 1.8)
  points(p[2], p[1], pch = 19, cex = 1.1, col = c1)
  points(q[2], q[1], pch = 19, cex = 1.1, col = mark)
  mtext(m$ttl, side = 3, line = 1.15, cex = 0.86, col = ink)
  mtext(m$v1, side = 1, line = 0.9, cex = 0.80, col = if (m$good) ok else mark)
  mtext(m$v2, side = 1, line = 1.85, cex = 0.80, col = if (m$good) ok else mark)
}
draw <- function() { par(mfrow = c(1, 3), oma = c(0.6, 0.4, 2.0, 0.4)); for (m in maps) panel(m)
  mtext("three involutions of the eternal hole, and only one can be the fold",
        outer = TRUE, side = 3, line = 0.4, cex = 0.92, col = ink) }
for (f in c("paper/fig_companion_involutions.pdf", "paper/fig_companion_involutions.png")) {
  if (grepl("pdf$", f)) pdf(f, width = 9.0, height = 3.8) else png(f, width = 1390, height = 590, res = 150)
  draw(); invisible(dev.off()); cat("  wrote", f, "\n")
}
