# FIGURE: "The fold." The picture a reader should see before any algebra.
#
# Two copies of the universe meet at the bang. Each has its own local arrow of time, pointing
# AWAY from the bang. CPT relates them. Nothing travels between them.
#
# Authored at 6.5 in, the width the manuscript places it at, so no type is shrunk on placement.
# Every label is positioned by hand with clearance checked by eye; nothing is auto-placed.

INK  <- "#1b1b1b"; TWO <- "#0b6e9e"; ONE <- "#c2461f"; GREY <- "#b0b0b0"
SKYU <- rgb(0.043,0.431,0.620, 0.10)     # our sheet wash
SKYD <- rgb(0.761,0.275,0.122, 0.10)     # mirror sheet wash

draw <- function() {
  par(mar=c(0.3,0.3,0.3,0.3), family="sans", xpd=NA)
  plot(NA, xlim=c(-1.34,1.52), ylim=c(-1.12,1.12), axes=FALSE, xlab="", ylab="", asp=1)

  # --- the two cones. Straight edges so they read as light cones.
  W <- 0.92; H <- 1.00
  polygon(c(0,-W, W), c(0, H, H), col=SKYU, border=TWO, lwd=2.0)
  polygon(c(0,-W, W), c(0,-H,-H), col=SKYD, border=ONE, lwd=2.0)

  # --- contents: "galaxies". Kept to the RIGHT of the time arrow and inside the cone edge,
  # which at height y sits at x = W*|y|/H, so nothing can stray over a boundary.
  # 2026-09-23: this used to call gal(+1) then gal(-1), each drawing its own runif, so the two
  # sheets carried INDEPENDENT random scatters. The caption says the far sheet is our universe
  # mirrored, and the picture said otherwise. Ben caught it. The positions are drawn once and
  # the mirror sheet is the same set reflected in time.
  set.seed(4)
  ys <- seq(0.34, 0.92, length.out=13) + runif(13, -0.02, 0.02)
  xs <- vapply(ys, function(y) { edge <- W*abs(y)/H; runif(1, 0.10*edge, 0.80*edge) }, 0)
  gal <- function(sgn, col) {
    for (i in seq_along(ys)) {
      y <- sgn*ys[i]; x <- xs[i]
      r <- 0.016 + 0.018*abs(y)
      points(x, y, pch=16, cex=r*44, col=adjustcolor(col, 0.45))
      points(x, y, pch=16, cex=r*17, col=adjustcolor(col, 0.95))
    }
  }
  gal(+1, TWO); gal(-1, ONE)

  # --- the bang
  points(0, 0, pch=16, cex=1.9, col=INK)
  points(0, 0, pch=1,  cex=3.4, col=INK, lwd=1.4)
  text(0.10, 0.002, "the bang", col=INK, cex=0.95, adj=0, font=2)

  # --- arrows of time. x = -0.25 is inside the cone for |y| > 0.27, so the shaft never
  # crosses an edge; the label sits beyond the head rather than on the shaft.
  ax <- -0.25
  arrows(ax, 0.36, ax, 0.80, length=0.10, lwd=2.4, col=TWO)
  text(ax, 0.88, "time", col=TWO, cex=0.95, adj=0.5)
  arrows(ax,-0.36, ax,-0.80, length=0.10, lwd=2.4, col=ONE)
  text(ax,-0.88, "time", col=ONE, cex=0.95, adj=0.5)

  # --- sheet labels, placed OUTSIDE the cone: at y=0.60 the edge is at x=0.55, so 0.86 is clear
  text(0.86, 0.62, "our sheet",             col=TWO, cex=1.06, adj=0, font=2)
  text(0.86, 0.50, "stars, us, now",        col=TWO, cex=0.90, adj=0)
  text(0.86,-0.50, "the mirror sheet",      col=ONE, cex=1.06, adj=0, font=2)
  text(0.86,-0.62, "every charge reversed", col=ONE, cex=0.90, adj=0)

  # --- the relation, on the left where the cones are narrow
  text(-1.28, 0.14, "CPT relates",     col=INK,  cex=0.98, adj=0)
  text(-1.28, 0.02, "the two sheets",  col=INK,  cex=0.98, adj=0)
  text(-1.28,-0.14, "nothing travels", col=GREY, cex=0.88, adj=0)
  text(-1.28,-0.25, "between them",    col=GREY, cex=0.88, adj=0)

  # --- the one sentence a skimming reader should leave with
  text(0, 1.08, "Both clocks run away from the bang. Neither sheet is the other's past.",
       col=INK, cex=1.00, adj=0.5)
}

dir.create("pub/paper2/figs", showWarnings=FALSE, recursive=TRUE)
pdf("pub/paper2/figs/fig_fold.pdf", width=6.5, height=5.2); draw(); dev.off()
png("pub/paper2/figs/fig_fold.png", width=1950, height=1560, res=300); draw(); dev.off()
cat("wrote fig_fold.{pdf,png}\n")
