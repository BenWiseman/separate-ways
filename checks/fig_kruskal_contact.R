#!/usr/bin/env Rscript
# fig_kruskal_contact.R -- the companion's one geometric claim, drawn.
#
# Kruskal diagram of an eternal hole, U = T - X and V = T + X, so UV = T^2 - X^2.
# The fold's wedge reflection J is (T,X) -> (-T,-X): it swaps the two exteriors and
# swaps the two interiors. Under it dU = 2U and dV = 2V, so the cross-sheet commutator
# goes as UV and vanishes identically wherever UV < 0, which is both exteriors.
#
# That is the two-dimensional statement. In four dimensions the fold also sends the unit
# vector n to -n, so a point and its image are ANTIPODAL on the sphere and a connecting
# curve must cover pi of angle. Causality caps the angular path at
#     Dphi_max(r) = 2 pi - 4 asin sqrt(r/2M),
# which reaches pi at r = M. So contact reaches the inner half of the hole, and for an
# uncharged one it is confined there; at charge the exclusion is open, as A.15 records. And
# the boundary drawn here is exact, not schematic. A.15 and pinch_4d_exact.R.

M <- 1
UV       <- function(r) (1 - r/(2*M))*exp(r/(2*M))
dphi_max <- function(r) 2*pi - 4*asin(sqrt(r/(2*M)))
r_c      <- 2*M*sin(pi/4)^2
UV_c     <- UV(r_c)
stopifnot(abs(r_c - M) < 1e-12, abs(dphi_max(r_c) - pi) < 1e-12,
          abs(UV_c - sqrt(exp(1))/2) < 1e-12)
cat(sprintf("  contact for r <= %.4f M, i.e. UV >= %.4f, singularity at UV = 1\n", r_c/M, UV_c))

L <- 1.34
ink <- "grey15"; horiz <- "grey45"; fill <- "#dce6f2"; band <- "#f2c9a8"
warm <- "#a8400f"; cool <- "#1f4e79"

hyp <- function(c, s) { x <- seq(-L, L, length.out=600)
                        t <- s*sqrt(c + x^2); k <- abs(t) <= L; list(x=x[k], t=t[k]) }

draw <- function() {
  par(mar=c(0.2,0.2,0.2,0.2), xpd=NA)
  plot(NA, xlim=c(-L,L), ylim=c(-L,L), asp=1, axes=FALSE, xlab="", ylab="")

  xs <- seq(-L, L, length.out=500)
  for (s in c(1,-1)) {                                   # UV > 0: both interiors
    tt <- s*sqrt(1 + xs^2); keep <- abs(tt) <= L
    polygon(c(xs[keep], rev(xs[keep])), c(tt[keep], rev(s*abs(xs[keep]))), col=fill, border=NA)
  }
  for (s in c(1,-1)) {                                   # the contact band, r <= M
    a <- hyp(UV_c, s); b <- hyp(1, s)
    polygon(c(a$x, rev(b$x)), c(a$t, rev(b$t)), col=band, border=NA)
    lines(a$x, a$t, col=warm, lwd=1.8)
  }

  segments(-L,-L, L, L, col=horiz, lwd=1.4)
  segments(-L, L, L,-L, col=horiz, lwd=1.4)
  for (s in c(1,-1)) { b <- hyp(1, s); lines(b$x, b$t, lwd=2.8, col=ink) }

  px <- 0.98; pt <- 0.20                                 # exterior pair: spacelike
  segments(px, pt, -px, -pt, col=cool, lty=3, lwd=1.5)
  points(c(px,-px), c(pt,-pt), pch=19, cex=1.2, col=cool)
  text(px+0.03, pt+0.17, expression(italic(x)), cex=1.15, col=cool)
  text(-px-0.03, -pt-0.20, expression(Theta*italic(x)), cex=1.15, col=cool)
  text(0.50, -0.015, "spacelike: no signal either way", cex=0.78, col=cool, srt=11.6)

  qx <- 0.42; qt <- sqrt(0.90 + qx^2)                    # interior pair, inside the band
  segments(qx, qt, -qx, -qt, col=warm, lwd=1.7)
  points(c(qx,-qx), c(qt,-qt), pch=19, cex=1.2, col=warm)
  text(qx+0.13, qt+0.05, expression(italic(y)), cex=1.15, col=warm)
  text(-qx-0.15, -qt-0.05, expression(Theta*italic(y)), cex=1.15, col=warm)

  text( 1.10, -0.46, "our sheet", cex=0.94, col="grey25")
  text(-1.10,  0.46, "the mirror sheet", cex=0.94, col="grey25")
  text( 0.00,  1.21, "singularity", cex=0.86, col=ink)
  text( 0.00, -1.23, "past singularity", cex=0.86, col=ink)
  text( 0.00,  1.075, "the sheets can touch", cex=0.78, col=warm)
  text( 0.00,  0.56, expression(italic(UV) > 0), cex=0.95, col="grey30")
  # These two are wider than the interior region is at their height, so the geodesic crosses
  # them wherever they sit. A mask was tried and looked worse: the labels straddle the edge of
  # the shaded region, so a single-colour patch shows as a box. They are on the accepted list
  # in checks/ink_overlaps_accepted.tsv instead, the crossing line being a thin one.
  text( 0.00,  0.42, "commutator nonzero, but the", cex=0.78, col="grey30")
  text( 0.00,  0.30, "angle cannot be covered here", cex=0.78, col="grey30")
  text( 0.00, -0.50, expression(italic(UV) > 0), cex=0.95, col="grey30")
  text( 0.72,  0.79, "horizon", cex=0.78, col=horiz, srt=45)
  text(-0.72,  0.79, "horizon", cex=0.78, col=horiz, srt=-45)
  text( 0.00,  0.845, expression(italic(r) == italic(M)), cex=0.84, col=warm)
  points(0, 0, pch=1, cex=1.5, col=warm, lwd=1.6)
  text(0.62, -0.22, "and only through here,", cex=0.78, col=warm)
  text(0.62, -0.33, "the bifurcation surface", cex=0.78, col=warm)
}

for (f in c("papers/2_over_the_horizon/fig_companion_kruskal.pdf","papers/2_over_the_horizon/fig_companion_kruskal.png")) {
  if (grepl("pdf$", f)) pdf(f, width=6.4, height=6.0) else png(f, width=960, height=900, res=150)
  draw(); invisible(dev.off()); cat("  wrote", f, "\n")
}
