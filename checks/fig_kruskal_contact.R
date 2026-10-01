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
  # Slid a tenth of a unit along its own line, which is where the dotted curve wants it anyway:
  # at the old start the descending horizon came down through the leading "s".
  text(0.598, 0.005, "spacelike: no signal either way", cex=0.78, col=cool, srt=11.6)

  qx <- 0.42; qt <- sqrt(0.90 + qx^2)                    # interior pair, inside the band
  segments(qx, qt, -qx, -qt, col=warm, lwd=1.7)
  points(c(qx,-qx), c(qt,-qt), pch=19, cex=1.2, col=warm)
  text(qx+0.13, qt+0.05, expression(italic(y)), cex=1.15, col=warm)
  text(-qx-0.15, -qt-0.05, expression(Theta*italic(y)), cex=1.15, col=warm)

  text( 1.10, -0.46, "our sheet", cex=0.94, col="grey25")
  text(-1.10,  0.46, "the mirror sheet", cex=0.94, col="grey25")
  # WHERE THE TEXT SITS, and why it is not where it used to be (2026-09-30). Every label below
  # was placed by solving for the clear space rather than by eye, because at r = 1 the interior
  # is only as wide as it is tall and four labels were being written through by the figure's own
  # lines. The three boundaries that matter, in user units:
  #   the horizons          x = +y and x = -y
  #   the singularities     t = +/- sqrt(1 + x^2)
  #   the contact geodesic  x = 0.4048 t, from (0.42, 1.0375) to (-0.42, -1.0375)
  # At height y the future interior runs from -y to +y and the geodesic splits it at 0.4048 y,
  # so the clear part is the LEFT 1.4048 y of it. Everything that has to be a sentence goes
  # there; short labels keep the axis.
  text( 0.00,  1.26,  "singularity", cex=0.86, col=ink)
  text( 0.00, -1.23,  "past singularity", cex=0.86, col=ink)
  # Above the singularity, not on it: the curve rises to t = 1.039 at the ends of this label,
  # which is what used to strike out the word "touch".
  text( 0.00,  1.125, "the sheets can touch", cex=0.78, col=warm)
  text( 0.00,  0.56,  expression(italic(UV) > 0), cex=0.95, col="grey30")
  # These two ran across the whole interior on the axis, so the geodesic crossed one and both
  # horizons crossed the other. Left of the geodesic, at the height the lower line sits, the
  # clear width is 0.939 of a unit against the 0.806 the wider line needs.
  text(-0.199, 0.805, "commutator nonzero, but the", cex=0.78, col="grey30")
  text(-0.199, 0.700, "angle cannot be covered here", cex=0.78, col="grey30")
  text( 0.00, -0.50,  expression(italic(UV) > 0), cex=0.95, col="grey30")
  # Centred 0.085 of a unit from the diagonal it names, measured perpendicular, and on the
  # exterior side of it: about 8 point of white. At the old position the line ran through it.
  text( 0.815, 0.695, "horizon", cex=0.78, col=horiz, srt=45)
  text(-0.815, 0.695, "horizon", cex=0.78, col=horiz, srt=-45)
  # Off the axis to the left, under the warm curve it names and away from the geodesic, which is
  # the same colour and would otherwise be the nearer line.
  text(-0.60,  0.97,  expression(italic(r) == italic(M)), cex=0.84, col=warm)
  points(0, 0, pch=1, cex=1.5, col=warm, lwd=1.6)
  # Clear of the descending horizon, which used to come down through the "t" of the second line.
  text(0.74, -0.22, "and only through here,", cex=0.78, col=warm)
  text(0.74, -0.33, "the bifurcation surface", cex=0.78, col=warm)
}

for (f in c("papers/2_over_the_horizon/fig_companion_kruskal.pdf","papers/2_over_the_horizon/fig_companion_kruskal.png")) {
  if (grepl("pdf$", f)) pdf(f, width=6.4, height=6.0) else png(f, width=960, height=900, res=150)
  draw(); invisible(dev.off()); cat("  wrote", f, "\n")
}
