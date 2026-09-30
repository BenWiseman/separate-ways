#!/usr/bin/env Rscript
# fig_nariai_convergence.R -- two unrelated quantities reaching one at the same point.
#
# Across the Schwarzschild-de Sitter family, parametrised by 9 Lambda M^2 from 0 to 1:
#   (a) the contact radius of Section 5, as a fraction of the black-hole horizon radius,
#       computed from the angular budget (contact_charged.R section 7);
#   (b) the ratio of the two horizons' surface gravities, kappa_c / kappa_b.
# The first is a causality calculation about where a point can reach its own fold image.
# The second is a statement about whether the two horizons are interchangeable, which is
# what decides whether one involution can serve both (fold_map_classification.R section 5).
# They were computed for different reasons and both reach one only at Nariai.

M <- 1
sds_roots <- function(lam) { z <- polyroot(c(2*M, -1, 0, lam/3))
                             sort(Re(z[abs(Im(z)) < 1e-8])) }        # r_n < 0 < r_b < r_c
Wfac <- function(r, lam, q) (lam/3)*r*(q[2]-r)*(q[3]-r)*(r-q[1])
leg <- function(r, lam) {
  q <- sds_roots(lam); rb <- q[2]; m <- rb/2
  lo <- if (r < m) integrate(function(x) 1/sqrt(Wfac(x,lam,q)), r, m,
                             rel.tol=1e-10, subdivisions=3000L)$value else 0
  hi <- integrate(function(u) { x <- rb - u*u
        2/sqrt((lam/3)*x*(q[3]-x)*(x-q[1])) }, 0, sqrt(rb - max(r,m)),
        rel.tol=1e-10, subdivisions=3000L)$value
  lo + hi }
contact_frac <- function(y) {                                        # y = 9 Lambda M^2
  if (y < 1e-9) return(0.5)
  lam <- y/9; rb <- sds_roots(lam)[2]
  uniroot(function(r) 2*leg(r,lam) - pi, c(1e-9, rb*(1-1e-10)), tol=1e-12)$root / rb }
kappa_ratio <- function(y) {
  if (y < 1e-12) return(0)
  lam <- y/9; q <- sds_roots(lam); rb <- q[2]; rc <- q[3]
  kb <- abs(2*M/rb^2 - 2*lam*rb/3)/2; kc <- abs(2*M/rc^2 - 2*lam*rc/3)/2
  kc/kb }

ys <- c(seq(1e-6, 0.97, length.out = 120), seq(0.972, 0.99995, length.out = 60))
cf <- sapply(ys, contact_frac); kr <- sapply(ys, kappa_ratio)
stopifnot(abs(cf[1] - 0.5) < 1e-6, cf[length(cf)] > 0.98, kr[length(kr)] > 0.98,
          all(diff(cf) > -1e-9), all(diff(kr) > -1e-9))
cat(sprintf("  contact fraction runs %.4f -> %.4f; kappa_c/kappa_b runs %.4f -> %.4f\n",
            cf[1], cf[length(cf)], kr[1], kr[length(kr)]))

ink <- "grey15"; c1 <- "#1f4e79"; c2 <- "#2e7d5b"; mark <- "#a8400f"
draw <- function() {
  par(mar=c(3.9,4.2,1.5,1.2), mgp=c(2.5,0.7,0), xpd=NA)
  plot(NA, xlim=c(0,1), ylim=c(0,1.045), axes=FALSE,
       xlab=expression(paste("9 ", Lambda, italic(M)^2, "     (0 = Schwarzschild, 1 = Nariai)")),
       ylab="ratio")
  axis(1, at=seq(0,1,0.25)); axis(2, at=seq(0,1,0.25), las=1)
  rect(0, 0, 1, 1.045, col="#fbfbfa", border=NA)
  abline(h=1, col="grey70", lwd=1, lty=2)
  lines(ys, cf, col=c1, lwd=2.4)
  lines(ys, kr, col=c2, lwd=2.4)
  text(0.40, contact_frac(0.40)+0.105, "contact radius, over the", col=c1, cex=0.82, adj=0)
  text(0.40, contact_frac(0.40)+0.062, "horizon radius", col=c1, cex=0.82, adj=0)
  text(0.55, kappa_ratio(0.55)-0.085, "the two horizons' surface", col=c2, cex=0.82, adj=0)
  text(0.55, kappa_ratio(0.55)-0.128, "gravities, smaller over larger", col=c2, cex=0.82, adj=0)
  points(1, 1, pch=19, cex=1.3, col=mark)
  segments(1, 0, 1, 1, col=mark, lty=3, lwd=1.2)
  text(0.985, 0.17, "Nariai: the only member of the family", col=mark, cex=0.80, adj=1)
  text(0.985, 0.125, "the fold can fix, and the only one whose", col=mark, cex=0.80, adj=1)
  text(0.985, 0.08, "interior is contact all the way out", col=mark, cex=0.80, adj=1)
  points(0, 0.5, pch=19, cex=1.1, col=c1)
  text(0.022, 0.462, expression(paste("Schwarzschild: ", italic(r) <= italic(M))), col=c1, cex=0.80, adj=0)
}
for (f in c("papers/2_over_the_horizon/fig_companion_nariai.pdf","papers/2_over_the_horizon/fig_companion_nariai.png")) {
  if (grepl("pdf$", f)) pdf(f, width=6.6, height=4.6) else png(f, width=1000, height=700, res=150)
  draw(); invisible(dev.off()); cat("  wrote", f, "\n") }
