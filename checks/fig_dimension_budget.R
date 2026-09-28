#!/usr/bin/env Rscript
# fig_dimension_budget.R -- the budget and the bill.
#
# In D-dimensional Schwarzschild-Tangherlini the maximum angular path a causal curve can
# accumulate on the two interior legs joining a point to its fold image is
#     Dphi_max(r) = (4/n)[pi/2 - arcsin (r/r_h)^{n/2}],   n = D - 3,
# rising to 2 pi/n at the singularity. The antipodal map bills pi, flat in r and in D.
# Contact happens where the curve lies above the line. Verified in contact_dimension.R.

dphi <- function(x, n) (4/n)*(pi/2 - asin(x^(n/2)))      # x = r/r_h
stopifnot(abs(dphi(0.5, 1) - pi) < 1e-12, dphi(1e-9, 2) < pi, dphi(0, 3) < pi)

ink <- "grey15"; bill <- "#a8400f"
cols <- c("#1f4e79", "#2e7d5b", "#8a6d1f", "#7a4a7a")    # D = 4,5,6,7
Ds <- 4:7

draw <- function() {
  par(mar=c(3.9,4.1,1.4,1.2), mgp=c(2.5,0.7,0), xpd=NA)
  x <- seq(0, 1, length.out=800)
  plot(NA, xlim=c(0,1), ylim=c(0, 2*pi), axes=FALSE,
       xlab=expression(paste("depth  ", italic(r)/italic(r)[h], "   (0 = singularity, 1 = horizon)")),
       ylab="angular path available")
  axis(1, at=seq(0,1,0.25), labels=c("0","0.25","0.5","0.75","1"))
  axis(2, at=c(0,pi/2,pi,3*pi/2,2*pi), las=1,
       labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))

  # the bill
  rect(0, pi, 1, 2*pi, col="#f7ece5", border=NA)
  abline(h=pi, col=bill, lwd=2.2)
  text(0.985, pi+0.17, "what the antipodal map asks for", col=bill, cex=0.82, adj=1)

  for (i in seq_along(Ds)) {
    n <- Ds[i]-3
    lines(x, dphi(x, n), col=cols[i], lwd=2.3)
    xl <- c(0.10, 0.16, 0.16, 0.16)[i]
    text(xl, dphi(xl, n) + 0.17, bquote(italic(D)==.(Ds[i])), col=cols[i], cex=0.95, adj=0)
  }

  # the D = 4 crossing
  points(0.5, pi, pch=19, cex=1.25, col=cols[1])
  segments(0.5, 0, 0.5, pi, col=cols[1], lty=3, lwd=1.2)
  text(0.535, 2.62, expression(paste(italic(r)==italic(r)[h]/2)), col=cols[1], cex=0.90, adj=0)
  arrows(0.5, 0.26, 0.015, 0.26, length=0.07, col=cols[1], lwd=1.4)
  text(0.255, 0.50, "the sheets can touch", col=cols[1], cex=0.82)

  points(0, pi, pch=19, cex=1.1, col=cols[2])
  # the D = 5 curve ran through this at pi-0.40; dropped clear of it
  # the D = 5, 6 and 7 curves are stacked close here, so there is no height that clears all
  # of them; masked instead of moved
  source("checks/fig_label.R")
  lab_on(0.055, pi-0.62, "reaches the line only at the singularity", bg = "white",
         col=cols[2], cex=0.78, adj=0)
  # below every curve: D=7, the lowest, only falls past 0.28 beyond r/r_h = 0.95
  text(0.48, 0.28, "never reach the line, at any depth", col=cols[3], cex=0.80, adj=0)
}

for (f in c("paper/fig_companion_dimension.pdf","paper/fig_companion_dimension.png")) {
  if (grepl("pdf$", f)) pdf(f, width=6.6, height=4.7) else png(f, width=1000, height=710, res=150)
  draw(); invisible(dev.off()); cat("  wrote", f, "\n")
}
