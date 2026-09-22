# FIGURE: "The fold supplies a family; the bound holds across all of it."
# The paper's headline result, which is easy to state and easy to misread as a state selection.
# Left : the band Theta-invariance leaves. Zero width in the deep infrared, opening with momentum.
# Right: the mass every member returns. The least-occupied member is the CEILING, not the pick.
INK<-"#1b1b1b"; BAND<-rgb(0.04,0.43,0.62,0.16); LINE<-"#0b6e9e"; HI<-"#c2461f"; GREY<-"#8a8a8a"
nmin <- function(x) (1-sqrt(1-exp(-x^2)))/2
nmax <- function(x) (1+sqrt(1-exp(-x^2)))/2
draw <- function() {
  par(mfrow=c(1,2), mar=c(4.6,4.9,3.4,1.4), mgp=c(2.8,0.75,0), family="sans",
      cex.axis=0.90, cex.lab=1.00, col.axis=INK, col.lab=INK, fg=INK, xpd=NA)

  ## ---- left: the band
  x <- seq(0,4,length.out=900)
  plot(NA, xlim=c(0,4), ylim=c(0,1.26), axes=FALSE,
       xlab=expression(paste("momentum  ", x)), ylab="occupation per mode")
  polygon(c(x,rev(x)), c(nmin(x),rev(nmax(x))), col=BAND, border=NA)
  lines(x, nmax(x), col=LINE, lwd=1.6, lty=2)
  lines(x, nmin(x), col=LINE, lwd=2.8)
  segments(0,0.5,4,0.5, col=GREY, lty=3)
  axis(1, at=0:4, lwd=0, lwd.ticks=1); axis(2, at=c(0,0.5,1), las=1, lwd=0, lwd.ticks=1)
  # the band fills to 1.0, so the only guaranteed-clear space is ABOVE it
  text(0.05, 1.17, "all states the fold allows", col=LINE, cex=0.90, adj=0)
  text(2.05, 0.125, "least occupied", col=LINE, cex=0.88, adj=0)
  points(0, 0.5, pch=16, cex=0.9, col=INK)
  # this one sits inside the band on the n = 1/2 line it is describing, where the band is
  # wide (at x=0.62 it spans 0.24 to 0.76) so it cannot touch either edge
  text(0.42, 0.575, "half in each component", cex=0.86, col=INK, adj=0)
  mtext("the fold leaves a band", side=3, line=1.3, cex=0.95, col=INK)

  ## ---- right: what each member returns, from band_price.R
  xc <- c(0.10,0.25,0.50,0.644,0.75,1.00,1.50)
  M  <- c(491.6,490.1,470.0,440.0,410.0,330.4,210.7)
  plot(NA, xlim=c(0,1.55), ylim=c(180,520), axes=FALSE,
       xlab=expression(paste("band used,  ", x[c])),
       ylab=expression(paste("dark-matter mass  ", M[1], "  (PeV)")))
  polygon(c(-0.1,1.6,1.6,-0.1), c(180,180,491.6,491.6), col=rgb(0.04,0.43,0.62,0.07), border=NA)
  segments(0,491.6,1.55,491.6, col=HI, lwd=2.4)
  lines(xc, M, col=LINE, lwd=2.8); points(xc, M, pch=16, cex=0.85, col=LINE)
  axis(1, at=c(0,0.5,1,1.5), lwd=0, lwd.ticks=1)
  axis(2, at=c(200,300,400,491.6), labels=c("200","300","400","491.6"), las=1, lwd=0, lwd.ticks=1)
  text(0.02, 511, "no admissible state above", col=HI, cex=0.88, adj=0)
  text(0.03, 238, "more band used\n-> bigger integral\n-> smaller mass", col=LINE, cex=0.84, adj=0)
  mtext("the bound holds across all of it", side=3, line=1.3, cex=0.95, col=INK)
}
dir.create("pub/paper2/figs", showWarnings=FALSE)
pdf("pub/paper2/figs/fig_band.pdf", width=6.5, height=3.45); draw(); dev.off()
png("pub/paper2/figs/fig_band.png", width=1950, height=1035, res=300); draw(); dev.off()
cat("wrote fig_band.{pdf,png}\n")
