# FIGURE: "The bang is an avoided crossing, and one sheet gets half of it."
# The paper's central mechanism, drawn so a skimming reader sees the argument in one look.
# Left : the two levels sweeping past one another, and which part of the sweep each cosmology has.
# Right: what each produces. Curves are closed forms, not fits.
#   two sheets (adopted): n = (1 - sqrt(1 - exp(-x^2)))/2   -> Gaussian tail, finite energy
#   one sheet  (bang-adiabatic): illustrative 1/2 / (1 + (x/a)^4) -> p^-4 tail, log-divergent energy
INK <- "#1b1b1b"; TWO <- "#0b6e9e"; ONE <- "#c2461f"; GREY <- "#b0b0b0"; WASH <- rgb(0.76,0.27,0.11,0.06)
draw <- function() {
  par(mfrow=c(1,2), mar=c(4.6,4.8,3.4,2.2), mgp=c(2.7,0.75,0), family="sans",
      cex.axis=0.90, cex.lab=1.00, col.axis=INK, col.lab=INK, fg=INK, xpd=NA)

  ## ---- panel 1: the crossing
  eta <- seq(-3,3,length.out=900); p <- 0.5
  Ep <- sqrt(eta^2+p^2)
  plot(NA, xlim=c(-3,3), ylim=c(-5.2,4.4), axes=FALSE,
       xlab=expression(paste("conformal time  ", eta)), ylab="mode energy")
  rect(-3, -5.2, 0, 4.4, col=WASH, border=NA)
  axis(1, at=c(-3,0,3), labels=c("before","bang","after"), lwd=0, lwd.ticks=1)
  axis(2, at=c(-3,0,3), las=1, lwd=0, lwd.ticks=1)
  lines(eta, eta, col=GREY, lty=3); lines(eta, -eta, col=GREY, lty=3)
  lines(eta, Ep, col=INK, lwd=2.8); lines(eta, -Ep, col=INK, lwd=2.8)
  segments(0,-p,0,p, col=INK, lwd=1.1)
  text(0.30, 0.02, "2p", cex=0.85, col=INK, adj=0)
  # sweeps
  arrows(-2.75, 2.72, 2.75, -2.72, col=TWO, lwd=2.4, length=0.085)
  arrows( 0.10,-0.05, 2.75, -2.72, col=ONE, lwd=2.4, length=0.085)
  # labels live in the empty bands beyond |E| = 3.04, never on a curve or an arrow
  text(-2.95, 3.95, "two sheets: the whole crossing", col=TWO, cex=0.90, adj=0)
  # the two red notes are on SEPARATE lines and on opposite sides, so neither can reach
  # the other: the left block ends near x=-0.6 and the right block starts near x=+0.7
  # three separate LINES, so width estimates cannot make them collide. The lower
  # hyperbola bottoms out at -3.04, so every one of these sits below the drawing.
  text( 2.95,-3.40, "one sheet: half a sweep", col=ONE, cex=0.90, adj=1)
  text(-2.95,-4.15, "this half does not exist", col=ONE, cex=0.84, adj=0)
  text(-2.95,-4.72, "for a one-sheet universe", col=ONE, cex=0.84, adj=0)
  mtext("every mode crosses an avoided crossing", side=3, line=1.3, cex=0.95, col=INK)

  ## ---- panel 2: the two outcomes
  x <- seq(0.02, 4.2, length.out=900)
  n_two <- (1-sqrt(1-exp(-x^2)))/2
  n_one <- 0.5/(1+(x/0.62)^4)
  plot(NA, xlim=c(0,4.2), ylim=c(1e-7,1), log="y", axes=FALSE,
       xlab=expression(paste("momentum  ", x)), ylab="occupation per mode")
  axis(1, at=0:4, lwd=0, lwd.ticks=1)
  axis(2, at=10^c(-7,-5,-3,-1), labels=expression(10^-7,10^-5,10^-3,10^-1), las=1,
       lwd=0, lwd.ticks=1)
  lines(x, n_one, col=ONE, lwd=2.8, lty=2)
  lines(x, n_two, col=TWO, lwd=2.8)
  # the red curve at x=2.6 is 1.4e-3 and the blue 3e-4, so this group sits well above both
  text(2.45, 2.6e-1, "one sheet", col=ONE, cex=0.92, adj=0)
  text(2.45, 6.0e-2, expression(paste("tail ", ~~ p^-4)), col=ONE, cex=0.84, adj=0)
  text(2.45, 1.4e-2, "energy diverges", col=ONE, cex=0.84, adj=0)
  # and this one far below both, where nothing is drawn
  text(0.12, 4.0e-5, "two sheets", col=TWO, cex=0.92, adj=0)
  text(0.12, 9.0e-6, "Gaussian tail", col=TWO, cex=0.84, adj=0)
  text(0.12, 2.0e-6, "energy finite", col=TWO, cex=0.84, adj=0)
  mtext("the two sweeps leave different states", side=3, line=1.3, cex=0.95, col=INK)
}
dir.create("pub/paper2/figs", showWarnings=FALSE)
pdf("pub/paper2/figs/fig_crossing.pdf", width=6.5, height=3.45); draw(); dev.off()
png("pub/paper2/figs/fig_crossing.png", width=1950, height=1035, res=300); draw(); dev.off()
cat("wrote fig_crossing.{pdf,png}\n")
