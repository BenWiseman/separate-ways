# FIGURE: "Two classes of black hole, and why it shows up in how fast they grow."
#
# Section 4.3 divides holes into one-sided (formed by collapse) and two-sided (a past
# singularity, two exteriors), then makes an accretion argument that reads as a non-sequitur
# without a picture. Draw the contrast, because the contrast IS the claim.
#
# Authored at 6.5 in, the width the manuscript places it at.

INK <- "#1b1b1b"; TWO <- "#0b6e9e"; ONE <- "#c2461f"; GREY <- "#8a8a8a"
GAS <- rgb(0.55,0.55,0.55,0.30)

hole <- function(cx, r=0.30) {
  ang <- seq(0, 2*pi, length.out=200)
  polygon(cx + r*cos(ang), r*sin(ang), col=INK, border=NA)
  lines(cx + (r*1.13)*cos(ang), (r*1.13)*sin(ang), col=INK, lwd=1.4, lty=3)
}

# Gas comes IN at shallow angles, light goes OUT at steep ones. Separating them by ANGLE
# rather than by position means nothing can congest the space between the two panels.
feed <- function(cx, from, col) {
  s <- if (from == "left") -1 else 1
  for (a in c(22, -22)) {
    th <- a*pi/180
    arrows(cx + s*1.02*cos(th), 1.02*sin(th),
           cx + s*0.44*cos(th), 0.44*sin(th), length=0.075, lwd=2.3, col=col)
  }
}
rad <- function(cx, to, col) {
  s <- if (to == "left") -1 else 1
  for (a in c(66, -66)) {
    th <- a*pi/180
    segments(cx + s*0.36*cos(th), 0.36*sin(th),
             cx + s*0.84*cos(th), 0.84*sin(th), col=col, lwd=1.9)
    points(cx + s*0.84*cos(th), 0.84*sin(th), pch=18, cex=0.75, col=col)
  }
}

draw <- function() {
  par(mar=c(0.2,0.2,0.2,0.2), family="sans", xpd=NA)
  plot(NA, xlim=c(-2.40,2.40), ylim=c(-1.90,1.42), axes=FALSE, xlab="", ylab="", asp=1)
  L <- -1.20; R <- 1.20

  hole(L); feed(L, "right", GREY); rad(L, "right", TWO)
  text(L, 1.30, "formed by collapse", col=INK, cex=1.02, adj=0.5, font=2)
  text(L, 1.14, "one exterior", col=GREY, cex=0.90, adj=0.5)
  text(L, -1.22, expression(kappa == 1), col=INK, cex=1.12, adj=0.5)

  hole(R); feed(R, "right", GREY); feed(R, "left", GREY)
  rad(R, "right", TWO); rad(R, "left", ONE)
  text(R, 1.30, "a past singularity", col=INK, cex=1.02, adj=0.5, font=2)
  text(R, 1.14, "two exteriors, one mass", col=GREY, cex=0.90, adj=0.5)
  text(R, -1.22, expression(kappa == 1 + dot(M)[other] / dot(M)[ours]),
       col=INK, cex=1.12, adj=0.5)

  # legend on its own line, then the sentence on its own line below it
  arrows(-1.62, -1.60, -1.36, -1.60, length=0.06, lwd=2.2, col=GREY)
  text(-1.30, -1.60, "gas in", col=GREY, cex=0.86, adj=0)
  segments(-0.70, -1.60, -0.46, -1.60, col=TWO, lwd=1.9)
  points(-0.46, -1.60, pch=18, cex=0.7, col=TWO)
  text(-0.40, -1.60, "light out, each side only", col=INK, cex=0.86, adj=0)

  text(0, -1.86,
       "Gravity is set by the shared mass, so each side may accrete at the full limit for the total.",
       col=INK, cex=0.94, adj=0.5)
}

dir.create("pub/paper2/figs", showWarnings=FALSE, recursive=TRUE)
pdf("pub/paper2/figs/fig_twosided.pdf", width=6.5, height=4.5); draw(); dev.off()
png("pub/paper2/figs/fig_twosided.png", width=1950, height=1350, res=300); draw(); dev.off()
cat("wrote fig_twosided.{pdf,png}\n")
