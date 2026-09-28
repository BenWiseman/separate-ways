#!/usr/bin/env Rscript
# fig_graphical_abstract.R -- the companion's argument in three panels, read at thumbnail size.
#
# Same three calculations as Figures 1, 5 and A.15, drawn with the detail stripped out:
#   1. outside any horizon a point and its fold image are spacelike, so the sheets are silent;
#   2. inside a hole with a past the separation turns timelike, and the angular budget puts the
#      contact region in the inner HALF of the interior;
#   3. the same budget in D dimensions clears the bill only in four.
# Nothing here is new arithmetic. It restates fig_desitter_silence, fig_contact_universal and
# fig_dimension_budget at a size a prize committee or an editor reads first.

turn  <- function(x) 4 * asin(sqrt(x))                 # two legs, depth fraction x
budget<- function(u, D) 2*pi/(D-3) - (4/(D-3))*asin(sqrt(u^(D-3)))   # u = r/r_h
stopifnot(abs(turn(0.5) - pi) < 1e-12)
for (D in 4:7) stopifnot(abs(budget(1, D)) < 1e-12)
cat("=== graphical abstract, the three statements it draws ===\n\n")
cat(sprintf("  contact opens at depth %.3f of the interior, where the two-leg budget reaches pi\n", 0.5))
cat(sprintf("  budget at the singularity: D=4 %.4f, D=5 %.4f, D=6 %.4f, bill %.4f\n",
            budget(1e-9,4), budget(1e-9,5), budget(1e-9,6), pi))
stopifnot(budget(1e-9,4) > pi, abs(budget(1e-9,5) - pi) < 1e-3, budget(1e-9,6) < pi)

ink<-"grey12"; blue<-"#1f4e79"; warm<-"#a8400f"; pale<-"#dde5ee"; sand<-"#f2e6dc"
draw <- function() {
  par(mfrow=c(1,3), mar=c(2.2,1.0,3.4,1.0), oma=c(2.6,0.6,3.4,0.6), xpd=FALSE)

  # 1. silent outside
  plot(NA, xlim=c(0,pi), ylim=c(0,pi), axes=FALSE, xlab="", ylab="")
  rect(0,0,pi,pi,col="#fcfcfb",border=NA)
  polygon(c(0,pi/2,0), c(pi/2,pi,pi), col=pale, border=NA)
  polygon(c(0,pi/2,0), c(pi/2,0,0),   col=pale, border=NA)
  segments(pi,0,pi,pi,col=warm,lwd=4)
  points(0,pi/2,pch=19,cex=2.0,col=blue); points(pi,pi/2,pch=19,cex=2.0,col=warm)
  text(0.13,pi/2+0.30,"a point",col=blue,cex=1.15,adj=0)
  text(pi-0.13,pi/2+0.30,"its image",col=warm,cex=1.15,adj=1)
  arrows(pi/2+0.10,pi/2,pi-0.12,pi/2,length=0.10,col=ink,lwd=1.6)
  text(pi/2-0.06,pi/2,"never reaches",col=ink,cex=1.05,adj=1)
  mtext("OUTSIDE any horizon",side=3,line=1.5,cex=1.05,col=ink,font=2)
  mtext("the sheets are silent, always",side=3,line=0.25,cex=1.0,col=ink)
  mtext("light cones fall short of the antipode",side=1,line=0.7,cex=0.95,col="#5b6b7d")

  # 2. open inside, in the inner half
  x <- seq(0,1,length.out=400)
  plot(NA, xlim=c(0,1), ylim=c(0,2*pi), axes=FALSE, xlab="", ylab="")
  rect(0,0,1,2*pi,col="#fcfcfb",border=NA)
  rect(0.5,0,1,2*pi,col=sand,border=NA)
  segments(0,pi,1,pi,col=warm,lwd=3)
  lines(x,turn(x),col=blue,lwd=4)
  points(0.5,pi,pch=19,cex=2.0,col=warm)
  text(0.52,pi-1.05,"THE SHEETS",col=warm,cex=1.15,adj=0,font=2)
  text(0.52,pi-1.72,"TOUCH HERE",col=warm,cex=1.15,adj=0,font=2)
  text(0.03,pi+0.55,"what the fold asks for",col=warm,cex=1.0,adj=0)
  text(0.53,2*pi-0.38,"what causality allows",col=blue,cex=1.0,adj=0)
  mtext("INSIDE a hole with a past",side=3,line=1.5,cex=1.05,col=ink,font=2)
  mtext("its inner half is reached at any charge",side=3,line=0.25,cex=1.0,col=ink)
  mtext("horizon on the left, singularity on the right",side=1,line=0.7,cex=0.95,col="#5b6b7d")

  # 3. only four dimensions
  u <- seq(1e-6,1,length.out=400)
  plot(NA, xlim=c(0,1), ylim=c(0,2*pi*1.06), axes=FALSE, xlab="", ylab="")
  rect(0,0,1,2*pi*1.06,col="#fcfcfb",border=NA)
  segments(0,pi,1,pi,col=warm,lwd=3)
  cols <- c(blue,"#2e7d5b","#8a7a1f","#7a4a7a")
  for (i in 1:4) lines(rev(u), budget(u,i+3), col=cols[i], lwd=if(i==1) 4 else 2.4)
  text(0.92,budget(1e-6,4)-0.18,expression(bold(italic(D)==4)),col=cols[1],cex=1.25,adj=1)
  text(0.92,pi-0.55,expression(italic(D)==5),col=cols[2],cex=1.05,adj=1)
  text(0.92,budget(1e-6,6)-0.42,expression(italic(D)==6),col=cols[3],cex=1.05,adj=1)
  text(0.06,pi+0.45,"what the fold asks for",col=warm,cex=1.0,adj=0)
  mtext("and ONLY in four dimensions",side=3,line=1.5,cex=1.05,col=ink,font=2)
  mtext("five touches the line and never crosses",side=3,line=0.25,cex=1.0,col=ink)
  mtext("nothing was tuned to arrange that",side=1,line=0.7,cex=0.95,col="#5b6b7d")

  # A label for the figure, not the paper's title. The title sits directly above this on page 1
  # and repeating it there reads as filler; the figure still stands alone if it is lifted out.
  mtext("Where the two sheets can touch, and where they cannot",
        outer=TRUE, side=3, line=1.5, cex=1.2, col=ink, font=2)
  mtext("A CPT-folded universe is sealed at every horizon by a theorem, and opens in one place only: inside a black hole that has a past, in the inner half of an uncharged one.",
        outer=TRUE, side=1, line=0.6, cex=0.98, col=ink)
}
for (f in c("paper/fig_companion_graphical_abstract.pdf",
            "paper/fig_companion_graphical_abstract.png")) {
  if (grepl("pdf$", f)) pdf(f, width=13.2, height=4.9) else png(f, width=1900, height=706, res=144)
  draw(); invisible(dev.off()); cat("  wrote", f, "\n")
}
