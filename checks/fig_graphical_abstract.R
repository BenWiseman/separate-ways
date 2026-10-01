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
#
# SIZE (2026-09-30). This was drawn 13.2 inches wide and placed at 100 per cent of a 6.5 inch
# text block, so the page scaled it to 0.49 and every label on it landed under seven point: the
# title at 6.9, the panel headings at 6.4, the labels inside the panels at 3.9. It is now drawn
# at 6.5 inches, which is the text block exactly, so one point on this canvas is one point on
# the page and the sizes below are the sizes a reader gets. par(mfrow) drops par("cex") to 0.66
# on a three-column layout, which is the other half of how the text got so small, so cex is set
# back to 1 straight afterwards and every size is stated in points against a 9 point base.
# Nothing here renders below 8 point.
#
# The sizes are whole points on purpose. R's pdf() device rounds a font to the nearest point
# before it writes it, so 7.6, 7.8, 8.2 and 8.4 all come out of the PDF as 8, while png() honours
# the fraction and comes out a few per cent smaller. Asking for the integer keeps the two files
# the same drawing. To read the sizes back out of the PDF, inflate its content stream and take
# the scale of each Tm: the device writes /F1 1 Tf and puts the size in the text matrix, so
# grepping for Tf returns 1 every time and tells you nothing.

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

WID <- 6.5; HGT <- 3.75; PS <- 9                       # inches, inches, points
pt  <- function(size) size/PS                          # a size in points, as a cex
TITLE<-pt(11); HEAD<-pt(9); SUB<-pt(8); FOOT<-pt(8); CAP<-pt(8)
BIG  <-pt(9);  SHOUT<-pt(8); SMALL<-pt(8); DBIG<-pt(9); DLAB<-pt(8)

XL <- 1.30          # panel 3 runs to u = 0 at x = 1 and keeps the rest for the four D labels.
XD <- 1.04          # Stacked at their own curve ends they cannot be read off the wrong curve.

draw <- function() {
  par(mfrow=c(1,3), mar=c(2.0,1.0,3.2,1.0), oma=c(3.2,0.6,3.0,0.6), xpd=FALSE)
  par(cex=1)        # mfrow has just set it to 0.66. Undo that, or every size above is a lie.

  # 1. silent outside
  plot(NA, xlim=c(0,pi), ylim=c(0,pi), axes=FALSE, xlab="", ylab="")
  rect(0,0,pi,pi,col="#fcfcfb",border=NA)
  polygon(c(0,pi/2,0), c(pi/2,pi,pi), col=pale, border=NA)
  polygon(c(0,pi/2,0), c(pi/2,0,0),   col=pale, border=NA)
  segments(pi,0,pi,pi,col=warm,lwd=4)
  points(0,pi/2,pch=19,cex=2.0,col=blue); points(pi,pi/2,pch=19,cex=2.0,col=warm)
  text(0.13,pi/2+0.34,"a point",col=blue,cex=BIG,adj=0)
  text(pi-0.13,pi/2+0.34,"its image",col=warm,cex=BIG,adj=1)
  arrows(pi/2+0.10,pi/2,pi-0.12,pi/2,length=0.07,col=ink,lwd=1.6)
  text(pi/2-0.06,pi/2,"never reaches",col=ink,cex=SMALL,adj=1)
  mtext("OUTSIDE any horizon",side=3,line=1.75,cex=HEAD,col=ink,font=2)
  mtext("the sheets are silent, always",side=3,line=0.55,cex=SUB,col=ink)
  mtext("light cones fall short of the antipode",side=1,line=0.75,cex=FOOT,col="#5b6b7d")

  # 2. open inside, in the inner half
  x <- seq(0,1,length.out=400)
  plot(NA, xlim=c(0,1), ylim=c(0,2*pi), axes=FALSE, xlab="", ylab="")
  rect(0,0,1,2*pi,col="#fcfcfb",border=NA)
  rect(0.5,0,1,2*pi,col=sand,border=NA)
  segments(0,pi,1,pi,col=warm,lwd=3)
  lines(x,turn(x),col=blue,lwd=4)
  points(0.5,pi,pch=19,cex=2.0,col=warm)
  text(0.50,pi-0.95,"THE SHEETS",col=warm,cex=SHOUT,adj=0,font=2)
  text(0.50,pi-1.55,"TOUCH HERE",col=warm,cex=SHOUT,adj=0,font=2)
  # Two lines, because at 6.5 inches the panel is no longer wide enough for one: on a single
  # line it runs to x = 0.65 and the causal curve is at x = 0.59 at the height it has to sit
  # at. Panel 3 is worse still. The words are the words they always were.
  text(0.02,pi+0.80,"what the fold",col=warm,cex=SMALL,adj=0)
  text(0.02,pi+0.40,"asks for",col=warm,cex=SMALL,adj=0)
  text(0.31,2*pi-0.38,"what causality allows",col=blue,cex=SMALL,adj=0)
  mtext("INSIDE a hole with a past",side=3,line=1.75,cex=HEAD,col=ink,font=2)
  mtext("its inner half is reached at any charge",side=3,line=0.55,cex=SUB,col=ink)
  mtext("horizon on the left, singularity on the right",side=1,line=0.75,cex=FOOT,col="#5b6b7d")

  # 3. only four dimensions
  u <- seq(1e-6,1,length.out=400)
  plot(NA, xlim=c(0,XL), ylim=c(0,2*pi*1.06), axes=FALSE, xlab="", ylab="")
  rect(0,0,XL,2*pi*1.06,col="#fcfcfb",border=NA)   # to XL, or the panel grows a seam at x = 1
  segments(0,pi,1,pi,col=warm,lwd=3)
  cols <- c(blue,"#2e7d5b","#8a7a1f","#7a4a7a")
  for (i in 1:4) lines(rev(u), budget(u,i+3), col=cols[i], lwd=if(i==1) 4 else 2.4)
  # FOUR curves are drawn and four are now named. D = 7 had no label at all, and D = 6 was set
  # 0.42 below its own curve, which put it on top of the unlabelled D = 7 one: a reader took
  # the purple curve for D = 6. Each label now sits at the exact height its own curve ends at,
  # immediately to the right of that end, in that curve's colour.
  for (i in 1:4)
    text(XD, budget(1e-6,i+3), col=cols[i], adj=0, cex=if(i==1) DBIG else DLAB,
         labels=switch(i, expression(bold(italic(D)==4)), expression(italic(D)==5),
                          expression(italic(D)==6), expression(italic(D)==7)))
  text(0.04,pi+0.80,"what the fold",col=warm,cex=SMALL,adj=0)
  text(0.04,pi+0.40,"asks for",col=warm,cex=SMALL,adj=0)
  mtext("and ONLY in four dimensions",side=3,line=1.75,cex=HEAD,col=ink,font=2)
  mtext("five touches the line and never crosses",side=3,line=0.55,cex=SUB,col=ink)
  mtext("nothing was tuned to arrange that",side=1,line=0.75,cex=FOOT,col="#5b6b7d")

  # A label for the figure, not the paper's title. The title sits directly above this on page 1
  # and repeating it there reads as filler; the figure still stands alone if it is lifted out.
  mtext("Where the two sheets can touch, and where they cannot",
        outer=TRUE, side=3, line=1.1, cex=TITLE, col=ink, font=2)
  # One line of this ran 9.2 inches at the size it needs to be read at, so it is two.
  mtext("A CPT-folded universe is sealed at every horizon by a theorem, and opens in one place only:",
        outer=TRUE, side=1, line=0.45, cex=CAP, col=ink)   # on side 1 a LARGER line sits LOWER,
  mtext("inside a black hole that has a past, in the inner half of an uncharged one.",
        outer=TRUE, side=1, line=1.55, cex=CAP, col=ink)   # so the second sentence takes it
}
for (f in c("papers/2_over_the_horizon/fig_companion_graphical_abstract.pdf",
            "papers/2_over_the_horizon/fig_companion_graphical_abstract.png")) {
  if (grepl("pdf$", f)) pdf(f, width=WID, height=HGT, pointsize=PS)
  else png(f, width=round(WID*300), height=round(HGT*300), res=300, pointsize=PS)
  draw(); invisible(dev.off()); cat("  wrote", f, "\n")
}
