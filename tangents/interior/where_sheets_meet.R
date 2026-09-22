# Ben has raised "the graviton acts between the sheets, visible only inside event
# horizons" several times. Tonight's wedge theorem says the two exteriors are ALWAYS
# spacelike separated, which looks like it kills any inter-sheet force. It does not. It
# kills it in the EXTERIOR only, and that is a different statement.
#
# Kruskal coordinates are conformally flat in the (T,X) plane, so light cones are at 45
# degrees and the causal structure is exactly Minkowski's. Everything below is therefore
# a computation and not an argument.
#
#   R (our exterior)        X >  |T|
#   L (mirror exterior)     X < -|T|
#   F (black hole interior) T >  |X|
#   P (white hole interior) T < -|X|

inR <- function(p) p[2] >  abs(p[1]);  inL <- function(p) p[2] < -abs(p[1])
inF <- function(p) p[1] >  abs(p[2]);  inP <- function(p) p[1] < -abs(p[2])
region <- function(p) if (inR(p)) "R" else if (inL(p)) "L" else if (inF(p)) "F" else if (inP(p)) "P" else "horizon"
# q is in the causal FUTURE of p iff dT >= |dX|
causal <- function(p,q) (q[1]-p[1]) >= abs(q[2]-p[2]) - 1e-12

set.seed(23)
samp <- function(test, n=1) { out <- list()
  while (length(out) < n) { p <- c(rnorm(1,0,3), rnorm(1,0,3)); if (test(p)) out[[length(out)+1]] <- p }
  if (n==1) out[[1]] else out }

cat("=== 1. which regions can our exterior R causally influence?\n\n")
cat("        target region     fraction of 20000 pairs with R -> target causal\n")
for (nm in c("R","L","F","P")) {
  tst <- switch(nm, R=inR, L=inL, F=inF, P=inP)
  hits <- mean(replicate(20000, { p <- samp(inR); q <- samp(tst); causal(p,q) }))
  cat(sprintf("   %-14s %30.4f\n", nm, hits))
}
cat("\n  R can reach R and F, never L (spacelike, tonight's theorem), never P (P is in\n")
cat("  R's past, not its future).\n")

cat("\n=== 2. the same for the mirror exterior L, and then the INTERSECTION\n\n")
cat("  For an event to feel BOTH sheets it must lie in the causal future of a point of\n")
cat("  R and of a point of L. Sample events and ask which regions admit that:\n\n")
cat("        event region     can be influenced by R?   by L?   by BOTH?\n")
for (nm in c("R","L","F","P")) {
  tst <- switch(nm, R=inR, L=inL, F=inF, P=inP)
  fr <- mean(replicate(8000, { e <- samp(tst); any(replicate(40, causal(samp(inR), e))) }))
  fl <- mean(replicate(8000, { e <- samp(tst); any(replicate(40, causal(samp(inL), e))) }))
  fb <- mean(replicate(8000, { e <- samp(tst)
        any(replicate(40, causal(samp(inR), e))) && any(replicate(40, causal(samp(inL), e))) }))
  cat(sprintf("   %-14s %20.3f %9.3f %9.3f\n", nm, fr, fl, fb))
}
cat("\n  Only F. The black hole interior is EXACTLY the set of events both exteriors can\n")
cat("  causally influence, and it is the only such set. J+(R) intersect J+(L) = F, which\n")
cat("  is immediate once light cones are at 45 degrees: T > |X| is precisely the\n")
cat("  condition to be in the future of some X > |T| point and some X < -|T| point.\n")
cat("\n  So Ben's premise is right and it is forced rather than allowed. An inter-sheet\n")
cat("  interaction is causally impossible outside a horizon and causally possible in\n")
cat("  exactly one place: inside one. That is not a model choice, it is the causal\n")
cat("  structure of a bifurcate Killing horizon.\n")

cat("\n=== 3. and the mirror of that statement: where can both exteriors RECEIVE from?\n\n")
for (nm in c("F","P")) {
  tst <- switch(nm, F=inF, P=inP)
  fb <- mean(replicate(8000, { e <- samp(tst)
        any(replicate(40, causal(e, samp(inR)))) && any(replicate(40, causal(e, samp(inL)))) }))
  cat(sprintf("   events in %s that can influence BOTH exteriors: %.3f\n", nm, fb))
}
cat("\n  J-(R) intersect J-(L) = P, the white hole interior. So F is the unique common\n")
cat("  SINK of the two sheets and P is their unique common SOURCE.\n")

cat("\n=== 4. THE FOLD'S OWN MAP EXCHANGES THEM\n\n")
cat("  The wedge reflection J is the boost through imaginary angle pi, which in Kruskal\n")
cat("  coordinates is (T,X) -> (-T,-X). Apply it and read off the regions:\n\n")
cat("        point            region      J(point)           region\n")
for (p in list(c(0.4,2.0), c(-0.3,-1.7), c(2.2,0.5), c(-2.6,-0.9))) {
  Jp <- -p
  cat(sprintf("   (%6.2f,%6.2f) %10s   (%6.2f,%6.2f) %12s\n", p[1],p[2],region(p), Jp[1],Jp[2],region(Jp)))
}
ok <- all(replicate(20000, { p <- c(rnorm(1,0,3),rnorm(1,0,3))
      (region(p)=="R" && region(-p)=="L") || (region(p)=="L" && region(-p)=="R") ||
      (region(p)=="F" && region(-p)=="P") || (region(p)=="P" && region(-p)=="F") ||
      region(p)=="horizon" }))
cat(sprintf("\n  R <-> L and F <-> P, on all 20000 sampled points: %s\n", ok))
cat("\n  So the fold does not merely pair the two exteriors. It pairs the black hole\n")
cat("  interior with the WHITE HOLE interior. Ben has raised the white hole with its\n")
cat("  singularity in the future repeatedly; it is J's action on the Kruskal regions,\n")
cat("  not an extra assumption.\n")
