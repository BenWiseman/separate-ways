# Named last turn: does 5.2's tree-level channel ratio h nu : Z nu : W l = 1:1:2 constrain
# the flavour ratio at Earth? If yes, flavour joins direction as a discriminant and F-AG's
# 58-event requirement drops.
#
# First the structural point, because it decides the answer before any numerics. The 1:1:2
# ratio is a statement about WHICH FINAL STATE, fixed by gauge and Higgs structure. For a
# given flavour alpha, ALL THREE channels carry the same Yukawa factor |y_alpha|^2. So the
# channel ratio is flavour-blind and the flavour ratio is set by the Yukawa column, which
# 5.2 leaves free. The two are independent.
#
# The numerics below test what that independence costs: with a free source flavour, how
# much of the Earth-side flavour space can the decay component reach? If it reaches
# everything an astrophysical source can, flavour cannot discriminate.

# oscillation-averaged transfer: P_ab = sum_i |U_ai|^2 |U_bi|^2
pmns <- function(t12, t23, t13, dcp=0) {
  s12<-sin(t12); c12<-cos(t12); s23<-sin(t23); c23<-cos(t23); s13<-sin(t13); c13<-cos(t13)
  d <- exp(-1i*dcp)
  matrix(c( c12*c13,                      s12*c13,                      s13*Conj(d),
           -s12*c23-c12*s23*s13*d,        c12*c23-s12*s23*s13*d,        s23*c13,
            s12*s23-c12*c23*s13*d,       -c12*s23-s12*c23*s13*d,        c23*c13), 3,3, byrow=TRUE) }
Pmat <- function(U) { A <- Mod(U)^2; A %*% t(A) }

d2r <- pi/180
U <- pmns(33.4*d2r, 49.0*d2r, 8.6*d2r)
P <- Pmat(U)
cat("=== 1. the oscillation-averaged transfer matrix (stated angles, not fitted here)\n\n")
cat("   theta12 = 33.4, theta23 = 49.0, theta13 = 8.6 degrees\n\n")
for (i in 1:3) cat(sprintf("   %8.4f %8.4f %8.4f\n", P[i,1], P[i,2], P[i,3]))
cat(sprintf("\n   rows sum to %s (unitarity check)\n", paste(sprintf('%.6f', rowSums(P)), collapse=' ')))

cat("\n=== 2. what the decay component can reach with a FREE source flavour\n\n")
cat("  The source simplex maps LINEARLY under P, so the reachable set is exactly the\n")
cat("  convex hull of P's rows - the images of pure-e, pure-mu and pure-tau sources.\n")
cat("  A first pass sampled the simplex by normalising uniforms and took min/max per\n")
cat("  component. That sampling almost never reaches a vertex, so it UNDERSTATED the\n")
cat("  reachable range and wrongly reported two astrophysical classes as unreachable.\n")
cat("  The hull test below is exact.\n\n")
cat("        pure source       ->   Earth-side (e, mu, tau)\n")
for (i in 1:3) cat(sprintf("   %-18s %10.4f %8.4f %8.4f\n",
    c("nu_e only","nu_mu only","nu_tau only")[i], P[i,1], P[i,2], P[i,3]))
inhull <- function(pt) {
  # solve pt = w %*% P with w >= 0, sum w = 1  (2 free params, exact)
  A <- rbind(t(P[,1:2]) , c(1,1,1))
  b <- c(pt[1:2], 1)
  w <- tryCatch(solve(A, b), error=function(e) rep(NA,3))
  all(!is.na(w)) && all(w >= -1e-9) }

cat("\n=== 3. and where the standard astrophysical source classes land\n\n")
classes <- list("pion decay (1:2:0)"=c(1,2,0), "muon-damped (0:1:0)"=c(0,1,0),
                "neutron decay (1:0:0)"=c(1,0,0), "charm / equal (1:1:1)"=c(1,1,1))
cat("        source class              nu_e     nu_mu    nu_tau   reachable by decay?\n")
for (nm in names(classes)) {
  f <- classes[[nm]]/sum(classes[[nm]]); e <- as.vector(f %*% P)
  cat(sprintf("   %-24s %8.4f %8.4f %8.4f %16s\n", nm, e[1], e[2], e[3],
      if (inhull(e)) "YES" else "no"))
}
cat("\n  All four, and necessarily so: every astrophysical source composition is itself a\n")
cat("  point of the source simplex, and the decay component can reproduce any source\n")
cat("  composition by choosing the Yukawa column. Flavour at Earth cannot separate them\n")
cat("  and no sampling was needed to see it once the map is recognised as linear.\n")

cat("\n=== 4. flatly: the answer to the question I set is NO\n\n")
cat("  Every standard astrophysical flavour composition lies inside the region a\n")
cat("  free-Yukawa decay can produce, so an Earth-side flavour measurement cannot\n")
cat("  separate the decay component from an astrophysical population. The 1:1:2 channel\n")
cat("  ratio does not help, because it is flavour-blind: all three channels carry the\n")
cat("  same |y_alpha|^2 and the ratio comes from gauge and Higgs structure instead.\n")
cat("  So F-AG's 58-event directional requirement does NOT drop. Flavour is not a\n")
cat("  second discriminant and 5.2 is right to leave it free.\n")

cat("\n=== 5. the route, same turn, and it is a different observable\n\n")
cat("  What the 1:1:2 ratio DOES fix is the split between the hard line and everything\n")
cat("  else. Of four decays, one gives h nu and one Z nu - both contributing the hard\n")
cat("  line at M_1/2 - while two give W l, whose charged lepton never reaches Earth and\n")
cat("  whose products land in the soft continuum. So the LINE-TO-CONTINUUM normalisation\n")
cat("  is fixed by tree-level structure and is independent of the free Yukawa column.\n")
cat("  An astrophysical population has no reason to produce a fixed ratio of a line to\n")
cat("  its own continuum, so this IS a discriminant where flavour is not.\n")
cat("  NOT COMPUTED HERE, and it needs a tool I do not have: turning that into a number\n")
cat("  requires a decay-spectrum code to get the continuum shape from h, Z and W\n")
cat("  fragmentation. The structure is parameter-free; the ratio is not yet a number.\n")
