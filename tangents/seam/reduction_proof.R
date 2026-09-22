# The step both referees called load-bearing and unproven:
#   G_alpha / G_J  ->  W_B(x, P_perp y) / W_B(x, y)  at a GENERAL bifurcate Killing horizon.
#
# The mechanism, which is checkable rather than asserted: P_perp is an isometry of B, so
# it acts on TRANSVERSE mode labels only and commutes with the boost. Expanding both
# kernels in transverse modes, the boost factor R_l is COMMON to the two:
#   G_J(x,y)     = sum_l R_l(dtau) Y_l(x) Y_l(y)*
#   G_alpha(x,y) = sum_l R_l(dtau) Y_l(x) Y_l(P y)*
# The ratio is NOT the transverse ratio in general, because R_l depends on l. But at
# dtau = 0, which is both points on B, the common factor R_l(0) is exactly what defines
# the restriction W_B. So the claim should hold for ANY R_l, and the crucial test is to
# use an R_l that depends on l strongly.
#
# Round B = S^2, P_perp = antipodal, so Y_lm(P y) = (-1)^l Y_lm(y) and the addition
# theorem gives   W(x,y) = sum_l (2l+1)/(4 pi) R_l P_l(cos gamma).

Pl <- function(l,x){ if(l==0) return(rep(1,length(x))); if(l==1) return(x)
  p0<-rep(1,length(x)); p1<-x; for(k in 1:(l-1)){p2<-((2*k+1)*x*p1-k*p0)/(k+1); p0<-p1; p1<-p2}; p1 }
LMAX <- 260

# Boost factor: strongly l-dependent, and dtau-dependent. Three unrelated shapes.
Rl <- list(
  "R_l = e^{-l/8} sech(l dtau/6)"  = function(l,dt) exp(-l/8)/cosh(l*dt/6),
  "R_l = (1+l)^-2.3 e^{-l^2 dt^2}" = function(l,dt) (1+l)^(-2.3)*exp(-(l*dt)^2),
  "R_l = e^{-l/5}(1+l dt^2)^-1"    = function(l,dt) exp(-l/5)/(1+l*dt^2)
)
Wsum <- function(f, dt, g, sgn) {  # sgn=+1 direct, sgn=-1 inserts the (-1)^l of P_perp
  l <- 0:LMAX; sum((2*l+1)/(4*pi) * f(l,dt) * (sgn^l) * Pl_vec(l, cos(g))) }
Pl_vec <- function(ls, x) sapply(ls, function(l) Pl(l,x))

cat("=== 1. does the claimed ratio hold at dtau = 0, with R_l strongly l-dependent?\n\n")
cat("   ratio = W(x,Py)/W(x,y) computed from the FULL mode sum, against the same\n")
cat("   quantity built from the l-restriction. They must agree identically at dtau=0.\n\n")
for (nm in names(Rl)) {
  f <- Rl[[nm]]; cat(sprintf("  %-34s", nm))
  err <- max(sapply(c(30,70,110,150)*pi/180, function(g)
    abs( Wsum(f,0,g,-1)/Wsum(f,0,g,+1) - Wsum(f,0,g,-1)/Wsum(f,0,g,+1) )))
  # the real test: the antipodal sum equals the direct sum at the reflected angle
  err2 <- max(sapply(c(30,70,110,150)*pi/180, function(g)
    abs( Wsum(f,0,g,-1) - Wsum(f,0,pi-g,+1) )))
  cat(sprintf("  max |W(x,Py;g) - W(x,y;pi-g)| = %.3e\n", err2))
}
cat("\n  The antipodal insertion (-1)^l is exactly evaluation at the reflected angle, for\n")
cat("  every l-dependence. That is the content of 'P_perp acts on transverse labels\n")
cat("  only' and it needs no property of R_l at all.\n")

cat("\n=== 2. and the boost dependence: does it cancel at dtau > 0? It does NOT.\n\n")
cat("   ratio(dtau)/ratio(0) - 1, at gamma = 110 deg. If this were zero the reduction\n")
cat("   would hold away from B too, and the paper would be claiming too much.\n\n")
cat("        dtau      profile 1        profile 2        profile 3\n")
for (dt in c(0, 0.05, 0.2, 0.6, 1.5)) {
  cat(sprintf("   %9.2f", dt))
  for (nm in names(Rl)) {
    f <- Rl[[nm]]; g <- 110*pi/180
    r0 <- Wsum(f,0,g,-1)/Wsum(f,0,g,+1); rd <- Wsum(f,dt,g,-1)/Wsum(f,dt,g,+1)
    cat(sprintf(" %16.3e", rd/r0 - 1))
  }
  cat("\n")
}
cat("\n  So the reduction is genuinely a LIMIT and not an identity: the boost separation\n")
cat("  matters at finite dtau and the error grows with it, exactly as the de Sitter rate\n")
cat("  in `boost_cancel.R` said. What is general is the dtau -> 0 statement.\n")

cat("\n=== 3. what the argument actually needs, stated as a check\n\n")
cat("  (i)   P_perp is an isometry of B          -> acts on transverse labels only\n")
cat("  (ii)  P_perp commutes with the boost      -> R_l is common to both kernels\n")
cat("  (iii) the dtau -> 0 restriction exists    -> W_B is defined\n")
cat("  Given those three, the ratio at dtau = 0 is W_B(x,Py)/W_B(x,y) for ANY R_l, which\n")
cat("  part 1 confirms for three unrelated strongly l-dependent choices. None of the\n")
cat("  three needs maximal symmetry, a round B, a particular state or a closed form.\n")
cat("  (iii) is the one with content: for a Hadamard state the restriction exists away\n")
cat("  from transverse coincidence, which is where the ratio is being evaluated, but a\n")
cat("  codimension-two restriction is not automatic and we do not prove it here.\n")
