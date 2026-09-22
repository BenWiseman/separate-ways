# ratio(135 deg) * ratio(45 deg) = 2.8559 * 0.3502 = 1.00013... Check whether that
# is exact, general, and what it implies.
#
# The claim to test:  R(pi - gamma) * R(gamma) = 1  identically,
# where R(gamma) = F(Z_alpha(gamma)) / F(Z_J(gamma)).
# If true it is profile-free and radius-free, because under gamma -> pi - gamma
# the two invariants simply swap: Z_J(pi-g) = Z_alpha(g).

Zj <- function(r,g) -(1-r^2) + r^2*cos(g)
Za <- function(r,g) -(1-r^2) - r^2*cos(g)
prof <- list("conformal"=function(Z) 1/(1-Z), "heavy"=function(Z)(1-Z)^(-1.7),
  "exponential"=function(Z) exp(2*Z), "screened"=function(Z) exp(-2*(1-Z))/(1-Z),
  "arctan"=function(Z) atan(4*Z)+2, "sqrt"=function(Z) sqrt(2+Z),
  "5th order"=function(Z)(3+Z)/(1-Z)^5, "logistic"=function(Z) 1/(1+exp(-6*Z)),
  "gaussian"=function(Z) exp(-(1-Z)^2/8), "two poles"=function(Z) 1/(1-Z)+0.3/(2-Z))

cat("=== 1. the swap identity Z_J(pi-g) = Z_alpha(g), exactly\n\n")
e <- max(abs(sapply(c(0.1,0.5,0.9,0.99), function(r)
       sapply(seq(0.01,pi-0.01,length.out=500), function(g) Zj(r,pi-g)-Za(r,g)))))
cat(sprintf("  max |Z_J(pi-g) - Z_alpha(g)| over 4 radii x 500 angles = %.3e\n", e))
cat("  So reflecting the angle about ninety degrees SWAPS the two kernels. That is\n")
cat("  the whole content, and it is geometry: the antipode of the antipode is the\n")
cat("  point itself.\n")

cat("\n=== 2. therefore R(pi-g) R(g) = 1. max |R(pi-g)R(g) - 1| per profile\n\n")
for (nm in names(prof)) {
  F <- prof[[nm]]
  v <- max(abs(sapply(c(0.1,0.5,0.9,0.99), function(r)
         sapply(seq(0.05,pi-0.05,length.out=400), function(g)
           (F(Za(r,pi-g))/F(Zj(r,pi-g)))*(F(Za(r,g))/F(Zj(r,g))) - 1))))
  cat(sprintf("  %-14s %.3e\n", nm, v))
}
cat("\n  Machine precision, every profile, every radius, every angle. This is an exact\n")
cat("  identity of the fold and not a numerical coincidence. Equivalently:\n")
cat("    log R is an ODD function about gamma = 90 degrees.\n")
cat("  The ninety-degree equality is the fixed point of that oddness, so it is a\n")
cat("  corollary rather than a separate fact.\n")

cat("\n=== 3. the multipole consequence. Legendre content of log R at rH = 0.9\n")
cat("   (coefficients a_l = (2l+1)/2 * integral log R(g) P_l(cos g) sin g dg)\n\n")
Pl <- function(l,x){ if(l==0) return(rep(1,length(x))); if(l==1) return(x)
  p0<-rep(1,length(x)); p1<-x; for(k in 1:(l-1)){p2<-((2*k+1)*x*p1-k*p0)/(k+1); p0<-p1; p1<-p2}; p1 }
for (nm in c("conformal","screened","arctan","gaussian")) {
  F <- prof[[nm]]; r <- 0.9
  a <- sapply(0:7, function(l) integrate(function(g){
        x<-cos(g); log(F(Za(r,g))/F(Zj(r,g)))*Pl(l,x)*sin(g) }, 0.001, pi-0.001,
        subdivisions=2000, rel.tol=1e-10)$value*(2*l+1)/2)
  cat(sprintf("  %-10s", nm)); cat(sprintf(" l=%d:%9.2e", 0:7, a)); cat("\n")
}
cat("\n  Every EVEN multipole vanishes to integration precision; only ODD l survive.\n")
cat("  That is forced by the oddness in part 2, and it sits exactly opposite A.11:\n")
cat("  the Kerr latitude weight A(theta) has only EVEN Legendre content and gives\n")
cat("  the l -> l+/-2 rule, while the kernel ratio itself carries only ODD content.\n")
cat("  Two different objects, two opposite parities, one geometry. Neither was put\n")
cat("  in by hand.\n")
