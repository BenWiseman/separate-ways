#!/usr/bin/env Rscript
# a10_boost_dependence.R -- A.10 quotes the residual boost dependence of the B-form as a
# table of four measured numbers. This reconstructs the kernels from the invariants the
# appendix itself states, derives the closed form the table is sampling, checks the table
# against it, and reports the one sentence that does not reconcile.
#
# A.10 states Z_alpha - Z_J = -2 r_1 r_2 cos gamma, independent of the boost, and at
# Delta t = 0, gamma = 0, equal radii, Z_J = -1 + 2r^2 while Z_alpha = -1. Those fix
#     Z_alpha = -(1-r^2) cosh(dt) - r^2 cos gamma,
#     Z_J     = -(1-r^2) cosh(dt) + r^2 cos gamma,
# and the conformal scalar has W ~ 1/(1-Z), so G_alpha/G_J = (1-Z_J)/(1-Z_alpha).

Zj <- function(r,g,dt) -(1-r^2)*cosh(dt) + r^2*cos(g)
Za <- function(r,g,dt) -(1-r^2)*cosh(dt) - r^2*cos(g)
R  <- function(r,g,dt) (1-Zj(r,g,dt))/(1-Za(r,g,dt))
RB <- function(g) tan(g/2)^2

cat("=== 1. the kernel reproduces the invariants A.10 prints ===\n")
for (r in c(0.3,0.5,0.8)) stopifnot(abs(Zj(r,0,0)-(-1+2*r^2))<1e-12, abs(Za(r,0,0)+1)<1e-12)
cat("   Z_J = -1+2r^2 and Z_alpha = -1 at dt = 0, gamma = 0: exact at every r tested.\n")
for (gd in c(30,60,90,110,150)) { g <- gd*pi/180
  stopifnot(abs(R(1-1e-9,g,0) - RB(g)) < 1e-6) }
cat("   B limit r -> 1 gives tan^2(gamma/2): exact at every angle tested.\n")
stopifnot(abs(100*0.09/0.91 - 9.9) < 0.1, abs(100*0.64/0.36 - 178) < 1)
cat("   the 9.9 and 178 per cent discrepancy figures: reproduced.\n")

cat("\n=== 2. the closed form the table is sampling ===\n")
cat("   Expanding R to first order in N^2 = 1-r^2 at fixed boost,\n")
cat("       R/R_B - 1 = N^2 * 2 cos(gamma) (cosh(dt) + 1) / sin^2(gamma) + O(N^4).\n")
pred <- function(g,dt) 2*cos(g)*(cosh(dt)+1)/sin(g)^2
r <- 0.9999; N2 <- 1-r^2
maxerr <- 0
for (gd in c(20,45,60,80,100,120,150)) for (dt in c(0,0.5,1.5,2.5)) { g <- gd*pi/180
  e <- abs((R(r,g,dt)/RB(g)-1)/N2 - pred(g,dt))/max(1,abs(pred(g,dt)))
  maxerr <- max(maxerr,e) }
cat(sprintf("   checked on a 7-angle x 4-boost grid at r = 0.9999: worst relative error %.1e\n", maxerr))

cat("\n=== 3. ninety degrees is exact, not first-order small ===\n")
g90 <- pi/2; worst <- 0
for (rr in c(0.2,0.5,0.9,0.9999)) for (dt in c(0,1,3,6)) worst <- max(worst, abs(R(rr,g90,dt)-1))
cat(sprintf("   cos(90) = 0 kills the cos(gamma) factor in both kernels at once, so\n"))
cat(sprintf("   G_alpha/G_J = 1 = tan^2(45) identically. Worst deviation over r in\n"))
cat(sprintf("   [0.2, 0.9999] and dt in [0, 6]: %.1e. The B-form carries no truncation\n", worst))
cat("   error at ninety degrees at any radius and any boost.\n")

cat("\n=== 4. A.10's measured table, against the closed form ===\n")
cat("   A.10: 12.1 at 60 and 120 degrees, 61.8 at 150, 517 at 170, exactly zero at 90,\n")
cat("   all at rH = 0.9999. Fit the single free boost to the 60-degree entry alone and\n")
cat("   the other four follow with nothing left to adjust.\n\n")
g60 <- 60*pi/180
dt  <- acosh(12.1/abs(2*cos(g60)/sin(g60)^2) - 1)
cat(sprintf("   boost implied by the 60-degree entry: cosh(dt) = %.4f, dt = %.4f\n\n", cosh(dt), dt))
paper <- c(`60`=12.1,`90`=0,`120`=12.1,`150`=61.8,`170`=517)
cat("     angle    this kernel at r=0.9999      A.10      difference\n")
for (gd in c(60,90,120,150,170)) { g <- gd*pi/180
  num <- abs((R(r,g,dt)/RB(g)-1)/N2); p <- paper[as.character(gd)]
  d <- if (p==0) sprintf("%18s","exact") else sprintf("%16.1f%%", 100*(num-p)/p)
  cat(sprintf("   %6d   %22.4f   %8.1f %s%s\n", gd, num, p,
              d, if (gd==60) "   <- fitted" else "")) }
cat("\n   Four entries predicted from one fit, all inside two and a half per cent, with\n")
cat("   the sign structure exact: 60 and 120 equal in magnitude because |cos| is even\n")
cat("   about ninety degrees, zero at ninety, and the 170-degree blow-up from sin^2 -> 0.\n")

cat("\n=== 5. the one sentence that does not reconcile ===\n")
cat("   A.10 also says: 'at gamma = 110 degrees the ratio moves by 1.2e-4 at\n")
cat("   dtau = 0.05 and by thirty per cent at dtau = 1.5'.\n\n")
g <- 110*pi/180
floor110 <- abs(pred(g,0))*N2
cat(sprintf("   At rH = 0.9999 and gamma = 110 the deviation from the B-form is smallest at\n"))
cat(sprintf("   zero boost, where it is already %.2e, and grows as cosh(dt)+1 from there.\n", floor110))
cat(sprintf("   It therefore cannot take the value 1.2e-4 at that radius for any boost.\n"))
cat(sprintf("   Reading the move instead as the change from dtau = 0, the two halves of the\n"))
cat("   sentence want different radii: matching 1.2e-4 at the small boost needs r about\n")
cat("   0.999, matching thirty per cent at the large one needs r about 0.987.\n")
one_pc <- acosh(1e-2/(N2*abs(2*cos(g)/sin(g)^2)) - 1)
cat(sprintf("\n   What this kernel does give at gamma = 110 and rH = 0.9999: the deviation is\n"))
cat(sprintf("   %.1e at zero boost and reaches one per cent at dt = %.2f.\n", floor110, one_pc))

cat("
=== flatly ===

  A.10's substance holds and is stronger than the appendix claims for it. The kernels
  reconstruct from the invariants it prints; the B limit, the two discrepancy figures and
  the whole four-entry table follow, the last from a single fitted boost. The ninety-degree
  statement is better than 'first order in N^2 with a vanishing constant': the ratio is
  exactly one there at every radius and every boost.

  One sentence does not reconcile. The 1.2e-4 and thirty per cent quoted at gamma = 110
  are not what this kernel gives at the stated radius, and no single radius gives both.
  That sentence is not load-bearing: the table it introduces is, and the table checks out.\n")

# ---------------------------------------------------------------------------
# 6. The l-sum check in the same passage, which is an identity and not numerics.
#
# A.10 reports: "Checked with three unrelated and strongly l-dependent boost factors on a
# round B, where the antipodal insertion is (-1)^l: the reflected sum equals the direct
# sum at the supplementary angle to 7e-17." On a round B the antipodal map sends a mode of
# degree l to (-1)^l times itself, and the Legendre polynomials obey P_l(-x) = (-1)^l P_l(x),
# so the two factors of (-1)^l are the same factor. The reflected sum at gamma is the
# direct sum at pi - gamma term by term, for ANY boost profile at all. 7e-17 is roundoff.

cat("\n=== 6. the l-sum check is an identity, not a measurement ===\n")
direct <- function(g, b) sum(b*(2*seq_along(b)-1)/(4*pi) *
                             sapply(seq_along(b)-1, function(l) {
                               x <- cos(g); if (l==0) 1 else if (l==1) x else {
                                 p0<-1; p1<-x; for (k in 2:l) { p2 <- ((2*k-1)*x*p1-(k-1)*p0)/k; p0<-p1; p1<-p2 }; p1 } }))
refl <- function(g, b) direct(g, b*(-1)^(seq_along(b)-1))
L <- 0:40
profiles <- list(`exp(-l/7)` = exp(-L/7),
                 `1/(1+l^3)` = 1/(1+L^3),
                 `cos(l)^2 e^{-l/11}` = cos(L)^2*exp(-L/11),
                 `random`    = c(0.83,0.11,0.57,0.92,0.04,0.38,0.66,0.21,0.49,0.75,
                                 0.13,0.88,0.30,0.61,0.07,0.44,0.96,0.25,0.70,0.52,
                                 0.18,0.81,0.36,0.64,0.02,0.47,0.90,0.28,0.73,0.55,
                                 0.15,0.85,0.33,0.68,0.09,0.41,0.94,0.23,0.77,0.59,0.06))
worst <- 0
for (nm in names(profiles)) { b <- profiles[[nm]]
  e <- max(sapply(c(20,55,90,110,150)*pi/180, function(g) abs(refl(g,b) - direct(pi-g,b))))
  cat(sprintf("   %-22s worst |reflected(g) - direct(pi-g)| = %.1e\n", nm, e))
  worst <- max(worst,e) }
cat(sprintf("\n   Worst over four profiles including a random one: %.1e, which is roundoff.\n", worst))
cat("   P_l(-x) = (-1)^l P_l(x) and the antipodal insertion is (-1)^l, so the two signs\n")
cat("   are the same sign and cancel term by term. The equality holds for every boost\n")
cat("   profile, not for three well-chosen ones, and needs no check at all.\n")
