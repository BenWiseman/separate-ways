# kernel_ratio_horizon.R -- the fold kernel against the thermal contour's, App. A.6.
# Base R only. Units H = 1. Conformal scalar on dS4, Bunch-Davies, G ~ 1/(1 - Z).
#
# Static-patch embedding, two points at radii r1, r2, angle gamma at the centre,
# time separation dt:
#     Z(x,y)   =  sqrt((1-r1^2)(1-r2^2)) cosh(dt) + r1 r2 cos(gamma)
# with J: (X0,X1) -> (-X0,-X1) and P_perp: Xvec -> -Xvec acting on y,
#     Z_J      = -sqrt((1-r1^2)(1-r2^2)) cosh(dt) + r1 r2 cos(gamma)
#     Z_alpha  = -sqrt((1-r1^2)(1-r2^2)) cosh(dt) - r1 r2 cos(gamma)
# hence the invariant difference Z_alpha - Z_J = -2 r1 r2 cos(gamma), independent of dt.
#
# THE HEADLINE. At the horizon (r -> 1) the radial factor vanishes and
#     1 - Z_J     -> 1 - cos(gamma),      1 - Z_alpha -> 1 + cos(gamma),
# so each kernel is singular exactly where its own map makes the two points coincide:
# the thermal contour's at gamma = 0, the fold's at gamma = pi. Away from those,
#     G_alpha / G_J  ->  tan^2(gamma/2),
# independent of radius and time separation, with the two equal at gamma = pi/2.
#
# An earlier reading of this appendix claimed the fold kernel is simply regular at the
# horizon. That is true only on the gamma = 0 slice and is not the general statement;
# the tan^2 law is. Both are checked below so the slice is not mistaken for the law.
#
# Source of the expected values: calc/tangents/scales/two_point_off_axis.py,
# which_kernel_is_singular.py and regularity_at_general_angle.py (2026-09-18).
source("helpers.R")

rad  <- function(d) d*pi/180
sfac <- function(r1, r2, dt) sqrt((1-r1^2)*(1-r2^2))*cosh(dt)
Z_J  <- function(r1, r2, g, dt) -sfac(r1,r2,dt) + r1*r2*cos(g)
Z_a  <- function(r1, r2, g, dt) -sfac(r1,r2,dt) - r1*r2*cos(g)
ratio<- function(r1, r2, g, dt) (1 - Z_J(r1,r2,g,dt)) / (1 - Z_a(r1,r2,g,dt))   # G_alpha/G_J

cat(sprintf("invariant difference Z_a - Z_J = -2 r1 r2 cos(gamma), independent of dt\n"))
cat(sprintf("horizon ratio G_alpha/G_J -> tan^2(gamma/2); equal at gamma = pi/2\n\n"))

cat("checks:\n")
# the invariant difference, including unequal radii and nonzero dt
for (p in list(c(0.3,0.3,0,0), c(0.8,0.8,0,2.5), c(0.9,0.2,30,0), c(0.6,0.6,120,1))) {
  r1 <- p[1]; r2 <- p[2]; g <- rad(p[3]); dt <- p[4]
  report(sprintf("Z_a - Z_J = -2 r1 r2 cos(g) at r=(%.1f,%.1f) g=%.0f dt=%.1f", r1, r2, p[3], dt),
         expected = -2*r1*r2*cos(g), reproduced = Z_a(r1,r2,g,dt) - Z_J(r1,r2,g,dt),
         tol = 1e-12, mode = "abs")
}
# A.6's two quoted percentages, fold-normalised
for (r in c(0.3, 0.8)) {
  report(sprintf("A.6 fold-normalised discrepancy at r = %.1f, dt = 0", r),
         expected = r^2/(1-r^2),
         reproduced = abs((Z_a(r,r,0,0) - Z_J(r,r,0,0)) / (1 - Z_J(r,r,0,0))))
}
# the gamma = 0 slice
report("1 - Z_alpha = 2 identically on the gamma = 0, dt = 0 slice (r = 0.9)",
       expected = 2, reproduced = 1 - Z_a(0.9,0.9,0,0), tol = 1e-12, mode = "abs",
       note = "true on this slice only; see the tan^2 law below")
report("G_alpha/G_J = 1 - r^2 on that slice (r = 0.9)",
       expected = 1 - 0.81, reproduced = ratio(0.9,0.9,0,0))
# the two singular angles at the horizon
report("thermal kernel singular at gamma = 0: 1 - Z_J -> 0 as r -> 1",
       expected = 0, reproduced = 1 - Z_J(1-1e-9, 1-1e-9, 0, 0), tol = 1e-8, mode = "abs")
report("fold kernel singular at gamma = pi: 1 - Z_alpha -> 0 as r -> 1",
       expected = 0, reproduced = 1 - Z_a(1-1e-9, 1-1e-9, pi, 0), tol = 1e-8, mode = "abs",
       note = "so neither kernel is the better behaved overall")
# the law itself
for (d in c(20, 45, 60, 90, 120, 150)) {
  report(sprintf("horizon ratio = tan^2(gamma/2) at gamma = %d deg", d),
         expected = tan(rad(d)/2)^2, reproduced = ratio(1-1e-9, 1-1e-9, rad(d), 0), tol = 1e-5)
}
report("the two kernels agree exactly at gamma = pi/2", expected = 1,
       reproduced = ratio(1-1e-9, 1-1e-9, pi/2, 0), tol = 1e-5)
report("the law is independent of time separation (gamma = 60 deg, dt = 3)",
       expected = tan(rad(60)/2)^2, reproduced = ratio(1-1e-9, 1-1e-9, rad(60), 3), tol = 1e-5)
