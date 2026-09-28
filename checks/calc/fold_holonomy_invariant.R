# The chart-independent version of the fold's pullback, and whether it predicts the answer.
#
# GR53 and GR54 between them showed the component pullback is meaningless on its own: the
# SAME fold reads diag(-1,-1,-1,-1) in flat Cartesian and +1 on the radial direction in
# spherical or Kruskal, because the basis is position-dependent. So the sign is chart noise.
#
# The invariant is P = (parallel transport from Theta x back to x) o dTheta. That maps T_x
# to ITSELF, so its eigenvalues are real invariants of the fold and the geodesic joining the
# pair. Compute it for the three cases whose answers are known and see whether it sorts them.
#
#   flat / Rindler (GR54):  exact, POSITIVE and divergent
#   Einstein static (GR52): exact, POSITIVE, delta^-4
#   A.19 constant radius:   exact, ZERO

cat("=== 1. parallel transport to the antipode on a sphere, computed not assumed ===\n")
# transport a vector along a great circle on the unit 2-sphere, by integrating the equation
transport <- function(v0, steps = 20000) {
  # great circle from north pole in the x-direction; transport v0 (in the embedding)
  v <- v0
  for (i in 1:steps) {
    s <- pi * (i - 0.5) / steps; ds <- pi / steps
    n <- c(sin(s), 0, cos(s))                    # position on the circle
    tdot <- c(cos(s), 0, -sin(s))                # unit tangent
    # parallel transport: dv/ds = -(v . tdot) n   (project out the normal component)
    v <- v + ds * (-sum(v * tdot) * n)
    v <- v - sum(v * n) * n                      # keep it tangent
  }
  v
}
along <- transport(c(1, 0, 0))     # tangent to the transport circle at the north pole
perp  <- transport(c(0, 1, 0))     # perpendicular to it
cat(sprintf("   tangent  (1,0,0) transported to the antipode: (%+.4f, %+.4f, %+.4f)\n",
            along[1], along[2], along[3]))
cat(sprintf("   perpend. (0,1,0) transported to the antipode: (%+.4f, %+.4f, %+.4f)\n",
            perp[1], perp[2], perp[3]))
cat("   so transport REVERSES the along-circle direction and PRESERVES the perpendicular.\n")
stopifnot(abs(along[1] + 1) < 1e-3, abs(perp[2] - 1) < 1e-3)

cat("\n=== 2. compose with the antipodal map, which is -Id in the embedding ===\n")
cat("   P = transport o dTheta, per direction:\n")
cat(sprintf("     along-circle :  dTheta -1, transport reverses  ->  P = %+d\n", +1))
cat(sprintf("     perpendicular:  dTheta -1, transport preserves ->  P = %+d\n", -1))

cat("\n=== 3. P for the three known cases, and their answers ===\n")
cases <- list(
  list(nm = "flat / Rindler", P = c(-1, -1, -1, -1), ans = "POSITIVE"),
  list(nm = "Einstein static", P = c(-1, +1, -1, -1), ans = "POSITIVE"),
  list(nm = "A.19 const. radius", P = c(-1, +1, +1, -1), ans = "ZERO"))
cat("      case                 P eigenvalues        +1 count   trace   answer\n")
for (c0 in cases) {
  cat(sprintf("   %-20s (%s)      %d        %+3d    %s\n", c0$nm,
              paste(sprintf("%+d", c0$P), collapse = ","),
              sum(c0$P > 0), sum(c0$P), c0$ans))
}
cat("\n   In flat Cartesian the connecting geodesic is the straight line through the origin\n")
cat("   and transport is trivial, so P = dTheta = -Id, no +1 at all.\n")
cat("   In the ESU the one +1 is the along-circle direction on S^3.\n")
cat("   In A.19's model there are two: the along-circle AND the radial direction, which\n")
cat("   that fold fixes outright and transport does not touch because r is flat and the\n")
cat("   geodesic has no radial component.\n")

cat("\n=== 4. the pattern, and how much weight it can carry ===\n")
cat("   Zero or one +1 gives a positive divergence; two gives the cancellation. The radial\n")
cat("   +1 is the extra one, and A.19's own algebra says why it matters: it is exactly the\n")
cat("   sign that makes the p^2 terms cancel between the radial operator and the mixed\n")
cat("   derivative. The along-circle +1 sits in the transverse sector where the parity is\n")
cat("   already doing the work, and does not produce a cancellation.\n")
cat("   Three points is a pattern and not a theorem. What makes it more than numerology is\n")
cat("   that the mechanism is identified rather than fitted.\n")

cat("\n=== 5. so what does a black hole have? ===\n")
cat("   The question is whether the radial direction carries +1 or -1 under P, which needs\n")
cat("   transport along the contact curve, the E = 0 null geodesic through the bifurcation\n")
cat("   surface. Not computed here. But the flat case settles one thing already: A.19's\n")
cat("   model fixes r because r is a FLAT product direction there and the connecting\n")
cat("   geodesic never moves in it. At a black hole the contact curve crosses the whole\n")
cat("   interior in r, so transport along it cannot leave the radial direction alone. The\n")
cat("   presumption is therefore against A.19's second +1 surviving, which puts a black\n")
cat("   hole with the positives. Presumption, not result.\n")
