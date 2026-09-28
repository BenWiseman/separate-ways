# Which of the two caustic models is a black hole? The answer is in the fold's differential,
# and A.19 used the wrong one.
#
# Two exactly solvable caustic geometries give opposite answers:
#   A.19's model, R_t x R_r x S^2(a):   image stress IDENTICALLY ZERO
#   the ESU,      R_t x S^3(a):         image stress POSITIVE, diverging as delta^-4
# Both have a caustic at the antipode. The difference has to be structural, and it is.
#
# A.19's fold is Theta(t, r, n) = (T - t, r, -n). It FIXES the radial coordinate and the
# radial direction: the pullback is (-1, +1, -1) on (time, radius, sphere). That +1 is what
# makes the p^2 terms cancel between the radial operator and the mixed derivative, which is
# the whole of A.19's zero.
#
# The real fold is J o P_perp with J: (U, V) -> (-U, -V). That reverses BOTH Kruskal null
# coordinates, so it is -Id on the entire two-dimensional Kruskal plane, and P_perp is -Id
# on the sphere. The real fold's differential is -Id on all four directions. It does not fix
# the radial direction. It looks like the ESU's fold, not like A.19's.
#
# If that is right, A.19's zero is an artefact of a fold that is not the fold, and the
# black-hole case belongs with the ESU's positive divergence.

cat("=== 1. what J does to the Kruskal directions, from the map itself ===\n")
J <- function(UV) -UV                      # (U,V) -> (-U,-V)
cat("   J is linear, so its differential is its matrix: diag(-1, -1) on (U, V).\n")
cat("   In (T, X) with T = (U+V)/2, X = (V-U)/2 that is also diag(-1, -1).\n")
for (v in list(c(1, 0), c(0, 1), c(1, 1), c(2, -3))) {
  cat(sprintf("     (U,V) = (%+.0f,%+.0f)  ->  (%+.0f,%+.0f)\n", v[1], v[2],
              J(v)[1], J(v)[2]))
}
cat("   So the fold reverses the radial direction. A.19's model does not.\n")

cat("\n=== 2. does r survive as a fixed coordinate? yes, and that is not the same thing ===\n")
M <- 1
UVofr <- function(r) (r / (2 * M) - 1) * exp(r / (2 * M))
cat("   UV is preserved by J, and r is a function of UV, so r is unchanged:\n")
for (r in c(0.5, 1.0, 1.5)) {
  cat(sprintf("     r = %.1f M   UV = %+8.5f   J(UV) = %+8.5f   same\n",
              r, UVofr(r), UVofr(r)))
}
cat("   A fixed COORDINATE is not a fixed DIRECTION. r is even under J because it is a\n")
cat("   function of the product UV, which is even; the radial vector field is odd. A.19\n")
cat("   read the first as the second and gave the radial pullback +1.\n")

cat("\n=== 3. the consequence, in A.19's own formula ===\n")
cat("   A.19: rho = (1/2) sum c_l (-1)^l [ f d_r^2 g + (f' + 2f/r) d_r g + s_r f d_r d_r' g ]\n")
cat("   with s_r the radial pullback. In WKB, d_r^2 g -> (-p^2 + ...) and d_r d_r' g ->\n")
cat("   (+p^2 + ...), so the p^2 terms cancel at s_r = +1 and DOUBLE at s_r = -1.\n\n")
p2 <- 1                                     # unit p^2, to show the structure
for (s_r in c(+1, -1)) {
  cat(sprintf("     s_r = %+d :  f(-p^2) + s_r f(+p^2) = %+.0f x f p^2   %s\n",
              s_r, -1 + s_r, ifelse(s_r > 0, "cancels, A.19's zero",
                                    "doubles, no cancellation")))
}

cat("\n=== 4. the cross-check: does the ESU fold have s_r = -1? ===\n")
cat("   On S^3 the antipodal map is x -> -x in the embedding R^4, so its differential is\n")
cat("   -Id on the whole tangent space, every direction included. The ESU therefore has\n")
cat("   s = (-1, -1, -1, -1), and it gave a positive divergence. The Kruskal fold has the\n")
cat("   same signature. A.19's model is the odd one out, with a +1 nobody chose on\n")
cat("   physical grounds: it was chosen because it made the transverse radius constant.\n")

cat("\n=== 5. so what does this say ===\n")
cat("   A.19's identically zero image stress is a property of a fold that fixes the radial\n")
cat("   direction. The fold of this paper does not fix it: J reverses both null coordinates\n")
cat("   and therefore reverses the radial direction, while leaving the radial COORDINATE\n")
cat("   unchanged because r depends on UV, which is even.\n")
cat("   The black-hole case therefore sits with the ESU, whose exact answer is a positive\n")
cat("   divergence going as the inverse fourth power of the distance from contact, and not\n")
cat("   with A.19's zero.\n")
cat("   That is a claim about a signature, so it is worth saying what would break it: if\n")
cat("   the correct assembly inside the horizon uses the radial pullback of the COORDINATE\n")
cat("   rather than of the VECTOR, the +1 is right and A.19 stands. Settling that needs the\n")
cat("   point-split done covariantly rather than in components, which is the next step.\n")

cat("\n=== 6. TESTING the claim in section 1-5, which is how it died ===\n")
cat("   The pullback of a VECTOR is not read off the map's matrix alone: the basis vector\n")
cat("   at the IMAGE point must be compared with, and it is built from the image point's\n")
cat("   own coordinates. Do it properly.\n\n")
cat("   Inside, constant t is the ray (U,V) -> (lam U, lam V), so d_r is proportional to\n")
cat("   U d_U + V d_V. Constant r is the hyperbola, and d_t is proportional to\n")
cat("   -U d_U + V d_V.\n\n")
dr <- function(U, V) c(U, V)
dt <- function(U, V) c(-U, V)
dphi <- function(v) -v                     # the differential of J is -Id on (U,V)

cat("      vector      at x=(U,V)      dJ(at x)        at Theta x        pullback\n")
for (uv in list(c(0.6, 0.9), c(1.3, 0.4))) {
  U <- uv[1]; V <- uv[2]
  for (nm in c("d_r", "d_t")) {
    vx  <- if (nm == "d_r") dr(U, V) else dt(U, V)
    vim <- if (nm == "d_r") dr(-U, -V) else dt(-U, -V)
    pb  <- dphi(vx) / vim
    cat(sprintf("      %-4s  (%+.2f,%+.2f)  (%+.2f,%+.2f)  (%+.2f,%+.2f)   %+.0f\n",
                nm, vx[1], vx[2], dphi(vx)[1], dphi(vx)[2], vim[1], vim[2], pb[1]))
  }
}
cat("\n   BOTH pullbacks are +1. The map's -Id and the image basis vector's own reversal\n")
cat("   cancel. So the radial pullback is +1 after all, A.19's cancellation stands, and\n")
cat("   sections 1-5 above are WRONG. Recorded rather than deleted.\n")

cat("\n=== 7. and the failure says something the success would not have ===\n")
cat("   Every coordinate pullback comes out +1, including the time one. Yet the fold\n")
cat("   demonstrably reverses time orientation: J carries the right wedge to the left,\n")
cat("   where the Killing field d_t points to the past. The reversal is in the CAUSAL\n")
cat("   CHARACTER of the Killing field between wedges, not in any component sign.\n")
cat("   Component bookkeeping cannot see that, which means every Group F assembly that\n")
cat("   assigned pullback signs by inspection, A.19's included, is resting on a method\n")
cat("   that is blind to the one property the fold is defined by. The point-split has to\n")
cat("   be done covariantly.\n")
