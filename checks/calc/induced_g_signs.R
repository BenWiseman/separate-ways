# Checking the obstruction M3 put in the way of the induced-gravity route.
#
# M3 walked that route and recorded an impasse with two load-bearing claims. The first is that the
# induced 1/G is dominated by the UV cutoff, so matching the observed Planck mass needs either
# Lambda_UV ~ M_Pl, which is circular, or a species count the fold does not supply. The second is
# that the fold's dark matter is fermionic and "the fermion contribution to 1/G has the wrong sign
# in any regularization", so the fold actively hurts induced gravity.
#
# The second claim is the one worth checking, because a FALSE obstruction left in a trail is worse
# than no trail: it would stop a later attempt for no reason. It is wrong.

cat("=== 1. the heat-kernel coefficients, written out rather than quoted ===\n")
cat("   Tr exp(-sA) = (4 pi s)^{-2} int sqrt(g) [ a_0 + a_1 s + ... ] for A = -box + X, with\n")
cat("        a_0 = tr 1,        a_1 = tr( R/6 - X ).\n")
cat("   The one-loop action is +(1/2) Tr ln A for a real boson and -(1/2) Tr ln A for a Dirac\n")
cat("   fermion, the sign being the Grassmann determinant. The R piece of a_1, times that sign,\n")
cat("   is what induces the Einstein-Hilbert term.\n\n")
cat("      field                 tr 1   X                    a_1's R piece   loop sign   product\n")
rows <- list(
  list("real scalar, xi",      1, "m^2 + xi R",     NA,         +1),
  list("real scalar, xi = 0",  1, "m^2",            1/6,        +1),
  list("real scalar, xi = 1/6",1, "m^2 + R/6",      0,          +1),
  list("real scalar, xi = 1/3",1, "m^2 + R/3",      1/6 - 1/3,  +1),
  list("Dirac fermion",        4, "m^2 + R/4",      4/6 - 4/4,  -1))
for (r in rows) {
  if (is.na(r[[4]])) {
    cat(sprintf("      %-21s %4d   %-20s   %-13s   %+d          %s\n",
                r[[1]], r[[2]], r[[3]], "(1/6 - xi)", r[[5]], "(1/6 - xi)"))
  } else {
    cat(sprintf("      %-21s %4d   %-20s   %+13.4f   %+d       %+10.4f\n",
                r[[1]], r[[2]], r[[3]], r[[4]], r[[5]], r[[4]] * r[[5]]))
  }
}
cat("\n   The Dirac row is the point. Its a_1 R piece is 4(1/6) - 4(1/4) = -1/3, and the loop\n")
cat("   sign is -1, so the product is +1/3: the SAME sign as a minimally coupled scalar, and\n")
cat("   twice its size. Per degree of freedom that is +1/12 against the scalar's +1/6.\n")
dirac <- (4/6 - 4/4) * (-1)
stopifnot(dirac > 0)
cat(sprintf("      Dirac product %+0.4f, minimal scalar product %+0.4f, both positive.\n", dirac, 1/6))

cat("\n=== 2. so the fermion claim is wrong, and the reason it is wrong matters ===\n")
cat("   Fermions do not push the induced 1/G the wrong way. What can push it the wrong way is a\n")
cat("   SCALAR with xi > 1/6, since its contribution carries (1/6 - xi):\n")
for (xi in c(0, 1/6, 1/4, 1/3, 1)) {
  v <- 1/6 - xi
  cat(sprintf("      xi = %5.3f:  contribution proportional to %+7.4f   %s\n", xi, v,
              ifelse(v > 0, "helps", ifelse(v == 0, "contributes nothing", "hurts"))))
}
cat("   That is also the same factor image_stress_conformal.R found controlling the fold's\n")
cat("   image stress, which is a coincidence worth noticing and not more than that.\n")
cat("   The standard species counting adds bosons and fermions with the same sign, which is why\n")
cat("   the species bound is a bound and not a cancellation.\n")

cat("\n=== 3. the first obstruction, which does survive ===\n")
MPl <- 1.220890e19      # GeV
M1  <- 491.6e6          # GeV, the fold's dark-matter ceiling: 491.6 PeV
cat(sprintf("      M_Pl = %.4e GeV, the fold's highest unambiguous scale M_1 = %.4e GeV\n", MPl, M1))
cat(sprintf("      (M_Pl/M_1)^2 = %.3e\n", (MPl / M1)^2))
cat(sprintf("      with the 6 pi of the naive count: N ~ %.2e species needed\n", 6 * pi * (MPl / M1)^2))
cat("   Twenty-one to twenty-two orders of magnitude of species, which the fold does not supply\n")
cat("   and which no reasonable completion supplies either. The route dies here, on the count,\n")
cat("   and not on a sign.\n")
stopifnot(6 * pi * (MPl / M1)^2 > 1e21)

cat("\n=== 4. the plant ===\n")
cat("   The sign test must fail if the loop sign is mishandled, or it is checking nothing.\n")
for (s in c(-1, +1)) {
  v <- (4/6 - 4/4) * s
  cat(sprintf("      Dirac with loop sign %+d:  product %+0.4f  ->  %s\n", s, v,
              ifelse(v > 0, "same sign as a boson", "the wrong sign, which is the claim being tested")))
}
cat("   Only the +1 loop sign, which is the wrong one for a fermion, reproduces M3's claim. So\n")
cat("   the claim is what you get by forgetting the Grassmann determinant.\n")
cat("\n=== 5. what this does to the trail ===\n")
cat("   One of the two load-bearing obstructions is removed and the other stands. The route is\n")
cat("   still dead, on the species count alone. That matters for the next attempt: anyone\n")
cat("   revisiting induced gravity here should not believe that the fold's fermionic dark matter\n")
cat("   is an obstacle, because it is not, and should go straight at the count.\n")
