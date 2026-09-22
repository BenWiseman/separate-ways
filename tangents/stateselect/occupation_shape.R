# F-AY asked "what fixes the LEVEL" as though the occupation were a single scale-independent
# number. It is not: it is a function n(k). Ben's standing point. So ask instead which PARTS
# of n(k) are fixed and by what, which is a scale-dependent question with a scale-dependent
# answer.
#
# Hadamard regularity is a SHORT-DISTANCE condition, so it constrains the LARGE-k tail of
# n(k): a state is Hadamard only if its two-point function differs from the vacuum by a
# smooth function, which fails if n_k carries a power-law tail. Test that directly by
# computing the difference kernel and looking at its smoothness.

# A first version summed n_k cos(k s)/k and watched it as s -> 0. That test cannot
# discriminate: the 1/k weight and the truncation make every tail converge, and all five
# rows came back "not diverging", so the prose was asserting what the table never showed.
#
# The decisive test is sharper and needs no truncation games. The m-th derivative of
# sum_k n_k cos(k s) at s = 0 is +/- sum_k n_k k^m. The difference from the vacuum is
# SMOOTH at coincidence iff that converges for EVERY m, which happens iff n_k decays
# faster than any power. That is exactly the Hadamard requirement.

kmax <- 200000; k <- 1:kmax
tails <- list(
  "power  n ~ k^-2"       = function(k) k^-2,
  "power  n ~ k^-4"       = function(k) k^-4,
  "power  n ~ k^-8"       = function(k) k^-8,
  "exponential e^{-k/20}" = function(k) exp(-k/20),
  "gaussian e^{-(k/20)^2}"= function(k) exp(-(k/20)^2)
)
mom <- function(f, m, K) sum(f(1:K)*(1:K)^m)
# A second version compared the SIZE of the moment at fixed m, which discriminates nothing:
# moments grow factorially even for a convergent tail, so the exponential looked as bad as
# the power law. The real question is CONVERGENCE, i.e. whether extending the cutoff changes
# the answer.

cat("=== 1. which occupations leave the short-distance structure alone?\n\n")
cat("   sum_k n_k k^m must CONVERGE for every m. Test by doubling the cutoff: a convergent\n")
cat("   sum stops changing, a divergent one keeps growing.\n\n")
cat("        tail                    m     K=50k        K=100k     ratio   converges?\n")
for (nm in names(tails)) {
  f <- tails[[nm]]
  for (m in c(6, 12)) {
    a <- mom(f,m,50000); b <- mom(f,m,100000)
    cat(sprintf("   %-22s %5d %11.2e %13.2e %7.1f %10s\n", if (m==6) nm else "", m, a, b, b/a,
        if (abs(b/a - 1) < 1e-3) "YES" else "no"))
  }
}
cat("\n  Read the ratios, not the flags. A power law converges while m is below its index\n")
cat("  and diverges beyond: k^-8 is fine at m = 6 (ratio 1.0) and gone at m = 12 (ratio\n")
cat("  32). So EVERY power law fails at some m, which is what matters, since Hadamard\n")
cat("  needs every m. The exponential and gaussian hold at 1.0 throughout.\n")
cat("  The exponential and gaussian\n")
cat("  stop changing entirely. Only faster-than-any-power decay leaves the short-distance\n")
cat("  structure alone, which is what Hadamard regularity demands - and it is a condition\n")
cat("  on the LARGE-k tail of the occupation, not on its overall size.\n")

cat("\n=== 2. so the occupation is fixed in PARTS, by different things\n\n")
cat("        feature of n(k)            fixed by                       status\n")
rows <- list(
 c("sheet symmetry n = m","Theta-invariance","FORCED (F-AY)"),
 c("large-k tail","Hadamard regularity","FORCED, must decay fast"),
 c("intermediate/IR shape","nothing here","IMPORTED from 5.1"))
for (r in rows) cat(sprintf("   %-26s %-30s %s\n", r[1], r[2], r[3]))

cat("\n=== 3. the answer is SCALE-DEPENDENT, which is the point\n\n")
cat("  Asking what fixes 'the level' has no answer because there is no single level. The\n")
cat("  algebra fixes the state at SHORT distance, through the admissibility condition it\n")
cat("  already uses on the seam, and fixes the sheet symmetry at every scale through the\n")
cat("  involution. What 5.1 imports is the INFRARED profile, and only that.\n")
cat("  So the live falsifier is narrower again: not 'the algebra fails on states' and not\n")
cat("  'it half succeeds', but 'it fixes the ultraviolet and the symmetry, and the\n")
cat("  infrared is imported'. That is a statement about which scales the pairing reaches.\n")

cat("\n=== 4. and it says where to look next, which is not another uniqueness argument\n\n")
cat("  The relic abundance of 5.1 is set by the occupation at the production scale, which\n")
cat("  is neither the deep ultraviolet nor the far infrared. Whether the algebra reaches\n")
cat("  THAT scale is a question about the intermediate regime, and the tools that work\n")
cat("  there are renormalisation-group rather than algebraic: what matters is how the\n")
cat("  imported infrared profile RUNS up to the production scale, not what it is at the\n")
cat("  bottom. We do not compute that here.\n")
