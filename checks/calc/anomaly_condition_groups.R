# The Boyle-Finn-Turok anomaly condition, solved once and run over candidate gauge groups.
#
# WHAT THIS COMPUTES. BFT (arXiv:2110.06258) add n'_0 dimension-zero scalars and write the two
# Weyl anomaly coefficients as
#     a ~ n_0 + (11/2) n_half + 62 n_1 - 28 n'_0
#     c ~ n_0 +      3 n_half + 12 n_1 -  8 n'_0
# with n_1 = dim(G). Their eq. (17) solves these with the vacuum-energy condition, eq. (13):
# n_half = 4 n_1, n'_0 = 3 n_1, n_0 = 0, and they note the Standard Model group then needs three
# generations. This eliminates n'_0 between the two rows, which gives one condition in the
# fermion and scalar counts, and runs it over the groups they do not run it over.
#
# WHAT IT DOES NOT SETTLE. The two coefficient rows are an INPUT. They were checked against
# eqs. (16a) and (16b) of arXiv:2110.06258v2 on 2026-10-01 and match. The condition is theirs
# and rests on their dimension-zero construction, which this does not test.
#
# Base R only. No package is loaded.

cat("\n== BFT's anomaly condition, solved and run over groups ==\n\n")

# a and c as (n_0, n_half, n_1, n'_0) coefficient rows
A <- c(1, 11/2, 62, -28)
C <- c(1,    3, 12,  -8)

# Eliminate n'_0: the combination that kills it is A * C[4] - C * A[4].
k <- A * C[4] - C * A[4]
cat(sprintf("  8a - 28c coefficients (n_0, n_half, n_1, n'_0): %s\n",
            paste(sprintf("%g", k), collapse = ", ")))
stopifnot(abs(k[4]) < 1e-12)
# k is proportional to (n_0 + 2 n_half - 8 n_1); normalise on n_0
kk <- k / k[1]
cat(sprintf("  so the condition is n_0 + %g n_half + %g n_1 = 0,\n", kk[2], kk[3]))
cat(sprintf("  that is   n_0 + %g n_half = %g dim(G)\n\n", kk[2], -kk[3]))
stopifnot(abs(kk[2] - 2) < 1e-12, abs(kk[3] + 8) < 1e-12)

# With an emergent Higgs, n_0 = 0, so n_half = 4 dim(G) and, from c = 0, n'_0 = 3 dim(G).
nhalf_of <- function(d) 4 * d
nprime_of <- function(d) {
  # c = 0 with n_0 = 0: 3 n_half + 12 d - 8 n'_0 = 0
  (C[2] * nhalf_of(d) + C[3] * d) / -C[4]
}
cat(sprintf("  emergent Higgs (n_0 = 0):  n_half = 4 dim(G),  n'_0 = %g dim(G)\n\n",
            nprime_of(1)))
stopifnot(abs(nprime_of(1) - 3) < 1e-12)

## ---------------------------------------------------------------- the groups
# gen = Weyl fermions in one generation of that group, so the generation count is comparable
g <- data.frame(
  group = c("SU(3)xSU(2)xU(1)", "SU(5) with RH nu", "SO(10)", "Pati-Salam SU(4)xSU(2)xSU(2)",
            "left-right", "trinification SU(3)^3", "E6", "SM + U(1)_B-L"),
  dim   = c(12, 24, 45, 21, 15, 24, 78, 13),
  gen   = c(16, 16, 16, 16, 16, 27, 27, 16),
  stringsAsFactors = FALSE)
g$n_half   <- nhalf_of(g$dim)
g$n_prime  <- nprime_of(g$dim)
g$gens     <- g$n_half / g$gen
# What the FULL condition asks for if the group is forced to three generations with
# fundamental scalars instead: n_0 = 8 dim(G) - 2 n_half(3 generations).
g$scal3    <- 8 * g$dim - 2 * (3 * g$gen)

cat("  group                          dim   n_half  n'_0   generations   real scalars at 3 gens\n")
for (i in seq_len(nrow(g))) {
  cat(sprintf("  %-29s %4d %7d %5d %12s %12d\n", g$group[i], g$dim[i], g$n_half[i],
              g$n_prime[i],
              if (abs(g$gens[i] - round(g$gens[i])) < 1e-9) sprintf("%d", round(g$gens[i]))
              else sprintf("%.2f", g$gens[i]),
              g$scal3[i]))
}

exact <- g[abs(g$gens - round(g$gens)) < 1e-9, ]
cat(sprintf("\n  groups fitting a whole number of generations with nothing added: %s\n",
            paste(sprintf("%s (%d)", exact$group, round(exact$gens)), collapse = ", ")))

## ---------------------------------------------------------------- checks, each able to fail
ok <- TRUE
chk <- function(label, good) {
  cat(sprintf("    %-56s %s\n", label, if (good) "yes" else "NO"))
  good
}
cat("\n  checks:\n")
ok <- chk("the SM needs 48 Weyl fermions, which is 3 generations of 16",
          g$n_half[1] == 48 && round(g$gens[1]) == 3) && ok
ok <- chk("the SM needs 36 dimension-zero scalars, BFT's own number",
          g$n_prime[1] == 36) && ok
ok <- chk("48 - 12 equals 36 only because 4 dim - dim = 3 dim",
          g$n_half[1] - g$dim[1] == g$n_prime[1]) && ok
ok <- chk("the SM needs no fundamental scalars at three generations",
          g$scal3[1] == 0) && ok
ok <- chk("exactly two groups fit a whole number of generations",
          nrow(exact) == 2) && ok
ok <- chk("SU(5) fits only at six generations",
          round(g$gens[2]) == 6) && ok
ok <- chk("an ordinary Higgs doublet breaks the SM case (4 + 96 != 96)",
          4 + 2 * 48 != 8 * 12) && ok
# E6: three 27s, not three 16s. An earlier hand table used 48 here and got 528.
ok <- chk("E6 at three 27s wants 462 real scalars, not 528",
          g$scal3[7] == 462) && ok

cat(if (ok) "\n  all checks passed\n\n" else "\n  A CHECK FAILED\n\n")
if (!ok) quit(status = 1)
