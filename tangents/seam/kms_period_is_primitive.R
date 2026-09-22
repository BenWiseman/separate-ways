# alpha_is_kms.R sweeps beta and concludes the match "holds at 2 pi/H and fails everywhere
# else". Its own table prints a second MATCH at 3 x (2 pi/H). Work out which betas match and
# what actually selects the fundamental.
H <- 1
Wt     <- function(t,b) 1/sinh(H*t/2 - 1i*H*b/4)^2
target <- function(t) -1/cosh(H*t/2)^2
ts <- c(0.4, 0.9, 1.7, 3.1)
dev <- function(f) max(sapply(ts, function(t) Mod(Wt(t, f*2*pi/H) - target(t))))

cat("  beta/(2pi/H)    max deviation      matches?\n")
for (f in c(1,2,3,4,5,6,7)) cat(sprintf("  %10d %18.3e   %s\n", f, dev(f),
      if (dev(f) < 1e-12) "YES" else "no"))
cat("\n  Every ODD multiple matches; no even one does. The reason is exact:\n")
cat("    sinh(z - i(2j+1)pi/2) = (-1)^j i cosh(z),  and squaring removes the sign,\n")
cat("  so W(t - i beta/2) returns -1/cosh^2(Ht/2) at beta = (2j+1) 2pi/H for every j.\n")
cat("  The match therefore determines beta only MODULO odd multiples. It does not by itself\n")
cat("  pick out 2pi/H, and the script's verdict line overstated what its own table showed.\n")

cat("\n  What selects the fundamental. A KMS state at inverse temperature beta requires the\n")
cat("  correlation function to be analytic in the strip 0 < Im t < beta and to satisfy the\n")
cat("  boundary condition across it. A function with period 2pi/H trivially also has period\n")
cat("  3 x 2pi/H, so the larger values are not independent temperatures: they are multiples of\n")
cat("  the PRIMITIVE period, and the KMS temperature is the primitive one. Check that the\n")
cat("  candidate strips below the fundamental contain no match, which is what primitivity needs:\n\n")
cat("     beta/(2pi/H) in (0,1)    max deviation\n")
for (f in c(0.2,0.4,0.6,0.8,0.95)) cat(sprintf("  %14.2f %22.3e\n", f, dev(f)))
cat("\n  No sub-fundamental match, so 2pi/H IS primitive and T_GH = H/2pi follows. The honest\n")
cat("  statement is that the match plus primitivity determines beta, not the match alone.\n")
