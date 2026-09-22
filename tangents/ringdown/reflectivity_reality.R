# The companion says a real reflectivity preserves the QNM mirror symmetry and any imaginary
# part breaks it. Conjugating F(w) + R G(w) = 0 gives F(-w*) + conj(R) G(-w*) = 0, so -w* solves
# the SAME equation when R(-w*) = conj(R(w)). For a CONSTANT R that is "R is real". For a
# frequency-dependent R it is not. Test the delayed response R(w) = r0 exp(i w tau).
cond <- function(R, w) Mod(R(-Conj(w)) - Conj(R(w)))     # zero iff the mirror symmetry survives
ws <- c(0.37-0.089i, 0.35-0.27i, 0.8+0i, 1.2-0.5i)

cat("  R(w)                         max |R(-w*) - conj(R(w))|   preserves mirror?\n")
tests <- list(
  "constant real 0.01"      = function(w) 0.01 + 0i,
  "constant imaginary 0.01i"= function(w) 0+0.01i,
  "delayed 0.01 exp(i w 2)" = function(w) 0.01*exp(1i*w*2),
  "delayed 0.01 exp(i w 7)" = function(w) 0.01*exp(1i*w*7),
  "0.01 (1 + i w)"          = function(w) 0.01*(1+1i*w),
  "0.01 (1 + i) w"          = function(w) 0.01*(1+1i)*w
)
for (nm in names(tests)) {
  R <- tests[[nm]]; d <- max(sapply(ws, function(w) cond(R,w)))
  cat(sprintf("  %-28s %22.3e   %s\n", nm, d, if (d < 1e-14) "YES" else "no"))
}
cat("\n  A delayed response is complex at every real frequency and still preserves the symmetry,\n")
cat("  so 'real reflectivity' is sufficient but NOT necessary. The companion's numerical\n")
cat("  illustration uses a constant R, where the two conditions coincide; the general statement\n")
cat("  it drew from that illustration does not hold.\n")
cat("\n  Note 0.01(1 + i w) also passes: any R with real Taylor coefficients in (i w) does.\n")
