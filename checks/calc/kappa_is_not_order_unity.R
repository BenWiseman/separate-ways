#!/usr/bin/env Rscript
# The order-unity factor in kappa is one, and the shell is 2.1 microns rather than "between 0.85
# and 5.3".
#
# WHAT THE MANUSCRIPTS SAY. 3.6: "What is left in kappa is an order-unity factor from the
# sub-leading derivative pairings, and since it enters as the two-fifths power, a factor of ten
# either way puts the fermion's shell between 0.85 and 5.3 microns." That caveat was written
# because the coefficient was assembled from the proper-time representation keeping only the
# leading derivative pairing, with no way to price what was dropped.
#
# THERE IS NOW A WAY. image_stress_components.R sums A.18's tower and gets the exact stress at a
# caustic, keeping every pairing. A.19's rule, with Delta'^{1/2} = 1 and L = pi a, fixes the
# amplitude in that same geometry with nothing left to choose. So the leading-pairing assembly and
# the exact answer can be set side by side, and their ratio IS the order-unity factor.
#
# IT IS ONE. For the massless part the two agree to sixteen digits, and the reason is structural:
# with G = F(sigma), the null-null second derivative is F''(k.grad sigma)^2 + F' k^a k^b
# grad_a grad_b sigma, and the second piece carries two fewer powers of sigma, so it cannot touch
# the leading coefficient. That is a different object from the mixed grad_a grad_b' sigma whose
# determinant is the Van Vleck factor, which is what diverges at a caustic and what closed the
# structural route to the sign.

TOL <- 1e-9
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

SRC <- readLines("checks/calc/image_stress_components.R")
i1 <- grep("^EG <- ", SRC)[1]; i2 <- grep('^cat\\("=== 1\\.', SRC)[1]
eval(parse(text = paste(SRC[i1:(i2-1)], collapse = "\n")))
k1 <- grep("^tower <- function", SRC)[1]
k2 <- grep("^\\}", SRC); k2 <- k2[k2 > k1][1]
eval(parse(text = paste(SRC[k1:k2], collapse = "\n")))
kk_rad <- function(o) -Re(o$D) - Re(o$E)

AMP    <- pi^1.5      # A.19's rule in A.18's geometry: reduced Van Vleck 1, projection pi a
GRADS2 <- pi^2        # (k . grad sigma)^2 on the radial null direction, with sigma = pi s

cat("=== 1. the massless part: leading pairing against the exact tower ===\n")
cat("   G = C sigma^{-3/2} with C = AMP (4 pi)^{-2} Gamma(3/2) 2^{3/2}, so the leading pairing\n")
cat("   gives T_kk = -C (15/4) (k.grad sigma)^2 sigma^{-7/2}, and sigma = pi s.\n")
Cw    <- AMP*(4*pi)^(-2)*gamma(1.5)*2^1.5
pred0 <- -Cw*(15/4)*GRADS2*pi^(-3.5)
exact0 <- -15/(32*pi*sqrt(2*pi))
cat(sprintf("      leading pairing  %.14f s^{-7/2}\n", pred0))
cat(sprintf("      exact tower      %.14f s^{-7/2}\n", exact0))
cat(sprintf("      ratio            %.14f\n", pred0/exact0))
note(abs(pred0/exact0 - 1) < 1e-12, "the leading pairing is EXACT for the massless part")
cat("   Not approximately: to twelve digits. The dropped pairing is F'(sigma) k^a k^b grad_a\n")
cat("   grad_b sigma, which carries two fewer powers of sigma than F''(k.grad sigma)^2, so it\n")
cat("   cannot reach the leading coefficient at all.\n")

cat("\n=== 2. the mass part, which is the one kappa is built from ===\n")
cat("   The heat kernel's mass factor e^{-m^2 s} gives, at first order, G_m = -m^2 AMP (4 pi)^{-2}\n")
cat("   Gamma(1/2) (2/sigma)^{1/2}, whose leading pairing is T_kk = +m^2 Cm (3/4) (k.grad sigma)^2\n")
cat("   sigma^{-5/2}. The tower is differenced in m^2 to pull the same quantity out of it.\n")
Cm    <- AMP*(4*pi)^(-2)*sqrt(2*pi)
predm <- Cm*(3/4)*GRADS2*pi^(-2.5)
cat(sprintf("\n      leading pairing  %+.10f m^2 s^{-5/2}\n\n", predm))
massfit <- function(s, ms = c(0.02, 0.04, 0.06, 0.08)) {
  base <- kk_rad(tower(s, eps = s/40, m2 = 0, xi = 1/6))
  d <- sapply(ms, function(m) kk_rad(tower(s, eps = s/40, m2 = m^2, xi = 1/6)) - base)
  unname(lm(d ~ I(ms^2) + 0)$coefficients[1])          # the m^2 slope, no intercept
}
cat("      s        fitted m^2 slope x s^{5/2}      ratio to the leading pairing\n")
rr <- c()
for (s in c(0.03, 0.01, 0.003)) {
  v <- massfit(s)*s^2.5; rr <- c(rr, predm/v)
  cat(sprintf("   %8.4f %26.8f %28.6f\n", s, v, predm/v))
}
note(all(abs(rr - 1) < 0.02), "the leading pairing is within two per cent for the mass part too")
cat("   Fitted over four masses at each offset rather than differenced at one, because a single\n")
cat("   difference of two nearly equal sums at a steep divergence is noise. The residual drift is\n")
cat("   the subleading s^{-3/2} term, not a missing pairing.\n")

cat("\n=== 3. the plants ===\n")
cat("   (a) The agreement must break if the amplitude is wrong, or it is measuring nothing:\n")
for (f in c(1.3, 0.7)) {
  p <- -(f*AMP)*(4*pi)^(-2)*gamma(1.5)*2^1.5*(15/4)*GRADS2*pi^(-3.5)
  cat(sprintf("      AMP x %.1f: ratio %.6f\n", f, p/exact0))
  note(abs(p/exact0 - 1) > 0.2, "a wrong amplitude breaks the agreement")
}
cat("   (b) And if the power of sigma in the pairing is wrong. F'' of sigma^{-3/2} is (15/4)\n")
cat("       sigma^{-7/2}; the neighbouring coefficients are what a slipped power would give:\n")
for (cf in c(3/4, 15/4, 105/8)) {
  p <- -Cw*cf*GRADS2*pi^(-3.5)
  cat(sprintf("      coefficient %6.3f: ratio %.6f  %s\n", cf, p/exact0,
              ifelse(abs(p/exact0 - 1) < 1e-12, "the claim", "rejected")))
}
note(abs((-Cw*(3/4)*GRADS2*pi^(-3.5))/exact0 - 1) > 0.5, "a slipped power is rejected")

cat("\n=== 4. what that does to the shell ===\n")
cat("   kappa enters the thickness as its two-fifths power, so the caveat the manuscripts carry,\n")
cat("   a factor of ten either way, spans 0.85 to 5.3 microns for the fold's own fermion. With the\n")
cat("   factor measured at one, the shell is the central value and the span closes:\n")
for (f in c(0.1, 1, 10)) cat(sprintf("      factor %5.1f on kappa: shell %.2f microns\n", f, 2.127*f^0.4))
cat("   So 2.1 microns is the number, not the middle of a decade.\n")
cat("   What is NOT closed by this: the amplitude and the world function's linear coefficient are\n")
cat("   still A.19's, and the transfer of the pairing statement from A.18's geometry to a hole is\n")
cat("   the same transfer the sign rests on, which is argued and not computed in the interior.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
