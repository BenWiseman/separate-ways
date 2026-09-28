# The amplitude of the quantum half, made a number.
#
# Two results sit in this construction without being joined. The fold's parity is the
# classical-quantum split: Phi_c = (Phi + Theta Phi)/2 is the fold-even part and Phi_q =
# Phi - Theta Phi the fold-odd one, and those are the Keldysh classical and quantum variables. And
# at a bifurcate Killing horizon the fold's map is the HALF-PERIOD shift, so the cross-sheet
# correlator is the direct one at t - i beta/2, a thermofield double.
#
# Put the two together and the ratio of the two halves is not qualitative. It is a function of one
# dimensionless variable, and the variable is beta omega.

cat("=== 1. the two correlators, from KMS rather than from memory ===\n")
cat("   For a single mode the Wightman function is G(t) = (1/2 omega)[(1+n) e^{-i omega t} +\n")
cat("   n e^{i omega t}] with n the Bose factor. Its symmetric part at coincidence and its\n")
cat("   half-period shift at coincidence are\n")
cat("        W(0) = coth(x)/2 omega,        W(-i beta/2) = 1/(2 omega sinh x),     x = beta omega/2.\n")
Wd <- function(x) 1 / tanh(x) / 2            # times 1/omega
Wc <- function(x) 1 / (2 * sinh(x))          # times 1/omega
cat("   Checked against the mode sum with the Bose factor put in by hand:\n")
cat("      beta omega    W(0) direct    from n        W(-i b/2) direct   from n\n")
for (bw in c(0.5, 1.5, 2.19722, 4)) {
  x <- bw / 2; n <- 1 / (expm1(bw))
  d1 <- Wd(x); d2 <- (1 + 2 * n) / 2
  c1 <- Wc(x); c2 <- ((1 + n) * exp(-bw / 2) + n * exp(bw / 2)) / 2
  cat(sprintf("   %10.5f  %13.8f  %11.8f  %17.8f  %10.8f\n", bw, d1, d2, c1, c2))
  stopifnot(abs(d1 - d2) < 1e-12, abs(c1 - c2) < 1e-12)
}

cat("\n=== 2. the ratio of the quantum half to the classical half ===\n")
cat("   <Phi_c Phi_c> = (W(0) + W_cross)/2 and <Phi_q Phi_q> = 2(W(0) - W_cross), so\n")
cat("        R = <Phi_q^2>/<Phi_c^2> = 4 (cosh x - 1)/(cosh x + 1) = 4 tanh^2(beta omega / 4).\n")
R_assembled <- function(bw) { x <- bw / 2; 2 * (Wd(x) - Wc(x)) / ((Wd(x) + Wc(x)) / 2) }
R_closed    <- function(bw) 4 * tanh(bw / 4)^2
cat("      beta omega    assembled      closed form     difference\n")
for (bw in c(0.2, 1.0, 2.19722, 5.0, 20.0)) {
  cat(sprintf("   %10.5f  %13.9f  %14.9f  %.1e\n",
              bw, R_assembled(bw), R_closed(bw), abs(R_assembled(bw) - R_closed(bw))))
  stopifnot(abs(R_assembled(bw) - R_closed(bw)) < 1e-12)
}
cat("\n   The limits are the physics. At low beta omega, meaning a mode far below the horizon\n")
cat("   temperature, R goes as (beta omega)^2/4 and the quantum half switches off: that is the\n")
cat("   classical limit, and it is a limit in temperature rather than in hbar. At high beta omega\n")
cat("   R saturates at 4.\n")
for (bw in c(0.01, 0.1, 1, 10, 100)) {
  cat(sprintf("      beta omega = %6.2f   R = %10.6f   (small-bw form %.6f)\n",
              bw, R_closed(bw), bw^2 / 4))
}

cat("\n=== 3. where the two halves are equal, which is a pure number ===\n")
cat("   R = 1 needs tanh(beta omega/4) = 1/2, so beta omega = 4 artanh(1/2) = 2 ln 3.\n")
bw_eq <- 4 * atanh(0.5)
cat(sprintf("      4 artanh(1/2) = %.10f      2 ln 3 = %.10f      difference %.1e\n",
            bw_eq, 2 * log(3), abs(bw_eq - 2 * log(3))))
stopifnot(abs(bw_eq - 2 * log(3)) < 1e-12)
cat(sprintf("      and R there is exactly %.10f\n", R_closed(bw_eq)))
stopifnot(abs(R_closed(bw_eq) - 1) < 1e-12)

cat("\n=== 4. and it lands on a number the construction already had ===\n")
cat("   The bang's crossing ceiling arcsinh(1/sqrt 2) meets the horizon squeeze at beta omega =\n")
cat("   ln 3, where tanh r = e^{-beta omega/2} = 1/sqrt 3. The equal-amplitude point sits at\n")
cat("   exactly twice that frequency, where the squeeze is the SQUARE of the ceiling's:\n")
cat(sprintf("      ceiling meets squeeze at beta omega = ln 3   = %.8f,  tanh r = %.8f = 1/sqrt3\n",
            log(3), exp(-log(3) / 2)))
cat(sprintf("      halves are equal at beta omega = 2 ln 3      = %.8f,  tanh r = %.8f = 1/3\n",
            2 * log(3), exp(-2 * log(3) / 2)))
stopifnot(abs(exp(-log(3)/2)^2 - exp(-log(3))) < 1e-14)
cat("      (1/sqrt3)^2 = 1/3 exactly, so one squeeze is the square of the other.\n")
cat("   Both are the same statement about occupancy. arcsinh(1/sqrt2) is where the crossing's\n")
cat("   occupancy is one half; 2 ln 3 is where the fold's two halves carry equal weight.\n")

cat("\n=== 5. what this joins ===\n")
cat("   The classical-quantum divide and the horizon temperature stop being two results. The\n")
cat("   divide is the fold's parity, the parity's two halves have amplitudes fixed by the\n")
cat("   cross-sheet correlator, and that correlator is thermal because the fold's map is a\n")
cat("   half-period shift. So the weight of the quantum half is 4 tanh^2(beta omega/4): a\n")
cat("   temperature, exactly, with no free parameter and no scale put in. The classical world is\n")
cat("   what a horizon looks like below its own temperature.\n")

cat("\n=== 6. the plants ===\n")
cat("   (a) the closed form must fail if the shift is not a HALF period.\n")
for (f in c(0.25, 0.5, 0.75)) {
  x <- 2.19722 / 2
  Wcf <- ((1 + 1/expm1(2.19722)) * exp(-2 * f * x) + (1/expm1(2.19722)) * exp(2 * f * x)) / 2
  Rf <- 2 * (Wd(x) - Wcf) / ((Wd(x) + Wcf) / 2)
  cat(sprintf("      shift = %.2f of a period:  R = %10.6f   %s\n", f, Rf,
              ifelse(abs(f - 0.5) < 1e-9, "<- the fold's value, and R = 1 here", "")))
  if (abs(f - 0.5) > 1e-9) stopifnot(abs(Rf - 1) > 0.05)
}
cat("   Only the half-period shift puts the equal-amplitude point at 2 ln 3. The result is the\n")
cat("   fold's, not an artefact of thermality alone.\n")
cat("   (b) the equal-amplitude condition must move if R's form is altered.\n")
for (k in c(2, 4, 8)) {
  r <- tryCatch(uniroot(function(b) k * tanh(b/4)^2 - 1, c(1e-6, 50))$root, error = function(e) NA)
  cat(sprintf("      R = %d tanh^2(bw/4):  equal at beta omega = %.6f   %s\n", k, r,
              ifelse(k == 4, "= 2 ln 3", "")))
}
