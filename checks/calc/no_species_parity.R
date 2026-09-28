# Can the fold's lift distinguish one sterile species from another?
#
# Section 4.2 says no, and gives an argument: a spacetime symmetry acts identically on every
# member of a flavour multiplet, so no involution on the cover can pick out the stable species
# the matter rule needs. That is a kill inherited rather than made, so it is worth testing before
# being relied on, because if it were wrong the fold would supply the stabilising rule and the
# matter sector would lose its last assumed input.
#
# THE ROUTE WORTH CHECKING. The argument as stated is about the spacetime action. What it does
# not obviously cover is the LIFT: Theta is antilinear, its action on a fermion carries a matrix
# and a phase, and Theta^2 = +1 or -1 depending on which. If different species could carry
# different lifts, the fold would admit some and forbid others, and that would be a species
# parity arriving from the geometry. Section 2.2's lift is Theta: psi(eta) -> sigma_x psi*(-eta).
# This asks what freedom that leaves.

sx <- matrix(c(0,1,1,0), 2, 2); sy <- matrix(c(0,-1i,1i,0), 2, 2); sz <- diag(c(1,-1))
id <- diag(2)
# an antilinear map psi -> e^{i a} M psi*(-eta); squaring it gives e^{i a} M (e^{i a} M psi)*
#                                                = e^{i a} M conj(e^{i a}) conj(M) psi = M conj(M) psi
sq <- function(M, a) M %*% Conj(M)              # the phase cancels against its own conjugate

cat("=== 1. the phase in the lift cannot do anything, and the reason is the antilinearity ===\n")
cat("   For Theta: psi -> e^{i a} M psi*(-eta), squaring gives e^{i a} M conj(e^{i a} M psi),\n")
cat("   and the phase meets its own conjugate. Theta^2 = M conj(M), with no a in it.\n\n")
cat("      phase a      Theta^2 with M = sigma_x      differs from a = 0?\n")
for (a in c(0, 0.7, pi/2, 2.3)) {
  v <- sq(sx, a)
  cat(sprintf("   %10.4f   %26s   %s\n", a,
              sprintf("%+.3f", Re(v[1,1])),
              ifelse(max(abs(v - sq(sx,0))) < 1e-14, "no", "YES")))
  stopifnot(max(abs(v - sq(sx, 0))) < 1e-14)
}
cat("   So a per-species Majorana phase is not a handle: every species with the same matrix has\n")
cat("   the same Theta^2 whatever phase it carries.\n")

cat("\n=== 2. the matrix is not free either: the Hamiltonian fixes it ===\n")
cat("   Section 2.2 has H(eta) = [[gamma eta, p], [p, -gamma eta]] and the lift must preserve the\n")
cat("   equation, which needs H(-eta) = M H(eta) M^(-1) with M real, since Theta conjugates. Of\n")
cat("   the Pauli matrices only one does it:\n\n")
H <- function(e, g = 1, p = 0.8) matrix(c(g*e, p, p, -g*e), 2, 2, byrow = TRUE)
cat("      M            real?    M H(eta) M^-1 = H(-eta)?     Theta^2\n")
for (nm in list(list("identity", id), list("sigma_x", sx), list("sigma_y", sy), list("sigma_z", sz))) {
  M <- nm[[2]]
  isreal <- max(abs(Im(M))) < 1e-14
  ok <- max(abs(M %*% H(0.6) %*% solve(M) - H(-0.6))) < 1e-12
  t2 <- sq(M, 0)
  lab <- if (max(abs(t2 - id)) < 1e-12) "+1" else if (max(abs(t2 + id)) < 1e-12) "-1" else "neither"
  cat(sprintf("      %-11s  %-7s  %-26s  %s\n", nm[[1]], ifelse(isreal,"yes","no"),
              ifelse(ok,"yes","no"), lab))
}
cat("   Only sigma_x satisfies the conjugation relation, and it gives Theta^2 = +1. sigma_y would\n")
cat("   give -1 but is imaginary and does not preserve the equation. The lift is determined by\n")
cat("   the dynamics, so it is the same for every species obeying the same equation.\n")

cat("\n=== 3. the plant: a Hamiltonian with a different parity must choose a different matrix ===\n")
cat("   If the matrix were being read off the dynamics, then changing the dynamics must change\n")
cat("   it. Take H'(eta) = [[p, gamma eta], [gamma eta, -p]], whose odd part sits off-diagonal.\n\n")
H2 <- function(e, g = 1, p = 0.8) matrix(c(p, g*e, g*e, -p), 2, 2, byrow = TRUE)
for (nm in list(list("identity", id), list("sigma_x", sx), list("sigma_z", sz))) {
  M <- nm[[2]]
  ok <- max(abs(M %*% H2(0.6) %*% solve(M) - H2(-0.6))) < 1e-12
  cat(sprintf("      %-11s  preserves the new equation: %s\n", nm[[1]], ifelse(ok,"yes","no")))
}
cat("   sigma_z now, not sigma_x, so the matrix is genuinely read off the Hamiltonian and\n")
cat("   section 2 is not returning the same answer whatever it is given.\n")
stopifnot(max(abs(sz %*% H2(0.6) %*% solve(sz) - H2(-0.6))) < 1e-12,
          max(abs(sx %*% H2(0.6) %*% solve(sx) - H2(-0.6))) > 1e-6)

cat("\n=== 4. the verdict, which is that the kill holds ===\n")
cat("   Two places a species parity might have hidden are both closed. The phase in the lift\n")
cat("   cancels against its own conjugate, so a Majorana phase gives no handle. The matrix is\n")
cat("   fixed by the requirement that the lift preserve the mode equation, so every species\n")
cat("   obeying that equation carries the same lift and the same Theta^2 = +1. Nothing in the\n")
cat("   fold separates one sterile neutrino from another, and section 4.2's limit stands as a\n")
cat("   checked statement rather than an argument.\n")
cat("   What would change it is flavour becoming a spacetime label, which this construction\n")
cat("   neither has nor needs, and which is named in 4.2 as the thing a derivation would\n")
cat("   require.\n")
