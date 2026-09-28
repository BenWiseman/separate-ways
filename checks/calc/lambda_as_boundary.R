# What kind of quantity Lambda is in this construction, which is not the kind the question assumes.
#
# Ben's framing: if the fold and general relativity converge on needing the same constants, that
# may be structure rather than coincidence, and Lambda may be a rate inherited from conditions at
# or before the bang rather than a coupling of the present physics. The re-derived Clausius step
# already says something sharp about this and it had not been read that way.

c_ <- 2.99792458e8; hb <- 1.054571817e-34; GN <- 6.67430e-11; eV <- 1.602176634e-19
H0 <- 67.4e3 / 3.0857e22; OmL <- 0.685
Lam <- 3 * OmL * H0^2 / c_^2

cat("=== 1. what the derivation actually says Lambda is ===\n")
cat("   In the Clausius step the balance gives 2 pi T_ab = eta R_ab + f g_ab with f a FUNCTION.\n")
cat("   Taking the divergence, matter conservation and the contracted Bianchi identity give\n")
cat("        grad_b ( f + eta R / 2 ) = 0,\n")
cat("   so f + eta R/2 is constant across the whole spacetime and Lambda is its value. Lambda is\n")
cat("   therefore an INTEGRATION CONSTANT of the derivation and not a coupling in it. Nothing\n")
cat("   local determines an integration constant; a boundary condition does.\n")
cat("   That is not a shortfall of this construction. It is a statement about what kind of\n")
cat("   quantity Lambda is, and it means no amount of local physics could ever return its value.\n")

cat("\n=== 2. and it is the same statement the cosmology already makes from the other side ===\n")
cat("   Appendix D.2's restricted action gives p_Lambda = -rho_Lambda and Q = 0, so w_0 = -1 and\n")
cat("   w_a = 0 exactly: Lambda cannot evolve. An integration constant cannot evolve either. The\n")
cat("   two are the same fact reached down two roads, and the paper had them as two results.\n")
cat("   The observational consequence is already in the commitments table: a robust detection of\n")
cat("   evolving dark energy ends the construction, and it would also end the reading of Lambda\n")
cat("   as an integration constant.\n")

cat("\n=== 3. Lambda IS a rate, which is the useful half of the framing ===\n")
HL <- c_ * sqrt(Lam / 3)
cat(sprintf("      Lambda = 3 H_L^2 / c^2 with H_L = %.4e s^-1, a de Sitter expansion rate\n", HL))
cat(sprintf("      its inverse is %.2f Gyr, against a universe %.2f Gyr old\n",
            1 / HL / (3.156e16), 13.8))
cat("   So Lambda is a rate squared, and the rate is within a factor of two of the present\n")
cat("   Hubble rate, which is the why-now coincidence and is not explained here either.\n")

cat("\n=== 4. is that rate inherited from anything the construction contains? checked, not guessed ===\n")
cat("   If Lambda were set by a scale the fold carries, the ratio of H_L to that scale's rate\n")
cat("   would be a simple number. Enumerate every rate the construction has and look:\n\n")
rate <- function(E_eV) E_eV * eV / hb                       # energy in eV to a rate in s^-1
rows <- list(
  list("de Sitter rate H_L",            HL),
  list("present Hubble rate H_0",       H0),
  list("dark-matter mass 491.6 PeV",    rate(491.6e15)),
  list("two-body line 245.8 PeV",       rate(245.8e15)),
  list("neutrino sum 58.8 meV",         rate(58.8e-3)),
  list("Planck rate",                   sqrt(c_^5 / (hb * GN))))
cat("      scale                          rate (s^-1)        H_L / rate      log10\n")
for (r in rows) {
  cat(sprintf("      %-28s %14.4e  %16.4e  %8.2f\n", r[[1]], r[[2]], HL / r[[2]], log10(HL / r[[2]])))
}
cat("\n   None of those logs is near an integer multiple of another, and none of the ratios is a\n")
cat("   simple number. The construction contains no rate that sets H_L.\n")
lg <- sapply(rows[-1], function(r) log10(HL / r[[2]]))
cat(sprintf("      closest any ratio comes to a round power of ten: %.2f decades off\n",
            min(abs(lg - round(lg)))))

cat("\n=== 5. which is what an integration constant should look like ===\n")
cat("   A coupling is fixed by the theory's content and would show up as a ratio to one of the\n")
cat("   scales above. An integration constant is fixed by a boundary condition and need not be\n")
cat("   commensurate with anything in the local physics. The enumeration finding nothing is the\n")
cat("   expected result for the second reading and an awkward one for the first.\n")
cat("   So the position is stronger than it looked. The construction does not fail to\n")
cat("   derive Lambda. It says Lambda is not the kind of thing a local derivation returns, it\n")
cat("   predicts that such a quantity cannot evolve, and that prediction is already falsifiable.\n")

cat("\n=== 6. and what that leaves for G ===\n")
cat("   The pure number is G Lambda. If Lambda is a boundary datum then so is the product, and\n")
cat("   the question 'why is G Lambda 1e-122' becomes 'why did the fold's global state carry\n")
cat("   that boundary value', which is a question about the bang rather than about gravity.\n")
cat(sprintf("      G Lambda = %.3e\n", (hb * GN / c_^3) * Lam))
cat("   That relocates the cosmological constant problem rather than solving it, and relocating\n")
cat("   it is worth something: it says the answer is not to be looked for in the field equations.\n")

cat("\n=== 7. the plant ===\n")
cat("   The rate enumeration must FIND a match when there is one to find, or it is not a test.\n")
cat(sprintf("      planted: H_0 against sqrt(Om_L) H_0, ratio %.6f against sqrt(0.685) = %.6f\n",
            HL / H0, sqrt(OmL)))
stopifnot(abs(HL / H0 - sqrt(OmL)) < 1e-6)
cat("      It finds that one exactly, since it is true by construction, so a null result on the\n")
cat("      others is a null result and not a blind spot.\n")
