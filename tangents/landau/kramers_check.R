# SUPERSEDED 2026-09-20 AND WRONG. This file reports Theta^2 = +1 on the whole Fock space.
# It is not. build() writes the transformed basis vectors into swapped columns, and the correct
# answer for c1 -> c2, c2 -> -c1 is Theta^2 = (-1)^F, which follows in one line from applying the
# map twice. Use kramers_recheck.R, which builds both exchange conventions explicitly. Kept so the
# error stays on the record rather than vanishing from it.
# Attack A.18 where it is weakest. The whole appendix rests on Theta being an involution on
# the states it constrains. For fermions an antiunitary CPT normally obeys Theta^2 = (-1)^F,
# and Kramers' theorem then FORBIDS any invariant state in the odd sector. If the paired
# vacuum sits there, A.18 is dead. Build the Fock space and check, do not argue.
#
# Two modes, one per sheet. Jordan-Wigner on |00>,|10>,|01>,|11>.
I2 <- diag(2); sz <- diag(c(1,-1)); sm <- matrix(c(0,0,1,0),2,2)   # sm = lowering
kron <- function(a,b) a %x% b
c1 <- kron(sm, I2)            # annihilates mode 1
c2 <- kron(sz, sm)            # mode 2, with the JW string
vac <- c(1,0,0,0)
d1 <- t(Conj(c1)); d2 <- t(Conj(c2))                                # creation operators
lbl <- c("|00>","|10>","|01>","|11>")
Fpar <- c(0,1,1,2)

cat("  anticommutators, as a sanity check on the construction:\n")
cat(sprintf("   {c1,c1*} = I ? %s    {c2,c2*} = I ? %s    {c1,c2} = 0 ? %s\n",
  all(abs(c1%*%d1 + d1%*%c1 - diag(4)) < 1e-12),
  all(abs(c2%*%d2 + d2%*%c2 - diag(4)) < 1e-12),
  all(abs(c1%*%c2 + c2%*%c1) < 1e-12)))

# Theta = U K, antilinear. Build U from the action on the creation operators.
# Convention A (naive swap):      c1 -> c2,  c2 -> c1
# Convention B (fermionic, T^2=-1): c1 -> c2,  c2 -> -c1
build <- function(sign2) {
  U <- matrix(0,4,4)
  U[,1] <- vac                                     # Theta|00> = |00>
  U[,2] <- as.vector(d2 %*% vac)                   # Theta c1* |00> = c2* |00>
  U[,3] <- as.vector(sign2 * (d1 %*% vac))         # Theta c2* |00> = +-c1* |00>
  U[,4] <- as.vector(d2 %*% (sign2 * d1) %*% vac)  # Theta c1*c2*|00> = c2*(+-c1*)|00>
  U
}
for (nm in c("A  naive  c2 -> +c1", "B  fermionic c2 -> -c1")) {
  s2 <- if (substr(nm,1,1)=="A") 1 else -1
  U <- build(s2); T2 <- U %*% Conj(U)
  cat(sprintf("\n  convention %s\n", nm))
  cat(sprintf("    Theta^2 = %s\n",
     if (all(abs(T2-diag(4))<1e-12)) "+1 on every state"
     else if (all(abs(T2+diag(4))<1e-12)) "-1 on every state"
     else paste0("diag(", paste(round(Re(diag(T2)),3), collapse=","), ")  by state ", paste(lbl,collapse=" "))))
  # does the paired vacuum survive?
  for (n in c(0.10, 0.25, 0.4999)) {
    psi <- rep(0,4); psi[1] <- sqrt(1-n); psi[4] <- sqrt(n)
    img <- as.vector(U %*% Conj(psi))
    ov  <- abs(sum(Conj(img)*psi))
    cat(sprintf("    n=%.4f  |<Theta psi | psi>| = %.6f   %s\n", n, ov,
        if (abs(ov-1) < 1e-10) "INVARIANT" else "NOT invariant"))
  }
}
cat("\n=== reading\n")
cat("  The pair state sqrt(1-n)|00> + sqrt(n)|11> has EVEN fermion number. If Theta^2 = -1\n")
cat("  on the odd sector only, Kramers forbids invariant one-particle states but says nothing\n")
cat("  about this one, and A.18 survives. If instead the naive convention is forced, the pair\n")
cat("  is not invariant for any n>0 and A.18 is dead. The table above decides it.\n")

cat("\n=== challenge from an independent ideation pass: does Theta^2 = +1 fail in the ODD\n")
cat("    fermion-parity sector? Print Theta^2 state by state rather than trusting a norm.\n\n")
for (nm in c("A naive", "B fermionic")) {
  s2 <- if (substr(nm,1,1)=="A") 1 else -1
  U <- build(s2); T2 <- U %*% Conj(U)
  cat(sprintf("  convention %s : Theta^2 diagonal by state\n", nm))
  for (k in 1:4)
    cat(sprintf("     %-6s  F=%d   (Theta^2)_kk = %+.3f   off-diagonal max %.1e\n",
        lbl[k], Fpar[k], Re(T2[k,k]), max(abs(T2[k,-k]))))
  ev <- eigen(T2)$values
  cat(sprintf("     eigenvalues: %s\n\n", paste(sprintf("%+.3f", Re(ev)), collapse=", ")))
}
cat("  If every diagonal entry is +1 and the off-diagonals vanish, Theta^2 = +1 on the WHOLE\n")
cat("  Fock space including the odd sector, and the graded-class claim does not arise here.\n")
cat("  Note what this does and does not settle: it is a statement about THIS two-mode model\n")
cat("  with the pairing A.18 uses, not about a general CPT operator on a full fermionic field,\n")
cat("  where Theta^2 = (-1)^F is standard. A.18 already scopes itself to the mode statement.\n")
