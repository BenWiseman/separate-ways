# A.13 forces the seam transparent by Hadamard regularity, which is a regularity argument.
# Is there an independent, information-theoretic characterisation of the same point?
#
# Condensed matter supplies the tool: a free-fermion chain with a DEFECT BOND is the same
# two-port scattering problem, with transmission T = 4 td^2/(1+td^2)^2 at the band centre,
# and the entanglement across the defect is computable exactly.
#
# FIRST ATTEMPT, RECORDED BECAUSE IT WAS WRONG. I computed the total entropy S of half the
# chain and expected it to be a function of T peaking at T = 1. It is not a function of T
# at all: S(td=0.5) = 0.62 and S(td=2.0) = 1.47 share T = 0.64. A strong bond builds a
# dimer across the cut whose local contribution swamps the transmission-dependent part,
# and S rises monotonically with td instead of peaking. Both conclusions I first drew from
# that table were false.
#
# The standard statement concerns the COEFFICIENT of log L, not the total, so extract that.

Sent <- function(L, td) {
  h <- matrix(0, L, L)
  for (i in 1:(L-1)) { t <- if (i == L/2) td else 1; h[i,i+1] <- -t; h[i+1,i] <- -t }
  e <- eigen(h, symmetric=TRUE); occ <- order(e$values)[1:(L/2)]
  V <- e$vectors[, occ, drop=FALSE]; C <- V %*% t(V)
  nu <- eigen(C[1:(L/2),1:(L/2)], symmetric=TRUE)$values
  nu <- pmin(pmax(nu,1e-14),1-1e-14); -sum(nu*log(nu)+(1-nu)*log(1-nu)) }
Tof  <- function(td) 4*td^2/(1+td^2)^2
Ls   <- c(64,128,256,512)
slope<- function(td) coef(lm(sapply(Ls, function(L) Sent(L,td)) ~ log(Ls)))[2]

cat("=== 1. the total entropy is NOT a function of transmission\n\n")
cat("        td        T          S(L=200)\n")
for (td in c(0.5, 2.0, 0.7, 1/0.7)) cat(sprintf("   %9.3f %9.5f %13.5f\n", td, Tof(td), Sent(200,td)))
cat("\n  Same T, different S, by more than a factor of two. The local dimer dominates.\n")

cat("\n=== 2. the log L coefficient IS a function of transmission\n\n")
s1 <- slope(1.0)
cat("        td        1/td         T        slope(td)/slope(1)   slope(1/td)/slope(1)\n")
for (td in c(0.3,0.5,0.7,0.9)) 
  cat(sprintf("   %9.2f %9.2f %10.5f %19.4f %21.4f\n", td, 1/td, Tof(td),
      slope(td)/s1, slope(1/td)/s1))
cat("\n  Symmetric under td -> 1/td to about one per cent, which is the inversion the seam's\n")
cat("  own rapidity parametrisation carries (A.13). So the scaling part depends on T alone.\n")

cat("\n=== 3. and it is maximal at the transparent point\n\n")
cat("        td          T        coefficient, normalised to T = 1\n")
for (td in c(0.8,0.9,0.95,1.0,1.05,1.1,1.25)) 
  cat(sprintf("   %9.2f %10.5f %24.5f\n", td, Tof(td), slope(td)/s1))
cat("\n  Peaks at td = 1 with both neighbours below it, so the transparent seam maximises\n")
cat("  the scaling part of the entanglement across the seam. The absolute normalisation\n")
cat("  is not quoted: half of an OPEN chain carries c/6 rather than c/3, so the raw slope\n")
cat("  is a convention. The shape is the result.\n")

cat("\n=== 4. how sharp is the selection? SOFT.\n\n")
cat("        |r|          T        fraction of the maximum\n")
for (td in c(1.0,0.95,0.9,0.8,0.5,0.3)) { Tt <- Tof(td)
  cat(sprintf("   %9.4f %10.5f %22.4f\n", sqrt(max(1-Tt,0)), Tt, slope(td)/s1)) }
cat("\n  A five per cent reflection keeps 99.7 per cent of the coefficient. So this\n")
cat("  identifies the transparent point as the maximum WITHOUT excluding its\n")
cat("  neighbourhood. It is a variational characterisation, not an exclusion, and A.13's\n")
cat("  regularity argument remains the thing that does the excluding.\n")

cat("\n=== 5. what it adds, and what it is not\n\n")
cat("  ADDS: a second and independent criterion landing on the same point. Hadamard\n")
cat("  regularity says a constant reflectivity is inadmissible; this says the transparent\n")
cat("  seam maximises the entanglement between the two sides. For a framework whose claim\n")
cat("  is that the sheets are joined, that is the natural way to say how strongly, and it\n")
cat("  is a computed quantity rather than a figure of speech.\n")
cat("  IS NOT: a gravitational calculation. A free-fermion chain with a defect bond shares\n")
cat("  the scattering problem and not the dynamics, and the entropy is the lattice\n")
cat("  model's, not the seam's.\n")
