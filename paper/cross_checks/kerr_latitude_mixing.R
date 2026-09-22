# kerr_latitude_mixing.R -- the latitude-dependent horizon response and its
# l -> l +/- 2 selection rule, App. A.11. Base R only, units M = 1.
#
# The Kerr bifurcation surface carries
#     ds^2 = rho_+^2 dth^2 + [(r_+^2+a^2)^2 sin^2 th / rho_+^2] dphi^2,
#     rho_+^2 = r_+^2 + a^2 cos^2 th,     r_+ = 1 + sqrt(1-a^2).
# Its local anisotropy, azimuthal proper length per unit coordinate over
# meridional, is
#     A(th) = sqrt(g_phph) / (sin th sqrt(g_thth)) = (r_+^2+a^2)/(r_+^2+a^2 cos^2 th),
# equal to 1 at every latitude only when a = 0.
#
# A horizon boundary condition weighted by A cannot be imposed mode by mode: it
# couples different l at fixed m. A is even in cos th, so its Legendre expansion
# has even terms only, and the coupling is l -> l +/- 2, l +/- 4, ... and never
# l +/- 1. The strength is set by c_2/c_0 and vanishes identically at a = 0.
#
# The same A also controls A.10: g_thth(pole)/g_thth(equator) = A(pi/2)/A(0),
# so the squashing that removes the closed-form kernel ratio and the mixing that
# moves multipoles are one function, not two coincidences.

rp <- function(a) 1 + sqrt(1 - a^2)
g_thth <- function(th, a) rp(a)^2 + a^2*cos(th)^2
g_phph <- function(th, a) (rp(a)^2 + a^2)^2*sin(th)^2 / (rp(a)^2 + a^2*cos(th)^2)
A_closed <- function(u, a) (rp(a)^2 + a^2) / (rp(a)^2 + a^2*u^2)      # u = cos th

cat("=== 1. the closed form is the metric anisotropy, checked directly\n\n")
ths <- c(0.15, 0.6, 1.1, pi/2, 2.2, 2.9)
for (a in c(0.3, 0.6, 0.9, 0.99)) {
  num <- sqrt(g_phph(ths,a)) / (sin(ths)*sqrt(g_thth(ths,a)))
  cat(sprintf("  a=%.2f   max|sqrt(g_pp)/(sin th sqrt(g_tt)) - A(cos th)| = %.3e\n",
              a, max(abs(num - A_closed(cos(ths), a)))))
}

cat("\n=== 2. Legendre content of A: even terms only\n\n")
# Legendre polynomials by recurrence; Gauss-Legendre nodes by Newton on P_n.
Pl <- function(n, x) { if (n == 0) return(rep(1, length(x))); if (n == 1) return(x)
  p0 <- rep(1, length(x)); p1 <- x
  for (k in 1:(n-1)) { p2 <- ((2*k+1)*x*p1 - k*p0)/(k+1); p0 <- p1; p1 <- p2 }
  p1 }
dPl <- function(n, x) n*(x*Pl(n,x) - Pl(n-1,x))/(x^2 - 1)
gauss <- function(n) { x <- cos(pi*(seq_len(n) - 0.25)/(n + 0.5))
  for (i in 1:100) x <- x - Pl(n,x)/dPl(n,x)
  list(x = x, w = 2/((1 - x^2)*dPl(n,x)^2)) }
G <- gauss(200)
coef <- function(n, a) sum(G$w * A_closed(G$x, a) * Pl(n, G$x)) * (2*n+1)/2

cat("      a       c_0       c_1       c_2       c_3       c_4    c_2/c_0\n")
for (a in c(0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)) {
  cs <- sapply(0:4, coef, a = a)
  cat(sprintf("  %5.2f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",
              a, cs[1], cs[2], cs[3], cs[4], cs[5], cs[3]/cs[1]))
}
cat("\n  every odd coefficient is zero to quadrature precision, so the rule is\n")
cat("  l -> l +/- 2 and never l +/- 1; c_2/c_0 is the strength and is 0 at a = 0.\n")

cat("\n=== 3. the mixing matrix <P_l|A|P_l'>, odd entries vanish\n\n")
M <- function(l, lp, a) sum(G$w * Pl(l,G$x) * A_closed(G$x,a) * Pl(lp,G$x))
a <- 0.9
cat(sprintf("  a = %.2f, entries <P_l|A|P_l'> for l,l' = 0..4\n\n", a))
tabm <- outer(0:4, 0:4, Vectorize(function(l, lp) M(l, lp, a)))
dimnames(tabm) <- list(paste0("l=",0:4), paste0("l'=",0:4))
print(round(tabm, 6))
odd <- max(abs(tabm[outer(0:4, 0:4, function(i,j) (i+j) %% 2 == 1)]))
cat(sprintf("\n  largest odd-parity entry = %.3e  (selection rule exact)\n", odd))

cat("\n=== 4. the same A controls the A.10 squashing\n\n")
for (a in c(0, 0.3, 0.5, 0.7, 0.9, 0.99))
  cat(sprintf("  a/M=%4.2f   g_thth(pole)/g_thth(eq) = %.6f   A(pi/2)/A(0) = %.6f   1+a^2/r_+^2 = %.6f\n",
      a, g_thth(0,a)/g_thth(pi/2,a), A_closed(0,a)/A_closed(1,a), 1 + a^2/rp(a)^2))
cat("\n  one function, two consequences: the closed-form kernel ratio is lost and\n")
cat("  the multipoles mix, both switched on by spin and both zero at a = 0.\n")
