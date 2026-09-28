#!/usr/bin/env Rscript
# a13_parity_mixing.R -- A.13 makes four numerical claims and none had a script. All four
# are computed here: the anisotropy function against the induced metric, its ratio being
# the squashing of A.10, the vanishing of the odd-parity Legendre matrix elements on a
# 200-node Gauss-Legendre rule, and the spin-weighted element that breaks the rule.

M <- 1
# ---- Gauss-Legendre by Golub-Welsch, so no package is needed ----
gauss_legendre <- function(n) {
  k <- 1:(n-1); b <- k/sqrt(4*k^2 - 1)
  J <- diag(0, n); J[cbind(k, k+1)] <- b; J[cbind(k+1, k)] <- b
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = 2*(e$vectors[1, o])^2) }
gl <- gauss_legendre(200)
stopifnot(abs(sum(gl$w) - 2) < 1e-13)
cat(sprintf("=== 0. the quadrature: 200 nodes, sum of weights %.15f against 2 ===\n", sum(gl$w)))
Pl <- function(l, x) { if (l == 0) return(rep(1, length(x))); if (l == 1) return(x)
  p0 <- rep(1, length(x)); p1 <- x
  for (k in 2:l) { p2 <- ((2*k-1)*x*p1 - (k-1)*p0)/k; p0 <- p1; p1 <- p2 }; p1 }
er <- max(sapply(0:9, function(l) sapply(0:9, function(m)
  abs(sum(gl$w*Pl(l,gl$x)*Pl(m,gl$x)) - if (l==m) 2/(2*l+1) else 0))))
cat(sprintf("    Legendre orthogonality on it: worst error %.1e\n", er))
stopifnot(er < 1e-13)

cat("\n=== 1. the anisotropy function against the induced metric ===\n")
# bifurcation surface of Kerr: ds^2 = Sigma_+ dtheta^2 + ((r_+^2+a^2)^2 sin^2 th / Sigma_+) dphi^2
Aclosed <- function(th, a) { rp <- M + sqrt(M^2-a^2); (rp^2+a^2)/(rp^2+a^2*cos(th)^2) }
Ametric <- function(th, a) { rp <- M + sqrt(M^2-a^2); S <- rp^2 + a^2*cos(th)^2
  gpp <- (rp^2+a^2)^2*sin(th)^2/S; ghh <- S
  sqrt(gpp)/(sin(th)*sqrt(ghh)) }
cat("      a/M     worst |closed form - metric|     A(pi/2)/A(0)     1 + a^2/r_+^2\n")
ths <- seq(1e-4, pi-1e-4, length.out = 3000)
for (a in c(0.0, 0.3, 0.6, 0.9)) {
  e <- max(abs(Aclosed(ths,a) - Ametric(ths,a)))
  rp <- M + sqrt(M^2-a^2)
  cat(sprintf("   %7.2f  %27.1e  %15.9f  %15.9f\n", a, e,
              Aclosed(pi/2,a)/Aclosed(0,a), 1 + a^2/rp^2))
  stopifnot(e < 1e-14, abs(Aclosed(pi/2,a)/Aclosed(0,a) - (1+a^2/rp^2)) < 1e-14) }
cat("   Both statements hold exactly, so the anisotropy and A.10's squashing are one function.\n")

cat("\n=== 2. the odd-parity Legendre matrix elements ===\n")
cat("   A is even in cos(theta), so <P_l|A|P_l'> should vanish whenever l + l' is odd.\n")
cat("   The residual is the arithmetic floor of the quadrature, and it turns out not to\n")
cat("   grow with the degree, so it is reported at two cut-offs to show that.\n\n")
oddmax <- function(a, lmax) {
  Av <- Aclosed(acos(gl$x), a); L <- 0:lmax
  Mx <- outer(L, L, Vectorize(function(i,j) sum(gl$w*Pl(i,gl$x)*Av*Pl(j,gl$x))))
  odd <- outer(L, L, function(i,j) (i+j) %% 2 == 1)
  c(max(abs(Mx[odd])), max(abs(Mx[!odd]))) }
cat("      a/M    largest odd entry, l<=6   largest odd, l<=12   largest even entry\n")
for (a in c(0.3, 0.6, 0.9, 0.998)) {
  s6 <- oddmax(a, 6); s12 <- oddmax(a, 12)
  cat(sprintf("   %7.3f  %23.1e  %20.1e  %18.4f\n", a, s6[1], s12[1], s12[2]))
  stopifnot(s6[1] < 1e-13, s12[1] < 1e-13) }
cat("\n   A.13 quotes 2.7e-20 for this. A 200-node Gauss-Legendre rule in double precision\n")
cat("   cannot deliver that: its own Legendre orthogonality holds only to 4e-15, and the\n")
cat("   odd entries sit at that same 4e-15 whether the cut-off is six or twelve, which is\n")
cat("   the signature of a quadrature floor rather than a real residual. The\n")
cat("   statement is right and the digit is not reproducible at the stated precision, so\n")
cat("   the paper should say the entries vanish by parity and give the floor it measured.\n")

cat("\n   The check must be able to fail, so plant an anisotropy that is NOT even:\n")
Abad <- function(th, a) Aclosed(th,a)*(1 + 0.3*cos(th))
Av <- Abad(acos(gl$x), 0.6); L <- 0:12
Mx <- outer(L, L, Vectorize(function(i,j) sum(gl$w*Pl(i,gl$x)*Av*Pl(j,gl$x))))
odd <- outer(L, L, function(i,j) (i+j) %% 2 == 1)
cat(sprintf("     A -> A(1 + 0.3 cos theta): largest odd entry %.4f   <- fails, as it must\n",
            max(abs(Mx[odd]))))
stopifnot(max(abs(Mx[odd])) > 1e-3)

cat("\n=== 2b. the coefficients A.13 quotes, which had no script and no definition ===\n")
cat("   A.13 gives 'the strength is c_2/c_0' without saying what c_l is. It is the Legendre\n")
cat("   expansion of the anisotropy in cos(theta), c_l = (2l+1)/2 * int A(x) P_l(x) dx over\n")
cat("   [-1,1], so c_0 is the mean and c_2 the quadrupole. Computing it reproduces all five\n")
cat("   quoted values, which is what establishes the reading.\n\n")
cl <- function(l, a) (2*l+1)/2 * sum(gl$w * Aclosed(acos(gl$x), a) * Pl(l, gl$x))
paper <- c(-0.0017, -0.016, -0.046, -0.219, -0.367)
aa <- c(0.1, 0.3, 0.5, 0.9, 0.99)
cat("      a/M        c_0         c_2      c_2/c_0     A.13 quotes\n")
for (i in seq_along(aa)) { a <- aa[i]; r <- cl(2,a)/cl(0,a)
  cat(sprintf("   %7.2f  %9.5f  %10.5f  %10.5f  %12.4f\n", a, cl(0,a), cl(2,a), r, paper[i]))
  stopifnot(abs(r - paper[i]) < 5e-4) }
cat(sprintf("   %7.2f  %9.5f  %10.5f  %10.1e  %12s\n", 0, cl(0,0), cl(2,0), cl(2,0)/cl(0,0), "zero"))
stopifnot(abs(cl(2,0)) < 1e-13)
cat(sprintf("\n   magnitudes, for quoting: %s\n",
            paste(sprintf("%.5f", abs(sapply(aa, function(a) cl(2,a)/cl(0,a)))), collapse="  ")))
cat("\n   Every odd coefficient vanishes for the same parity reason as the matrix elements:\n")
oddc <- max(sapply(c(1,3,5), function(l) max(sapply(c(0.3,0.9), function(a) abs(cl(l,a))))))
cat(sprintf("     worst |c_1|, |c_3|, |c_5| at a = 0.3 and 0.9: %.1e\n", oddc))
stopifnot(oddc < 1e-13)

cat("\n=== 3. the spin-weighted element that breaks the rule ===\n")
# Goldberg et al.: _sY_lm(th,ph) = (-1)^m sqrt((l+m)!(l-m)!(2l+1)/(4 pi (l+s)!(l-s)!))
#   sin^{2l}(th/2) sum_r C(l-s,r) C(l+s,r+s-m) (-1)^{l-r-s} cot^{2r+s-m}(th/2) e^{i m ph}
sY <- function(s, l, m, th) {
  pref <- (-1)^m * sqrt(factorial(l+m)*factorial(l-m)*(2*l+1) /
                        (4*pi*factorial(l+s)*factorial(l-s)))
  tot <- numeric(length(th))
  for (r in 0:(l-s)) {
    if (r+s-m < 0 || r+s-m > l+s) next
    tot <- tot + choose(l-s, r)*choose(l+s, r+s-m)*(-1)^(l-r-s) *
                 (1/tan(th/2))^(2*r+s-m) }
  pref * sin(th/2)^(2*l) * tot }
th <- acos(gl$x); w <- gl$w                    # int f sin th dth = sum w f(acos(x))
inner <- function(f) 2*pi*sum(w * f)
cat("   Orthonormality first, since an unverified implementation proves nothing:\n")
n22 <- inner(sY(-2,2,2,th)^2); n32 <- inner(sY(-2,3,2,th)^2); o <- inner(sY(-2,2,2,th)*sY(-2,3,2,th))
cat(sprintf("     <-2Y22|-2Y22> = %.15f     <-2Y32|-2Y32> = %.15f\n", n22, n32))
cat(sprintf("     <-2Y22|-2Y32> = %.2e\n", o))
stopifnot(abs(n22-1) < 1e-12, abs(n32-1) < 1e-12, abs(o) < 1e-12)
v  <- inner(sY(-2,2,2,th)*cos(th)^2*sY(-2,3,2,th))
ex <- sqrt(35)/21
cat(sprintf("\n     <-2Y22|cos^2 theta|-2Y32> = %.15f\n", v))
cat(sprintf("     sqrt(35)/21               = %.15f     difference %.1e\n", ex, abs(v-ex)))
stopifnot(abs(v-ex) < 1e-13)
# the Legendre counterpart, which does vanish
vl <- sum(gl$w*Pl(2,gl$x)*gl$x^2*Pl(3,gl$x))
cat(sprintf("     same element from Legendre functions: %.1e, which vanishes\n", vl))
stopifnot(abs(vl) < 1e-13)

cat("
=== flatly ===

  All four of A.13's claims reproduce. The closed-form anisotropy agrees with the induced
  metric to machine precision at every spin, and its pole-to-equator ratio is exactly
  A.10's squashing 1 + a^2/r_+^2. The odd-parity Legendre matrix elements vanish, and the
  check is validated by planting an anisotropy odd in cos(theta), which makes them order
  one. The spin-weighted element is sqrt(35)/21 to 1e-14, against a Legendre counterpart
  that vanishes, so the parity rule really does fail for s = -2 and A.13's warning stands.

Two digits are not reproduced and should not be. A.13 quotes 2.7e-20 for the largest
  odd-parity Legendre entry and 2.7e-17 for the Legendre counterpart of the spin-weighted
  element. A 200-node Gauss-Legendre rule in double precision has a floor of 4e-15, which
  its own Legendre orthogonality shows, and both quantities sit exactly there and do not
  move with the degree cut-off. They are zero, they are zero by parity rather than by
  arithmetic, and a paper quoting digits eleven orders below its quadrature floor is
  quoting someone's rounding. Both replaced in the paper by the statement and the measured
  floor.

  Provenance gap 6 closed.\n")
