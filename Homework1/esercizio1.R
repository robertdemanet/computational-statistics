# =============================================================================
# Homework 1 - Statistica Computazionale
# Esercizio 1: mistura Poisson/Gamma, stime Monte Carlo
#
#   Z ~ Bern(p)
#   Y|Z = 0 ~ Poisson(lambda)
#   Y|Z = 1 ~ Gamma(a, b)
# =============================================================================

set.seed(343240)

# ---- Parametri --------------------------------------------------------------

a <- runif(1, 2.5, 5.5)   # 4.706901
lambda <- runif(1, 2, 4)  # 3.049805
b <- 1
p <- 0.5
n <- 1000

cat("a =", a, "\n")
cat("lambda =", lambda, "\n")


# =============================================================================
# (a) Campioni dalla marginale di Y e stima MC della CDF
# =============================================================================

z <- rbinom(n, 1, p)

y_sample <- rep(NA, n)     # campioni dalla marginale di Y
y_sample_z0 <- rep(NA, n)  # campioni da Y|Z = 0 (-1 dove Z = 1)
y_sample_z1 <- rep(NA, n)  # campioni da Y|Z = 1 (-1 dove Z = 0)

for (i in 1:n)
{
  if (z[i] == 0)
  {
    y_sample[i] <- rpois(1, lambda)
    y_sample_z0[i] <- y_sample[i]
    y_sample_z1[i] <- -1
  }
  else
  {
    y_sample[i] <- rgamma(1, shape = a, rate = b)
    y_sample_z1[i] <- y_sample[i]
    y_sample_z0[i] <- -1
  }
}

# Osservazioni effettive delle due componenti (il valore -1 e' solo riempitivo)
y_z0 <- y_sample_z0[y_sample_z0 != -1]
y_z1 <- y_sample_z1[y_sample_z1 != -1]
n1 <- length(y_z0)  # numero di osservazioni da Y|Z = 0
n2 <- length(y_z1)  # numero di osservazioni da Y|Z = 1

# Stima MC della CDF: F(t) = E[I{Y <= t}] ~ media delle indicatrici
y_seq <- seq(min(y_sample), max(y_sample), length.out = n)
cdf <- rep(NA, n)

for (i in 1:n)
{
  cdf[i] <- sum(y_sample <= y_seq[i]) / n
}

plot(y_seq, cdf, type = "l", col = "blue", lwd = 2,
     main = "Cumulative Distribution of Y",
     xlab = "y_seq", ylab = "cdf")


# =============================================================================
# (b) Stima MC di P(Y in [3.5, 4.5]) e P(Y in [3.5, 4.5] | Z = 0)
#     con varianza degli stimatori
# =============================================================================

ind_y <- (y_sample >= 3.5 & y_sample <= 4.5)
ind_y_z0 <- (y_z0 >= 3.5 & y_z0 <= 4.5)

p_y <- sum(ind_y) / n
p_y_z0 <- sum(ind_y_z0) / n1

cat("P(Y in [3.5,4.5]) =", p_y, "\n")
cat("P(Y in [3.5,4.5] | Z = 0) =", p_y_z0, "\n")

# Varianza dello stimatore MC: v_n = (1/n^2) * somma degli scarti al quadrato
v_n_y <- sum((ind_y - p_y)^2) / n^2
v_n_yz0 <- sum(((y_sample_z0 >= 3.5 & y_sample_z0 <= 4.5) - p_y_z0)^2) / n1^2

cat("var stimatore su Y =", v_n_y, "\n")
cat("var stimatore su Y|Z=0 =", v_n_yz0, "\n")


# =============================================================================
# (c) Stima MC di P(Z = 0 | Y in [1.5, 3.5]) e P(Z = 1 | Y in [1.5, 3.5])
#     tramite il teorema di Bayes
# =============================================================================

p_A <- sum(y_sample >= 1.5 & y_sample <= 3.5) / n       # P(Y in A)
p_A_z0 <- sum(y_z0 >= 1.5 & y_z0 <= 3.5) / n1           # P(Y in A | Z = 0)
p_A_z1 <- sum(y_z1 >= 1.5 & y_z1 <= 3.5) / n2           # P(Y in A | Z = 1)

p_z0 <- n1 / n  # stima MC di P(Z = 0)
p_z1 <- n2 / n  # stima MC di P(Z = 1)

post_z0 <- p_A_z0 * p_z0 / p_A
post_z1 <- p_A_z1 * p_z1 / p_A

cat("P(Z = 0 | Y in [1.5,3.5]) =", post_z0, "\n")
cat("P(Z = 1 | Y in [1.5,3.5]) =", post_z1, "\n")


# =============================================================================
# (d) Quantili empirici tramite l'inversa generalizzata
#     F^-1(u) = inf{ x : u <= F(x) }
# =============================================================================

u_levels <- c(0.1, 0.2, 0.5, 0.75)

# --- Quantili di Y: si usa la CDF stimata al punto (a) ---
quantili_y <- rep(NA, length(u_levels))

for (k in 1:length(u_levels))
{
  idx <- which(u_levels[k] <= cdf)[1]
  quantili_y[k] <- y_seq[idx]
}

cat("Quantili empirici di Y:\n")
print(data.frame(livello = u_levels, valore = quantili_y))

# --- Quantili di Z: la CDF viene costruita manualmente ---
z_seq <- seq(min(z), max(z), length.out = n)
cdf_z <- rep(NA, n)

for (i in 1:n)
{
  cdf_z[i] <- sum(z <= z_seq[i]) / n
}

plot(z_seq, cdf_z, type = "l", col = "blue", lwd = 2,
     main = "Cumulative Distribution of Z",
     xlab = "z_seq", ylab = "cdf")

quantili_z <- rep(NA, length(u_levels))

for (k in 1:length(u_levels))
{
  idx <- which(u_levels[k] <= cdf_z)[1]
  quantili_z[k] <- z_seq[idx]
}

cat("Quantili empirici di Z:\n")
print(data.frame(livello = u_levels, valore = quantili_z))


# =============================================================================
# (e) Stima MC della parte discreta (Z = 0) e della parte continua (Z = 1)
# =============================================================================

x <- seq(0, 10, by = 0.25)  # 41 valori su cui stimare le due parti
tol <- 0.5                  # tolleranza necessaria per la parte continua

stima_y_z0 <- rep(NA, length(x))
stima_y_z1 <- rep(NA, length(x))

for (i in 1:length(x))
{
  # parte discreta: frequenza relativa del valore esatto
  stima_y_z0[i] <- sum(y_z0 == x[i]) / n1

  # parte continua: frequenza relativa nell'intorno (x - tol, x + tol)
  stima_y_z1[i] <- sum(y_z1 > x[i] - tol & y_z1 < x[i] + tol) / n2
}

# --- Parte discreta vs Poisson teorica ---
plot(x, stima_y_z0, type = "l", col = "red", lwd = 2,
     ylim = c(0, max(stima_y_z0, dpois(round(x), lambda))),
     main = "Empirical vs Theoretical Distribution of Y | Z = 0",
     xlab = "x values", ylab = "Density")
lines(x, dpois(round(x), lambda), col = "green", lwd = 2)
legend("topright", legend = c("Empirical", "Theoretical (Poisson)"),
       col = c("red", "green"), lwd = 2)

# --- Parte continua vs Gamma teorica ---
plot(x, stima_y_z1, type = "l", col = "blue", lwd = 2,
     ylim = c(0, max(stima_y_z1, dgamma(x, shape = a, rate = b))),
     main = "Empirical vs Theoretical Distribution of Y | Z = 1",
     xlab = "x values", ylab = "Density")
lines(x, dgamma(x, shape = a, rate = b), col = "green", lwd = 2)
legend("topright", legend = c("Empirical", "Theoretical (Gamma)"),
       col = c("blue", "green"), lwd = 2)
