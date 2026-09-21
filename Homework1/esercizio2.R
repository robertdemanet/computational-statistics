# =============================================================================
# Homework 1 - Statistica Computazionale
# Esercizio 2: distribuzione logistic-normal
#   - simulazione con metodo Accept-Reject
#   - derivazione e campionamento della a-posteriori di mu
#   - confronto prior/posterior
#   - stima della media con Importance Sampling
# =============================================================================

set.seed(343240)

# ---- Parametri --------------------------------------------------------------

mu <- runif(1, -1.5, 1.5)      # 0.7069013
sigma2 <- runif(1, 0.5, 1.5)   # 1.0249025

cat("mu =", mu, "\n")
cat("sigma2 =", sigma2, "\n")

# ---- Densita' della logistic-normal -----------------------------------------

logit <- function(x)
{
  log(x / (1 - x))
}

f_ln <- function(x, mu, sigma2)
{
  1 / (x * (1 - x)) * 1 / sqrt(2 * pi * sigma2) *
    exp(-(logit(x) - mu)^2 / (2 * sigma2))
}


# =============================================================================
# (a) Simulazione dei dati con il metodo Accept-Reject
# =============================================================================

# Massimo della densita' sul suo dominio (0, 1)
x_grid <- seq(0.0001, 0.9999, length.out = 10000)
M <- max(f_ln(x_grid, mu, sigma2))
cat("M =", M, "\n")

n_star <- 1000  # numero di proposte

# Ascisse dei campioni proposti e corrispondenti ordinate
Y <- runif(n_star, 0, 1)
U <- runif(n_star, 0, M)

# Si accetta la coppia (x*, u) se u < f(x*)
accepted <- U < f_ln(Y, mu, sigma2)

perc_acceptance <- sum(accepted) / n_star
cat("perc_acceptance =", perc_acceptance, "\n")

X <- Y[accepted]    # ascisse accettate
U_X <- U[accepted]  # ordinate accettate

# La richiesta e' di n = 10 dati: si prendono i primi 10 accettati
n_obs <- 10
X10_sample <- X[1:n_obs]
U_X10_sample <- U_X[1:n_obs]

plot(x_grid, f_ln(x_grid, mu, sigma2), type = "l", col = "blue", lwd = 2,
     ylim = c(0, M),
     main = "Density with Accepted Samples (N = 10)",
     xlab = "X values", ylab = "Density")
points(Y[!accepted][1:n_obs], U[!accepted][1:n_obs], col = "black", pch = 19)
points(X10_sample, U_X10_sample, col = "red", pch = 19)
lines(density(X10_sample, from = 0, to = 1), col = "green", lwd = 2)
legend("topleft", legend = c("Density", "Accepted Samples", "Empiric Density"),
       col = c("blue", "red", "green"), lwd = 2)


# =============================================================================
# (b) A-posteriori di mu con prior mu ~ N(0, 100) e sigma2 noto
#
#     pi(mu | x_1,...,x_n) ~ N( sum(logit(x_i)) / (n + sigma2/100),
#                               sigma2 / (n + sigma2/100) )
# =============================================================================

posterior_par <- function(x, sigma2)
{
  n_x <- length(x)
  media <- sum(logit(x)) / (n_x + sigma2 / 100)
  varianza <- sigma2 / (n_x + sigma2 / 100)

  return(c(media, varianza))
}

par10 <- posterior_par(X10_sample, sigma2)
cat("posterior (N = 10): media =", par10[1], " varianza =", par10[2], "\n")

n_post <- 1000
post10_samples <- rnorm(n_post, par10[1], sqrt(par10[2]))

hist(post10_samples, breaks = 30, freq = FALSE, col = "lightblue",
     main = "Posterior Distribution of Mu (N = 10)",
     xlab = "Mu values", ylab = "Density")
lines(density(post10_samples), col = "red", lwd = 2)


# =============================================================================
# (c) Confronto tra prior e posterior al crescere delle osservazioni
# =============================================================================

n_obs2 <- 100
X100_sample <- X[1:n_obs2]
U_X100_sample <- U_X[1:n_obs2]

# --- Confronto tra le due simulazioni accept-reject ---
par(mfrow = c(1, 2))

plot(x_grid, f_ln(x_grid, mu, sigma2), type = "l", col = "blue", lwd = 2,
     ylim = c(0, M), main = "Accepted Samples (N = 10)",
     xlab = "X values", ylab = "Density")
points(X10_sample, U_X10_sample, col = "red", pch = 19)
lines(density(X10_sample, from = 0, to = 1), col = "green", lwd = 2)

plot(x_grid, f_ln(x_grid, mu, sigma2), type = "l", col = "blue", lwd = 2,
     ylim = c(0, M), main = "Accepted Samples (N = 100)",
     xlab = "X values", ylab = "Density")
points(X100_sample, U_X100_sample, col = "red", pch = 19)
lines(density(X100_sample, from = 0, to = 1), col = "green", lwd = 2)

par(mfrow = c(1, 1))

# --- Prior, posterior con 10 osservazioni e posterior con 100 osservazioni ---
par100 <- posterior_par(X100_sample, sigma2)
cat("posterior (N = 100): media =", par100[1], " varianza =", par100[2], "\n")

prior_samples <- rnorm(n_post, 0, 10)
post100_samples <- rnorm(n_post, par100[1], sqrt(par100[2]))

par(mfrow = c(1, 3))

plot(density(prior_samples), col = "purple", lwd = 2, xlim = c(-5, 5),
     main = "Prior Distribution", xlab = "Mu values", ylab = "Density")

plot(density(post10_samples), col = "orange", lwd = 2, xlim = c(-5, 5),
     main = "Posterior (N = 10)", xlab = "Mu values", ylab = "Density")
abline(v = mu, col = "red", lwd = 2)

plot(density(post100_samples), col = "green", lwd = 2, xlim = c(-5, 5),
     main = "Posterior (N = 100)", xlab = "Mu values", ylab = "Density")
abline(v = mu, col = "red", lwd = 2)

par(mfrow = c(1, 1))


# =============================================================================
# (d) Stima di E(X) con il metodo Importance Sampling
#
#     E(X) = int x f(x) dx = E_g[ X f(X) / g(X) ],  g = densita' Beta(1.2, 1.2)
# =============================================================================

n_is <- 1000
a_beta <- 1.2
b_beta <- 1.2

x_beta <- rbeta(n_is, a_beta, b_beta)

weights <- x_beta * f_ln(x_beta, mu, sigma2) / dbeta(x_beta, a_beta, b_beta)
E_X <- sum(weights) / n_is

cat("E(X) stimato con Importance Sampling =", E_X, "\n")

plot(x_grid, f_ln(x_grid, mu, sigma2), type = "l", col = "purple", lwd = 2,
     main = "Logistic-Normal Density and Importance Sampling Mean",
     xlab = "X values", ylab = "Density")
abline(v = E_X, col = "red", lwd = 2)
legend("topleft", legend = c("Density", "E(X)"),
       col = c("purple", "red"), lwd = 2)
