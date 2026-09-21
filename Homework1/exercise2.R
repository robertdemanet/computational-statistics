# =============================================================================
# Homework 1 - Computational Statistics
# Exercise 2: logistic-normal distribution
#   - data simulation through the Accept-Reject method
#   - derivation of and sampling from the posterior of mu
#   - prior/posterior comparison
#   - estimation of the mean through Importance Sampling
# =============================================================================

set.seed(343240)

# ---- Parameters -------------------------------------------------------------

mu <- runif(1, -1.5, 1.5)      # 0.7069013
sigma2 <- runif(1, 0.5, 1.5)   # 1.0249025

cat("mu =", mu, "\n")
cat("sigma2 =", sigma2, "\n")

# ---- Logistic-normal density ------------------------------------------------

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
# (a) Data simulation through the Accept-Reject method
# =============================================================================

# Maximum of the density over its support (0, 1)
x_grid <- seq(0.0001, 0.9999, length.out = 10000)
M <- max(f_ln(x_grid, mu, sigma2))
cat("M =", M, "\n")

n_star <- 1000  # number of proposals

# Abscissae of the proposed samples and corresponding ordinates
Y <- runif(n_star, 0, 1)
U <- runif(n_star, 0, M)

# The pair (x*, u) is accepted if u < f(x*)
accepted <- U < f_ln(Y, mu, sigma2)

perc_acceptance <- sum(accepted) / n_star
cat("perc_acceptance =", perc_acceptance, "\n")

X <- Y[accepted]    # accepted abscissae
U_X <- U[accepted]  # accepted ordinates

# The exercise asks for n = 10 data points: the first 10 accepted ones are kept
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
# (b) Posterior of mu under the prior mu ~ N(0, 100) with sigma2 known
#
#     pi(mu | x_1,...,x_n) ~ N( sum(logit(x_i)) / (n + sigma2/100),
#                               sigma2 / (n + sigma2/100) )
# =============================================================================

posterior_par <- function(x, sigma2)
{
  n_x <- length(x)
  post_mean <- sum(logit(x)) / (n_x + sigma2 / 100)
  post_var <- sigma2 / (n_x + sigma2 / 100)

  return(c(post_mean, post_var))
}

par10 <- posterior_par(X10_sample, sigma2)
cat("posterior (N = 10): mean =", par10[1], " variance =", par10[2], "\n")

n_post <- 1000
post10_samples <- rnorm(n_post, par10[1], sqrt(par10[2]))

hist(post10_samples, breaks = 30, freq = FALSE, col = "lightblue",
     main = "Posterior Distribution of Mu (N = 10)",
     xlab = "Mu values", ylab = "Density")
lines(density(post10_samples), col = "red", lwd = 2)


# =============================================================================
# (c) Prior against posterior as the number of observations grows
# =============================================================================

n_obs2 <- 100
X100_sample <- X[1:n_obs2]
U_X100_sample <- U_X[1:n_obs2]

# --- Comparison between the two accept-reject simulations ---
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

# --- Prior, posterior with 10 observations and posterior with 100 observations ---
par100 <- posterior_par(X100_sample, sigma2)
cat("posterior (N = 100): mean =", par100[1], " variance =", par100[2], "\n")

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
# (d) Estimation of E(X) through Importance Sampling
#
#     E(X) = int x f(x) dx = E_g[ X f(X) / g(X) ],  g = Beta(1.2, 1.2) density
# =============================================================================

n_is <- 1000
a_beta <- 1.2
b_beta <- 1.2

x_beta <- rbeta(n_is, a_beta, b_beta)

weights <- x_beta * f_ln(x_beta, mu, sigma2) / dbeta(x_beta, a_beta, b_beta)
E_X <- sum(weights) / n_is

cat("E(X) estimated through Importance Sampling =", E_X, "\n")

plot(x_grid, f_ln(x_grid, mu, sigma2), type = "l", col = "purple", lwd = 2,
     main = "Logistic-Normal Density and Importance Sampling Mean",
     xlab = "X values", ylab = "Density")
abline(v = E_X, col = "red", lwd = 2)
legend("topleft", legend = c("Density", "E(X)"),
       col = c("purple", "red"), lwd = 2)
