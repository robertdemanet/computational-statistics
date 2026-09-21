# =============================================================================
# Homework 2 - Computational Statistics
# Exercise 1: hierarchical spatial model
#
#   Y(s) | W(s) ~ GP(W(s), tau2)
#   W(s)        ~ GP(m(s), C(||s - s'||; phi, sigma2))
#   m(s) = beta0 + beta1 * s1
#   C(||s - s'||; phi, sigma2) = sigma2 * exp(-phi * ||s - s'||)
# =============================================================================

library(MASS)
library(mvtnorm)
library(coda)
library(ggplot2)

set.seed(343240)


# =============================================================================
# (1) Simulation of 100 observations from the model
# =============================================================================

n <- 100
s1 <- runif(n, 0, 10)  # x coordinates
s2 <- runif(n, 0, 10)  # y coordinates

points <- cbind(s1, s2)  # matrix holding the coordinates of the points (x, y)

# Model parameters
tau2 <- 0.5
sigma2 <- 0.5
phi <- 3 / 10
beta0 <- 2
beta1 <- 0.1

# Mean and covariance matrix of W(s)
m_s <- beta0 + s1 * beta1
distance_matrix <- as.matrix(dist(points))  # distance matrix
covariance_matrix <- sigma2 * exp(-phi * distance_matrix)

# Simulation of the Gaussian process W(s)
W_s <- mvrnorm(1, mu = m_s, Sigma = covariance_matrix)

# Y(s) = W(s) + epsilon,  epsilon ~ N(0, tau2)
epsilon <- rnorm(n, 0, sqrt(tau2))
Y_s <- W_s + epsilon

data_sim <- data.frame(x = s1, y = s2, W_s = W_s, Y_s = Y_s)
print(head(data_sim, 3))

# Spatial representation of W(s) and Y(s)
print(
  ggplot(data_sim, aes(x = x, y = y, color = W_s)) +
    geom_point(size = 3) +
    scale_color_gradientn(colours = c("blue", "green", "yellow", "red")) +
    labs(title = "Simulation of the Gaussian process W(s)",
         x = "x coordinate", y = "y coordinate", color = "W(s)")
)

print(
  ggplot(data_sim, aes(x = x, y = y, color = Y_s)) +
    geom_point(size = 3) +
    scale_color_gradientn(colours = c("blue", "green", "yellow", "red")) +
    labs(title = "Simulation of the Gaussian process Y(s)",
         x = "x coordinate", y = "y coordinate", color = "Y(s)")
)


# =============================================================================
# (2) Scatterplot of the coordinates with 4 groups defined by the quartiles of y
# =============================================================================

quantiles_y <- quantile(data_sim$Y_s, probs = c(0, 0.25, 0.5, 0.75, 1))

data_sim$group <- cut(data_sim$Y_s,
                      breaks = quantiles_y,
                      include.lowest = FALSE,
                      include.highest = TRUE,
                      labels = c('0% <= q_i <= 25%', '25% < q_i <= 50%',
                                 '50% < q_i <= 75%', '75% < q_i <= 100%'))

print(head(data_sim, 3))

print(
  ggplot(data_sim, aes(x = x, y = y, color = group)) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c("blue", "green", "orange", "red")) +
    labs(title = "Scatterplot of the coordinates with 4 groups",
         x = "x coordinate", y = "y coordinate", color = "Quantile group")
)


# =============================================================================
# (3) MCMC algorithm for f(w, beta0, beta1, tau2, sigma2, phi | y_o)
#
#     Priors: beta ~ N_2(0, 100 * I_2),  sigma2 ~ IG(1,1),  tau2 ~ IG(1,1),
#             phi ~ Gamma(1,1),  w ~ N(X beta, C(phi, sigma2))
# =============================================================================

# --- Random selection of the observations ------------------------------------

n_new <- runif(1, 10, 90)
n_new <- floor(n_new)

indices <- sample(1:100, n_new)
indices <- sort(indices)

y_o <- data_sim$Y_s[indices]
D_o <- data_sim[indices, c('x', 'y')]

y_u <- data_sim$Y_s[-indices]
D_u <- data_sim[-indices, c('x', 'y')]

cat("n_new =", n_new, "\n")

# --- Initialisation of the algorithm -----------------------------------------

n_iter <- 60000
burn_in <- 20000

beta0_samples <- rep(NA, n_iter)
beta0_samples[1] <- 0
beta1_samples <- rep(NA, n_iter)
beta1_samples[1] <- 0
tau2_samples <- rep(NA, n_iter)
tau2_samples[1] <- 1
sigma2_samples <- rep(NA, n_iter)
sigma2_samples[1] <- 1
phi_samples <- rep(NA, n_iter)
phi_samples[1] <- 0.1
W_samples <- matrix(NA, nrow = n_iter, ncol = n_new)
W_samples[1, ] <- rep(0, n_new)

# x coordinates and distance matrix of the selected observations only
D_o_x <- data_sim[indices, c('x')]
distance_matrix_o <- as.matrix(dist(D_o))

# Design matrix
ones <- rep(1, n_new)
X <- cbind(ones, D_o_x)

# Prior precision matrix of beta
diag_matrix <- diag(1 / 100, 2)

# --- Full conditional of phi (up to a constant) ------------------------------

log_fc_phi <- function(phi, w, beta, sigma2, X, C)
{
  # Determinant
  det_C <- determinant(C, logarithm = TRUE)$modulus  # log-determinant

  det_C <- abs(det_C)  # abs is used in case the determinant is negative

  # Log-likelihood
  res <- w - X %*% beta
  log_likelihood <- -0.5 * (det_C + t(res) %*% solve(C) %*% res)

  # Log-prior of phi, Gamma distributed
  log_prior <- dgamma(phi, 1, 1, log = TRUE)

  return(log_likelihood + log_prior)
}

# --- Settings of the adaptive Metropolis step for phi ------------------------

sd_proposal <- 0.1
alpha_phi <- 0
nbatch <- 50
A <- 100
B <- 1000
alpha_target <- 0.234

# Bounds for phi derived from the spatial distances
phi_lower <- 3 / max(dist(points))
phi_upper <- 3 / min(dist(points))

# --- MCMC loop ---------------------------------------------------------------

for (i in 2:n_iter)
{
  # Current covariance matrix and its inverse
  C_o <- sigma2_samples[i - 1] * exp(-phi_samples[i - 1] * distance_matrix_o)
  inv_C_o <- solve(C_o)

  ## Step (1): update of beta
  # Vp
  Vp <- solve(t(X) %*% inv_C_o %*% X + diag_matrix)

  # Mp
  Mp <- Vp %*% (t(X) %*% inv_C_o %*% W_samples[i - 1, ])

  # Sample beta
  beta <- mvrnorm(1, Mp, Vp)
  beta0_samples[i] <- beta[1]
  beta1_samples[i] <- beta[2]

  ## Step (2): update of sigma2
  a_sigma2 <- n_new / 2 + 1
  quad_form <- t(W_samples[i - 1, ] - X %*% beta) %*% inv_C_o %*%
    (W_samples[i - 1, ] - X %*% beta)
  b_sigma2 <- quad_form / 2 + 1

  # Sample sigma2
  sigma2_samples[i] <- 1 / rgamma(1, a_sigma2, b_sigma2)

  ## Step (3): update of tau2
  a_tau2 <- n_new / 2 + 1
  b_tau2 <- 0.5 * t(y_o - W_samples[i - 1, ]) %*% (y_o - W_samples[i - 1, ]) + 1

  # Sample tau2
  tau2_samples[i] <- 1 / rgamma(1, a_tau2, b_tau2)

  ## Step (4): update of phi (Metropolis)
  phi <- phi_samples[i - 1]
  phi_prop <- rnorm(1, phi_samples[i - 1], sd_proposal)

  if (phi_prop >= phi_lower && phi_prop <= phi_upper)
  {
    # The covariances must be recomputed at the value of phi they refer to
    C_prop <- sigma2_samples[i] * exp(-phi_prop * distance_matrix_o)
    C_curr <- sigma2_samples[i] * exp(-phi_samples[i - 1] * distance_matrix_o)

    log_ratio <- log_fc_phi(phi_prop, W_samples[i - 1, ], beta,
                            sigma2_samples[i], X, C_prop) -
      log_fc_phi(phi_samples[i - 1], W_samples[i - 1, ], beta,
                 sigma2_samples[i], X, C_curr)

    alpha_phi <- alpha_phi + min(exp(log_ratio), 1)

    if (runif(1, 0, 1) < exp(log_ratio))
    {
      phi <- phi_prop
    }
  }

  # Adaptation of the proposal sd towards the target acceptance rate
  if (i %% nbatch == 0)
  {
    alpha_phi <- alpha_phi / nbatch
    sd_proposal <- exp(log(sd_proposal) + A / (B + i) * (alpha_phi - alpha_target))
    alpha_phi <- 0
  }

  phi_samples[i] <- phi

  ## Step (5): update of W
  Q <- inv_C_o + diag(1 / tau2_samples[i], n_new)
  b <- inv_C_o %*% X %*% beta + y_o / tau2_samples[i]
  muW_post <- solve(Q) %*% b
  sigmaW_post <- solve(Q)

  # Sample W
  W_samples[i, ] <- mvtnorm::rmvnorm(1, as.vector(muW_post), sigmaW_post)
}

# --- Traces and effective sample size ----------------------------------------

keep <- (burn_in + 1):n_iter

plot(beta0_samples[keep], type = "l", main = "beta0", xlab = "Iterations", ylab = "beta0")
plot(beta1_samples[keep], type = "l", main = "beta1", xlab = "Iterations", ylab = "beta1")
plot(tau2_samples[keep], type = "l", main = "tau2", xlab = "Iterations", ylab = "tau2")
plot(sigma2_samples[keep], type = "l", main = "sigma2", xlab = "Iterations", ylab = "sigma2")
plot(phi_samples[keep], type = "l", main = "phi", xlab = "Iterations", ylab = "phi")

eff_samples <- data.frame(
  parameter = c("beta0", "beta1", "tau2", "sigma2", "phi"),
  effective_samples = c(effectiveSize(beta0_samples[keep]),
                        effectiveSize(beta1_samples[keep]),
                        effectiveSize(tau2_samples[keep]),
                        effectiveSize(sigma2_samples[keep]),
                        effectiveSize(phi_samples[keep]))
)
print(eff_samples)


# =============================================================================
# (4) Posterior of epsilon(s) = Y(s) - (beta0 + beta1 * s1),  s in D_o
# =============================================================================

epsilon_samples <- matrix(NA, nrow = n_iter - burn_in, ncol = n_new)

for (i in (burn_in + 1):n_iter)
{
  m_post <- beta0_samples[i] + beta1_samples[i] * D_o$x
  epsilon_samples[i - burn_in, ] <- y_o - m_post
}

epsilon_mean <- apply(epsilon_samples, 2, mean)
D_o$epsilon_mean <- epsilon_mean

print(
  ggplot(D_o, aes(x = x, y = y, color = epsilon_mean)) +
    geom_point(size = 3) +
    scale_color_gradientn(colours = c("blue", "green", "yellow", "red")) +
    labs(title = "Spatial distribution of the mean of epsilon(s)",
         x = "x coordinate", y = "y coordinate", color = "mean")
)


# =============================================================================
# (5) Posterior prediction f(y(s) | y_o) at the unobserved locations
# =============================================================================

n_Du_star <- 20
indices_Du_star <- sample(1:nrow(D_u), n_Du_star)
D_u_star <- D_u[indices_Du_star, ]
true_values <- y_u[indices_Du_star]

distance_matrix_u <- as.matrix(dist(D_u_star))

y_s_samples <- matrix(NA, nrow = n_iter - burn_in, ncol = n_Du_star)

for (j in (burn_in + 1):n_iter)
{
  m_s_star <- beta0_samples[j] + D_u_star$x * beta1_samples[j]

  covariance_matrix_star <- sigma2_samples[j] * exp(-phi_samples[j] * distance_matrix_u)
  W_s_star <- mvrnorm(1, mu = m_s_star, Sigma = covariance_matrix_star)
  epsilon_star <- rnorm(n_Du_star, 0, sqrt(tau2_samples[j]))

  y_s_samples[j - burn_in, ] <- W_s_star + epsilon_star
}

# 95% credible intervals
conf_intervals <- matrix(NA, n_Du_star, 2)

for (i in 1:n_Du_star)
{
  conf_intervals[i, 1] <- quantile(y_s_samples[, i], 0.025)
  conf_intervals[i, 2] <- quantile(y_s_samples[, i], 0.975)
}

is_inside <- sapply(1:n_Du_star, function(i) {
  true_values[i] >= conf_intervals[i, 1] && true_values[i] <= conf_intervals[i, 2]
})

cat("True values falling inside the interval:", sum(is_inside), "out of", n_Du_star, "\n")

df_pred <- data.frame(
  location = 1:n_Du_star,
  true_value = true_values,
  lower = conf_intervals[, 1],
  upper = conf_intervals[, 2],
  post_mean = apply(y_s_samples, 2, mean),
  inside = is_inside
)

print(
  ggplot(df_pred, aes(x = location, y = post_mean)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, color = "grey40") +
    geom_point(aes(y = true_value, color = inside), size = 2.5) +
    scale_color_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
    labs(title = "95% credible intervals and true values",
         x = "Location in D_u*", y = "y(s)", color = "Inside")
)
