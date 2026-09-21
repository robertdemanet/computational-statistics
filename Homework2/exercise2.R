# =============================================================================
# Homework 2 - Computational Statistics
# Exercise 2: Poisson mixture model with latent variable
#
#   Y_i | z_i ~ Poisson(lambda_{z_i}),   P(z_i = k) = pi_k
#   i = 1,...,200,  z_i in {1,2,3},  pi_k = 1/3,
#   lambda_1 = 1, lambda_2 = 10, lambda_3 = 25
# =============================================================================

library(gtools)   # rdirichlet
library(ggplot2)
library(tidyr)

set.seed(343240)


# =============================================================================
# (1) Simulation from the model and graphical representation of the data
# =============================================================================

n <- 200
Y <- rep(NA, n)
lambda_vec <- c(1, 10, 25)
z <- rep(NA, n)

for (i in 1:n)
{
  sim <- sample(1:3, 1)
  z[i] <- sim
  Y[i] <- rpois(1, lambda_vec[sim])
}

data_sim <- data.frame(index = 1:n, Y = Y, group = factor(z))

# Observations split by group
print(
  ggplot(data_sim, aes(x = index, y = Y, color = group)) +
    geom_point(size = 2) +
    scale_color_manual(values = c("blue", "green", "red")) +
    labs(title = "Observations simulated from Y", x = "Index", y = "Y",
         color = "Group")
)

# Empirical density of Y
print(
  ggplot(data_sim, aes(x = Y)) +
    geom_density(fill = "lightblue", alpha = 0.6) +
    labs(title = "Empirical density of Y", x = "Y", y = "Density")
)

# Overlaid group densities
print(
  ggplot(data_sim, aes(x = Y, fill = group)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("blue", "green", "red")) +
    labs(title = "Distribution of the observations by group",
         x = "Y", y = "Density", fill = "Group")
)

# Boxplots and densities by group
print(
  ggplot(data_sim, aes(x = group, y = Y, fill = group)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.15, alpha = 0.8) +
    scale_fill_manual(values = c("blue", "green", "red")) +
    labs(title = "Boxplots and distribution of the observations by group",
         x = "Group", y = "Y", fill = "Group")
)


# =============================================================================
# (2) Gibbs sampler for f(lambda_1, lambda_2, lambda_3, pi, z_1,...,z_n | y)
#
#     Priors: lambda_k ~ Gamma(a_prior, b_prior),  pi ~ Dir(alpha_prior),
#             z_i | pi ~ Discrete(pi)
# =============================================================================

n_iter <- 10000
k <- 3
a_gamma_prior <- 1
b_gamma_prior <- 1
alpha_prior <- c(1, 1, 1)

lambda_samples <- matrix(NA, nrow = n_iter, ncol = k)
pi_samples <- matrix(NA, nrow = n_iter, ncol = k)
zi_samples <- matrix(NA, nrow = n_iter, ncol = n)

# Starting values
lambda <- rep(mean(Y), k)
pi <- rep(1 / k, k)
z <- sample(1:k, n, replace = TRUE)

for (i in 1:n_iter)
{
  ## Step (1): update of lambda_k
  ##   lambda_k | y, z ~ Gamma(a_prior + sum_{i: z_i = k} y_i, b_prior + n_k)
  for (j in 1:k)
  {
    Y_j <- Y[z == j]
    lambda[j] <- rgamma(1, a_gamma_prior + sum(Y_j), b_gamma_prior + length(Y_j))
  }

  lambda_samples[i, ] <- lambda

  ## Step (2): update of pi
  ##   pi | z, y ~ Dir(alpha_k + n_k)
  n_k <- table(factor(z, levels = 1:k))
  pi <- rdirichlet(1, alpha_prior + n_k)
  pi_samples[i, ] <- pi

  ## Step (3): update of z_i
  ##   P(z_i = k | y_i, lambda, pi) propto pi_k * Pois(y_i | lambda_k)
  for (j in 1:n)
  {
    log_prob_non_norm <- c(0, 0, 0)
    log_prob_non_norm[1] <- log(pi[1]) + dpois(Y[j], lambda[1], log = TRUE)
    log_prob_non_norm[2] <- log(pi[2]) + dpois(Y[j], lambda[2], log = TRUE)
    log_prob_non_norm[3] <- log(pi[3]) + dpois(Y[j], lambda[3], log = TRUE)

    # Normalisation on the log scale (log-sum-exp)
    c <- max(log_prob_non_norm)
    prob <- exp(log_prob_non_norm - (c + log(sum(exp(log_prob_non_norm - c)))))

    z[j] <- sample(1:k, 1, prob = prob)
  }

  zi_samples[i, ] <- z
}


# =============================================================================
# Results
# =============================================================================

estimates <- data.frame(
  parameter = c("lambda_1", "lambda_2", "lambda_3", "pi_1", "pi_2", "pi_3"),
  posterior_mean = c(colMeans(lambda_samples), colMeans(pi_samples)),
  true_value = c(lambda_vec, rep(1 / 3, k))
)
print(estimates)

# MCMC traces of the parameters
traces <- data.frame(
  iteration = 1:n_iter,
  lambda_1 = lambda_samples[, 1],
  lambda_2 = lambda_samples[, 2],
  lambda_3 = lambda_samples[, 3],
  pi_1 = pi_samples[, 1],
  pi_2 = pi_samples[, 2],
  pi_3 = pi_samples[, 3],
  zi = zi_samples[, 1]
)

traces_long <- pivot_longer(traces, cols = -iteration,
                            names_to = "parameter", values_to = "value")
traces_long$parameter <- factor(traces_long$parameter,
                                levels = c("lambda_1", "lambda_2", "lambda_3",
                                           "pi_1", "pi_2", "pi_3", "zi"))

print(
  ggplot(traces_long, aes(x = iteration, y = value)) +
    geom_line(linewidth = 0.2) +
    facet_wrap(~ parameter, scales = "free_y") +
    labs(title = "MCMC traces of the parameters", x = "Iterations", y = "Value")
)
