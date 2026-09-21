# =============================================================================
# Homework 1 - Computational Statistics
# Exercise 1: Poisson/Gamma mixture, Monte Carlo estimation
#
#   Z ~ Bern(p)
#   Y|Z = 0 ~ Poisson(lambda)
#   Y|Z = 1 ~ Gamma(a, b)
# =============================================================================

set.seed(343240)

# ---- Parameters -------------------------------------------------------------

a <- runif(1, 2.5, 5.5)   # 4.706901
lambda <- runif(1, 2, 4)  # 3.049805
b <- 1
p <- 0.5
n <- 1000

cat("a =", a, "\n")
cat("lambda =", lambda, "\n")


# =============================================================================
# (a) Samples from the marginal of Y and Monte Carlo estimate of the CDF
# =============================================================================

z <- rbinom(n, 1, p)

y_sample <- rep(NA, n)     # samples from the marginal of Y
y_sample_z0 <- rep(NA, n)  # samples from Y|Z = 0 (-1 where Z = 1)
y_sample_z1 <- rep(NA, n)  # samples from Y|Z = 1 (-1 where Z = 0)

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

# Actual draws of the two components (the value -1 is only a placeholder)
y_z0 <- y_sample_z0[y_sample_z0 != -1]
y_z1 <- y_sample_z1[y_sample_z1 != -1]
n1 <- length(y_z0)  # number of draws from Y|Z = 0
n2 <- length(y_z1)  # number of draws from Y|Z = 1

# Monte Carlo estimate of the CDF: F(t) = E[I{Y <= t}] ~ mean of the indicators
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
# (b) Monte Carlo estimate of P(Y in [3.5, 4.5]) and P(Y in [3.5, 4.5] | Z = 0),
#     together with the variance of the estimators
# =============================================================================

ind_y <- (y_sample >= 3.5 & y_sample <= 4.5)
ind_y_z0 <- (y_z0 >= 3.5 & y_z0 <= 4.5)

p_y <- sum(ind_y) / n
p_y_z0 <- sum(ind_y_z0) / n1

cat("P(Y in [3.5,4.5]) =", p_y, "\n")
cat("P(Y in [3.5,4.5] | Z = 0) =", p_y_z0, "\n")

# Variance of the Monte Carlo estimator: v_n = (1/n^2) * sum of squared deviations
v_n_y <- sum((ind_y - p_y)^2) / n^2
v_n_yz0 <- sum(((y_sample_z0 >= 3.5 & y_sample_z0 <= 4.5) - p_y_z0)^2) / n1^2

cat("variance of the estimator on Y =", v_n_y, "\n")
cat("variance of the estimator on Y|Z=0 =", v_n_yz0, "\n")


# =============================================================================
# (c) Monte Carlo estimate of P(Z = 0 | Y in [1.5, 3.5]) and
#     P(Z = 1 | Y in [1.5, 3.5]) through Bayes' theorem
# =============================================================================

p_A <- sum(y_sample >= 1.5 & y_sample <= 3.5) / n       # P(Y in A)
p_A_z0 <- sum(y_z0 >= 1.5 & y_z0 <= 3.5) / n1           # P(Y in A | Z = 0)
p_A_z1 <- sum(y_z1 >= 1.5 & y_z1 <= 3.5) / n2           # P(Y in A | Z = 1)

p_z0 <- n1 / n  # Monte Carlo estimate of P(Z = 0)
p_z1 <- n2 / n  # Monte Carlo estimate of P(Z = 1)

post_z0 <- p_A_z0 * p_z0 / p_A
post_z1 <- p_A_z1 * p_z1 / p_A

cat("P(Z = 0 | Y in [1.5,3.5]) =", post_z0, "\n")
cat("P(Z = 1 | Y in [1.5,3.5]) =", post_z1, "\n")


# =============================================================================
# (d) Empirical quantiles through the generalised inverse
#     F^-1(u) = inf{ x : u <= F(x) }
# =============================================================================

u_levels <- c(0.1, 0.2, 0.5, 0.75)

# --- Quantiles of Y: the CDF estimated in point (a) is used ---
quantiles_y <- rep(NA, length(u_levels))

for (k in 1:length(u_levels))
{
  idx <- which(u_levels[k] <= cdf)[1]
  quantiles_y[k] <- y_seq[idx]
}

cat("Empirical quantiles of Y:\n")
print(data.frame(level = u_levels, value = quantiles_y))

# --- Quantiles of Z: the CDF is built explicitly ---
z_seq <- seq(min(z), max(z), length.out = n)
cdf_z <- rep(NA, n)

for (i in 1:n)
{
  cdf_z[i] <- sum(z <= z_seq[i]) / n
}

plot(z_seq, cdf_z, type = "l", col = "blue", lwd = 2,
     main = "Cumulative Distribution of Z",
     xlab = "z_seq", ylab = "cdf")

quantiles_z <- rep(NA, length(u_levels))

for (k in 1:length(u_levels))
{
  idx <- which(u_levels[k] <= cdf_z)[1]
  quantiles_z[k] <- z_seq[idx]
}

cat("Empirical quantiles of Z:\n")
print(data.frame(level = u_levels, value = quantiles_z))


# =============================================================================
# (e) Monte Carlo estimate of the discrete part (Z = 0) and of the
#     continuous part (Z = 1) of Y
# =============================================================================

x <- seq(0, 10, by = 0.25)  # 41 points on which the two parts are estimated
tol <- 0.5                  # tolerance required by the continuous part

est_y_z0 <- rep(NA, length(x))
est_y_z1 <- rep(NA, length(x))

for (i in 1:length(x))
{
  # discrete part: relative frequency of the exact value
  est_y_z0[i] <- sum(y_z0 == x[i]) / n1

  # continuous part: relative frequency within (x - tol, x + tol)
  est_y_z1[i] <- sum(y_z1 > x[i] - tol & y_z1 < x[i] + tol) / n2
}

# --- Discrete part against the theoretical Poisson ---
plot(x, est_y_z0, type = "l", col = "red", lwd = 2,
     ylim = c(0, max(est_y_z0, dpois(round(x), lambda))),
     main = "Empirical vs Theoretical Distribution of Y | Z = 0",
     xlab = "x values", ylab = "Density")
lines(x, dpois(round(x), lambda), col = "green", lwd = 2)
legend("topright", legend = c("Empirical", "Theoretical (Poisson)"),
       col = c("red", "green"), lwd = 2)

# --- Continuous part against the theoretical Gamma ---
plot(x, est_y_z1, type = "l", col = "blue", lwd = 2,
     ylim = c(0, max(est_y_z1, dgamma(x, shape = a, rate = b))),
     main = "Empirical vs Theoretical Distribution of Y | Z = 1",
     xlab = "x values", ylab = "Density")
lines(x, dgamma(x, shape = a, rate = b), col = "green", lwd = 2)
legend("topright", legend = c("Empirical", "Theoretical (Gamma)"),
       col = c("blue", "green"), lwd = 2)
