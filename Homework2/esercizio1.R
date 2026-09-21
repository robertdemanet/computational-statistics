# =============================================================================
# Homework 2 - Statistica Computazionale
# Esercizio 1: modello spaziale gerarchico
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
# (1) Simulazione di 100 osservazioni dal modello
# =============================================================================

n <- 100
s1 <- runif(n, 0, 10)  # Coordinate x
s2 <- runif(n, 0, 10)  # Coordinate y

points <- cbind(s1, s2)  # Matrice con le coordinate dei punti (x, y)

# Parametri del modello
tau2 <- 0.5
sigma2 <- 0.5
phi <- 3 / 10
beta0 <- 2
beta1 <- 0.1

# Media e matrice di covarianza di W(s)
m_s <- beta0 + s1 * beta1
distance_matrix <- as.matrix(dist(points))  # Matrice delle distanze
covariance_matrix <- sigma2 * exp(-phi * distance_matrix)

# Simulazione del processo gaussiano W(s)
W_s <- mvrnorm(1, mu = m_s, Sigma = covariance_matrix)

# Y(s) = W(s) + epsilon,  epsilon ~ N(0, tau2)
epsilon <- rnorm(n, 0, sqrt(tau2))
Y_s <- W_s + epsilon

data_sim <- data.frame(x = s1, y = s2, W_s = W_s, Y_s = Y_s)
print(head(data_sim, 3))

# Rappresentazione spaziale di W(s) e di Y(s)
print(
  ggplot(data_sim, aes(x = x, y = y, color = W_s)) +
    geom_point(size = 3) +
    scale_color_gradientn(colours = c("blue", "green", "yellow", "red")) +
    labs(title = "Simulazione del processo gaussiano W(s)",
         x = "Coordinate x", y = "Coordinate y", color = "W(s)")
)

print(
  ggplot(data_sim, aes(x = x, y = y, color = Y_s)) +
    geom_point(size = 3) +
    scale_color_gradientn(colours = c("blue", "green", "yellow", "red")) +
    labs(title = "Simulazione del processo gaussiano Y(s)",
         x = "Coordinate x", y = "Coordinate y", color = "Y(s)")
)


# =============================================================================
# (2) Scatterplot delle coordinate con 4 gruppi definiti dai quartili di y
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
    labs(title = "Scatterplot delle Coordinate con 4 Gruppi",
         x = "Coordinate x", y = "Coordinate y", color = "Gruppo di Quantile")
)


# =============================================================================
# (3) Algoritmo MCMC per f(w, beta0, beta1, tau2, sigma2, phi | y_o)
#
#     Prior:  beta ~ N_2(0, 100 * I_2),  sigma2 ~ IG(1,1),  tau2 ~ IG(1,1),
#             phi ~ Gamma(1,1),  w ~ N(X beta, C(phi, sigma2))
# =============================================================================

# --- Selezione casuale delle osservazioni ------------------------------------

n_new <- runif(1, 10, 90)
n_new <- floor(n_new)

indices <- sample(1:100, n_new)
indices <- sort(indices)

y_o <- data_sim$Y_s[indices]
D_o <- data_sim[indices, c('x', 'y')]

y_u <- data_sim$Y_s[-indices]
D_u <- data_sim[-indices, c('x', 'y')]

cat("n_new =", n_new, "\n")

# --- Inizializzazione dell'algoritmo -----------------------------------------

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

# Coordinate x e matrice delle distanze delle sole osservazioni selezionate
D_o_x <- data_sim[indices, c('x')]
distance_matrix_o <- as.matrix(dist(D_o))

# Matrice delle covariate
ones <- rep(1, n_new)
X <- cbind(ones, D_o_x)

# Matrice di precisione della prior su beta
diag_matrix <- diag(1 / 100, 2)

# --- Full conditional (a meno di costanti) di phi ----------------------------

log_fc_phi <- function(phi, w, beta, sigma2, X, C)
{
  # Calcolo del determinante
  det_C <- determinant(C, logarithm = TRUE)$modulus  # logaritmo del determinante

  det_C <- abs(det_C)  # Usato abs nel caso il determinante sia negativo

  # Calcolo della log-likelihood
  res <- w - X %*% beta
  log_likelihood <- -0.5 * (det_C + t(res) %*% solve(C) %*% res)

  # Calcolo del log-prior di phi usando la distribuzione Gamma
  log_prior <- dgamma(phi, 1, 1, log = TRUE)

  return(log_likelihood + log_prior)
}

# --- Parametri del passo Metropolis adattivo per phi -------------------------

sd_proposal <- 0.1
alpha_phi <- 0
nbatch <- 50
A <- 100
B <- 1000
alpha_target <- 0.234

# Limiti per phi ricavati dalle distanze spaziali
phi_lower <- 3 / max(dist(points))
phi_upper <- 3 / min(dist(points))

# --- Ciclo MCMC ---------------------------------------------------------------

for (i in 2:n_iter)
{
  # Matrice di covarianza corrente e sua inversa
  C_o <- sigma2_samples[i - 1] * exp(-phi_samples[i - 1] * distance_matrix_o)
  inv_C_o <- solve(C_o)

  ## Passo (1): aggiornamento beta
  # Calcolo Vp
  Vp <- solve(t(X) %*% inv_C_o %*% X + diag_matrix)

  # Calcolo Mp
  Mp <- Vp %*% (t(X) %*% inv_C_o %*% W_samples[i - 1, ])

  # Campiono beta
  beta <- mvrnorm(1, Mp, Vp)
  beta0_samples[i] <- beta[1]
  beta1_samples[i] <- beta[2]

  ## Passo (2): aggiornamento sigma2
  a_sigma2 <- n_new / 2 + 1
  prod <- t(W_samples[i - 1, ] - X %*% beta) %*% inv_C_o %*%
    (W_samples[i - 1, ] - X %*% beta)
  b_sigma2 <- prod / 2 + 1

  # Campiono sigma2
  sigma2_samples[i] <- 1 / rgamma(1, a_sigma2, b_sigma2)

  ## Passo (3): aggiornamento tau2
  a_tau2 <- n_new / 2 + 1
  b_tau2 <- 0.5 * t(y_o - W_samples[i - 1, ]) %*% (y_o - W_samples[i - 1, ]) + 1

  # Campiono tau2
  tau2_samples[i] <- 1 / rgamma(1, a_tau2, b_tau2)

  ## Passo (4): aggiornamento phi (Metropolis)
  phi <- phi_samples[i - 1]
  phi_prop <- rnorm(1, phi_samples[i - 1], sd_proposal)

  if (phi_prop >= phi_lower && phi_prop <= phi_upper)
  {
    # Le covarianze vanno ricalcolate con il valore di phi a cui si riferiscono
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

  # Adattamento della sd della proposal verso il tasso di accettazione target
  if (i %% nbatch == 0)
  {
    alpha_phi <- alpha_phi / nbatch
    sd_proposal <- exp(log(sd_proposal) + A / (B + i) * (alpha_phi - alpha_target))
    alpha_phi <- 0
  }

  phi_samples[i] <- phi

  ## Passo (5): aggiornamento W
  Q <- inv_C_o + diag(1 / tau2_samples[i], n_new)
  b <- inv_C_o %*% X %*% beta + y_o / tau2_samples[i]
  muW_post <- solve(Q) %*% b
  sigmaW_post <- solve(Q)

  # Campiono W
  W_samples[i, ] <- mvtnorm::rmvnorm(1, as.vector(muW_post), sigmaW_post)
}

# --- Tracce e numero di campioni indipendenti --------------------------------

keep <- (burn_in + 1):n_iter

plot(beta0_samples[keep], type = "l", main = "beta0", xlab = "Iterazioni", ylab = "beta0")
plot(beta1_samples[keep], type = "l", main = "beta1", xlab = "Iterazioni", ylab = "beta1")
plot(tau2_samples[keep], type = "l", main = "tau2", xlab = "Iterazioni", ylab = "tau2")
plot(sigma2_samples[keep], type = "l", main = "sigma2", xlab = "Iterazioni", ylab = "sigma2")
plot(phi_samples[keep], type = "l", main = "phi", xlab = "Iterazioni", ylab = "phi")

eff_samples <- data.frame(
  parametro = c("beta0", "beta1", "tau2", "sigma2", "phi"),
  campioni_indipendenti = c(effectiveSize(beta0_samples[keep]),
                            effectiveSize(beta1_samples[keep]),
                            effectiveSize(tau2_samples[keep]),
                            effectiveSize(sigma2_samples[keep]),
                            effectiveSize(phi_samples[keep]))
)
print(eff_samples)


# =============================================================================
# (4) A-posteriori di epsilon(s) = Y(s) - (beta0 + beta1 * s1),  s in D_o
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
    labs(title = "Distribuzione spaziale della media di epsilon(s)",
         x = "Coordinate x", y = "Coordinate y", color = "media")
)


# =============================================================================
# (5) Predizione a-posteriori f(y(s) | y_o) nei punti non osservati
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

# Intervalli di credibilita' al 95%
conf_intervals <- matrix(NA, n_Du_star, 2)

for (i in 1:n_Du_star)
{
  conf_intervals[i, 1] <- quantile(y_s_samples[, i], 0.025)
  conf_intervals[i, 2] <- quantile(y_s_samples[, i], 0.975)
}

is_inside <- sapply(1:n_Du_star, function(i) {
  true_values[i] >= conf_intervals[i, 1] && true_values[i] <= conf_intervals[i, 2]
})

cat("Valori veri contenuti nell'intervallo:", sum(is_inside), "su", n_Du_star, "\n")

df_pred <- data.frame(
  punto = 1:n_Du_star,
  vero = true_values,
  lower = conf_intervals[, 1],
  upper = conf_intervals[, 2],
  media = apply(y_s_samples, 2, mean),
  dentro = is_inside
)

print(
  ggplot(df_pred, aes(x = punto, y = media)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, color = "grey40") +
    geom_point(aes(y = vero, color = dentro), size = 2.5) +
    scale_color_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
    labs(title = "Intervalli di credibilita' al 95% e valori veri",
         x = "Punto di D_u*", y = "y(s)", color = "Nell'intervallo")
)
