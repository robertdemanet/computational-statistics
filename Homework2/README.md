# Homework 2 — MCMC methods

## `exercise1.R`

Hierarchical spatial model on `s = (s1, s2) in R^2`:

```
Y(s) | W(s) ~ GP(W(s), tau2)
W(s)        ~ GP(m(s), C(||s - s'||; phi, sigma2))
m(s) = beta0 + beta1 * s1
C(||s - s'||; phi, sigma2) = sigma2 * exp(-phi * ||s - s'||)
```

with `tau2 = 0.5`, `sigma2 = 0.5`, `phi = 3/10`, `beta0 = 2`, `beta1 = 0.1` and
`n = 100` locations sampled on `U(0,10) x U(0,10)`.

Points covered:

1. simulation of the process and spatial representation of `W(s)` and `Y(s)`;
2. scatterplot of the coordinates with four groups defined by the quartiles of
   `y`;
3. MCMC algorithm for `f(w, beta0, beta1, tau2, sigma2, phi | y_o)` on a random
   subset of `n_new in [10, 90]` observations. The scheme is a Gibbs sampler with
   one Metropolis step:

   | Step | Parameter | Full conditional |
   |---|---|---|
   | 1 | `beta` | conjugate `N_2(Mp, Vp)` |
   | 2 | `sigma2` | `IG(n_new/2 + a, (w - X beta)' C^-1 (w - X beta)/2 + b)` |
   | 3 | `tau2` | `IG(n_new/2 + a, (y_o - w)'(y_o - w)/2 + b)` |
   | 4 | `phi` | random-walk Metropolis with adaptive sd (target 0.234) |
   | 5 | `w` | `N(Q^-1 b, Q^-1)`, with `Q = I/tau2 + C^-1` |

   The effective sample size is computed with `coda::effectiveSize`;
4. posterior of the residual `epsilon(s) = Y(s) - (beta0 + beta1 s1)` at the
   observed locations, together with the spatial distribution of its mean;
5. posterior prediction `f(y(s) | y_o)` at 20 unobserved locations, with 95%
   credible intervals and a check of the coverage of the true values.

The Metropolis step on `phi` recomputes the covariance matrix `C` both at the
current and at the proposed value: this is what makes the acceptance ratio depend
on the likelihood and not on the prior alone.

## `exercise2.R`

Poisson mixture model with a latent variable:

```
Y_i | z_i ~ Poisson(lambda_{z_i}),   P(z_i = k) = pi_k
i = 1,...,200,   z_i in {1,2,3},   pi_k = 1/3
lambda_1 = 1,  lambda_2 = 10,  lambda_3 = 25
```

Points covered:

1. simulation from the model and graphical representation of the data
   (observations by group, empirical density, overlaid densities, boxplots and
   violin plots);
2. Gibbs sampler for `f(lambda, pi, z | y)` under the priors
   `lambda_k ~ Gamma(1,1)`, `pi ~ Dir(1,1,1)`, `z_i | pi ~ Discrete(pi)`:

   | Step | Parameter | Full conditional |
   |---|---|---|
   | 1 | `lambda_k` | `Gamma(a + sum_{i: z_i = k} y_i, b + n_k)` |
   | 2 | `pi` | `Dir(alpha_k + n_k)` |
   | 3 | `z_i` | categorical, `P(z_i = k) propto pi_k Pois(y_i | lambda_k)` |

   The update of `z_i` is normalised on the log scale (log-sum-exp) to avoid
   underflow when the `lambda` are well separated.

## `report/`

LaTeX report (`main.tex`) together with the figures produced by the scripts, in
`report/images`. It contains the DAG of the spatial model, the derivation of the
full conditionals and the commented numerical results.
