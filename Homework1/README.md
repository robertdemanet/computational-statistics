# Homework 1 — Monte Carlo methods

## `exercise1.R`

Mixture model with a latent variable:

```
Z ~ Bern(p),   Y|Z = 0 ~ Poisson(lambda),   Y|Z = 1 ~ Gamma(a, b)
```

with `a ~ U(2.5, 5.5)`, `lambda ~ U(2, 4)`, `b = 1`, `p = 0.5`, `n = 1000`.

Points covered:

1. simulation of the samples from the marginal of `Y` and Monte Carlo estimate of
   the CDF (`F(t) = E[I{Y <= t}]`);
2. estimation of `P(Y in [3.5, 4.5])` and `P(Y in [3.5, 4.5] | Z = 0)`, together
   with the variance of the two estimators;
3. `P(Z = 0 | Y in [1.5, 3.5])` and `P(Z = 1 | Y in [1.5, 3.5])` through Bayes'
   theorem;
4. empirical quantiles of `Y` and `Z` at levels 0.1, 0.2, 0.5, 0.75 through the
   generalised inverse `F^-1(u) = inf{x : u <= F(x)}`;
5. estimation of the discrete part (`Z = 0`) and of the continuous part (`Z = 1`)
   of `Y` on the grid `seq(0, 10, 0.25)`, compared against the theoretical
   densities.

## `exercise2.R`

**Logistic-normal** distribution on `[0, 1]`:

```
f(x | mu, sigma2) = 1/(x(1-x)) * 1/sqrt(2 pi sigma2) * exp(-(logit(x) - mu)^2 / (2 sigma2))
```

with `mu ~ U(-1.5, 1.5)` and `sigma2 ~ U(0.5, 1.5)`.

Points covered:

1. data simulation through the **Accept-Reject** method (uniform proposal on the
   rectangle `[0,1] x [0,M]`, where `M` is the maximum of the density);
2. derivation of the posterior of `mu` under the prior `mu ~ N(0, 100)` with
   `sigma2` known. The kernel turns out to be Gaussian:

   ```
   pi(mu | x) ~ N( sum(logit(x_i)) / (n + sigma2/100),  sigma2 / (n + sigma2/100) )
   ```

   so sampling is direct;
3. comparison between prior and posterior as the number of observations grows
   (10 vs 100);
4. estimation of `E(X)` through **Importance Sampling**, using a `Beta(1.2, 1.2)`
   as the instrumental density.

## `report/`

LaTeX report (`main.tex`) together with the figures produced by the scripts, in
`report/images`. It contains the derivation of the posterior of `mu` and the
commented numerical results.
