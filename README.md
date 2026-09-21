# Computational Statistics — Homework

R code for the two homework assignments of the Computational Statistics course
(MSc).

| Folder | Content |
|---|---|
| [`Homework1/`](Homework1) | Monte Carlo methods: Poisson/Gamma mixture, Accept-Reject, Importance Sampling |
| [`Homework2/`](Homework2) | MCMC methods: hierarchical spatial model, Poisson mixture with latent variable |

## Layout

```
Homework1/
  exercise1.R         Poisson/Gamma mixture: CDF, probabilities, quantiles, discrete and continuous parts
  exercise2.R         logistic-normal: Accept-Reject, posterior of mu, Importance Sampling
  report/main.tex     report, with the figures in report/images
Homework2/
  exercise1.R         spatial Gaussian process: Gibbs + adaptive Metropolis, spatial prediction
  exercise2.R         three-component Poisson mixture: Gibbs sampler with latent variable
  report/main.tex     report, with the figures in report/images
```

## Running the code

Each script is self-contained and is run from its own folder:

```bash
Rscript Homework1/exercise1.R
```

The seed is fixed at the top of every script (`set.seed(343240)`), so the results
are reproducible.

The reports are compiled from the corresponding `report` folder:

```bash
cd Homework1/report && pdflatex main.tex && pdflatex main.tex
```

The second pass is needed to resolve the table of contents.

## Required packages

- Homework 1: none (base R only)
- Homework 2: `MASS`, `mvtnorm`, `coda`, `gtools`, `ggplot2`, `tidyr`

```r
install.packages(c("MASS", "mvtnorm", "coda", "gtools", "ggplot2", "tidyr"))
```

## Notes

`Homework2/exercise1.R` runs 60000 MCMC iterations (burn-in 20000) and takes a
few minutes.
