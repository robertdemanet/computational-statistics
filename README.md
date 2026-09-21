# Statistica Computazionale — Homework

Codice R dei due homework del corso di Statistica Computazionale (Laurea Magistrale).

| Cartella | Contenuto |
|---|---|
| [`Homework1/`](Homework1) | Metodi Monte Carlo: mistura Poisson/Gamma, Accept-Reject, Importance Sampling |
| [`Homework2/`](Homework2) | Metodi MCMC: modello spaziale gerarchico, mistura di Poisson con variabile latente |

## Struttura

```
Homework1/
  esercizio1.R   mistura Poisson/Gamma: CDF, probabilita', quantili, parti discreta e continua
  esercizio2.R   logistic-normal: Accept-Reject, a-posteriori di mu, Importance Sampling
Homework2/
  esercizio1.R   processo gaussiano spaziale: Gibbs + Metropolis adattivo, predizione spaziale
  esercizio2.R   mistura di Poisson a 3 componenti: Gibbs sampler con variabile latente
```

## Esecuzione

Ogni script e' indipendente e si esegue dalla propria cartella:

```bash
Rscript Homework1/esercizio1.R
```

Il seed e' fissato all'inizio di ogni script (`set.seed(343240)`), quindi i risultati
sono riproducibili.

## Pacchetti richiesti

- Homework 1: nessuno (solo funzioni di base)
- Homework 2: `MASS`, `mvtnorm`, `coda`, `gtools`, `ggplot2`, `tidyr`

```r
install.packages(c("MASS", "mvtnorm", "coda", "gtools", "ggplot2", "tidyr"))
```

## Note

Lo script `Homework2/esercizio1.R` esegue 60000 iterazioni MCMC (burn-in 20000) e
richiede qualche minuto.
