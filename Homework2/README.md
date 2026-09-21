# Homework 2 — Metodi MCMC

## `esercizio1.R`

Modello spaziale gerarchico su `s = (s1, s2) in R^2`:

```
Y(s) | W(s) ~ GP(W(s), tau2)
W(s)        ~ GP(m(s), C(||s - s'||; phi, sigma2))
m(s) = beta0 + beta1 * s1
C(||s - s'||; phi, sigma2) = sigma2 * exp(-phi * ||s - s'||)
```

con `tau2 = 0.5`, `sigma2 = 0.5`, `phi = 3/10`, `beta0 = 2`, `beta1 = 0.1` e
`n = 100` punti campionati su `U(0,10) x U(0,10)`.

Punti svolti:

1. simulazione del processo e rappresentazione spaziale di `W(s)` e `Y(s)`;
2. scatterplot delle coordinate con quattro gruppi definiti dai quartili di `y`;
3. algoritmo MCMC per `f(w, beta0, beta1, tau2, sigma2, phi | y_o)` su un
   sottoinsieme casuale di `n_new in [10, 90]` osservazioni. Lo schema e' un
   Gibbs sampler con un passo Metropolis:

   | Passo | Parametro | Full conditional |
   |---|---|---|
   | 1 | `beta` | `N_2(Mp, Vp)` coniugata |
   | 2 | `sigma2` | `IG(n_new/2 + a, (w - X beta)' C^-1 (w - X beta)/2 + b)` |
   | 3 | `tau2` | `IG(n_new/2 + a, (y_o - w)'(y_o - w)/2 + b)` |
   | 4 | `phi` | Metropolis random-walk con sd adattiva (target 0.234) |
   | 5 | `w` | `N(Q^-1 b, Q^-1)`, con `Q = I/tau2 + C^-1` |

   Il numero di campioni indipendenti e' calcolato con `coda::effectiveSize`;
4. a-posteriori del residuo `epsilon(s) = Y(s) - (beta0 + beta1 s1)` sui punti
   osservati, con la relativa distribuzione spaziale della media;
5. predizione a-posteriori `f(y(s) | y_o)` su 20 punti non osservati, con
   intervalli di credibilita' al 95% e verifica della copertura dei valori veri.

Il passo Metropolis su `phi` ricalcola la matrice di covarianza `C` sia nel valore
corrente sia in quello proposto: e' la condizione perche' il rapporto di
accettazione dipenda effettivamente dalla verosimiglianza e non solo dal prior.

## `esercizio2.R`

Modello mistura di Poisson con variabile latente:

```
Y_i | z_i ~ Poisson(lambda_{z_i}),   P(z_i = k) = pi_k
i = 1,...,200,   z_i in {1,2,3},   pi_k = 1/3
lambda_1 = 1,  lambda_2 = 10,  lambda_3 = 25
```

Punti svolti:

1. simulazione dal modello e rappresentazione grafica dei dati (osservazioni per
   gruppo, densita' empirica, densita' sovrapposte, boxplot e violin plot);
2. Gibbs sampler per `f(lambda, pi, z | y)` con prior
   `lambda_k ~ Gamma(1,1)`, `pi ~ Dir(1,1,1)`, `z_i | pi ~ Discrete(pi)`:

   | Passo | Parametro | Full conditional |
   |---|---|---|
   | 1 | `lambda_k` | `Gamma(a + sum_{i: z_i = k} y_i, b + n_k)` |
   | 2 | `pi` | `Dir(alpha_k + n_k)` |
   | 3 | `z_i` | categorica con `P(z_i = k) propto pi_k Pois(y_i | lambda_k)` |

   L'aggiornamento di `z_i` e' normalizzato in scala logaritmica (log-sum-exp)
   per evitare underflow quando le `lambda` sono ben separate.
