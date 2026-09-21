# Homework 1 — Metodi Monte Carlo

## `esercizio1.R`

Modello mistura con variabile latente:

```
Z ~ Bern(p),   Y|Z = 0 ~ Poisson(lambda),   Y|Z = 1 ~ Gamma(a, b)
```

con `a ~ U(2.5, 5.5)`, `lambda ~ U(2, 4)`, `b = 1`, `p = 0.5`, `n = 1000`.

Punti svolti:

1. simulazione dei campioni dalla marginale di `Y` e stima Monte Carlo della CDF
   (`F(t) = E[I{Y <= t}]`);
2. stima di `P(Y in [3.5, 4.5])` e `P(Y in [3.5, 4.5] | Z = 0)`, con la varianza
   dei due stimatori;
3. `P(Z = 0 | Y in [1.5, 3.5])` e `P(Z = 1 | Y in [1.5, 3.5])` via teorema di Bayes;
4. quantili empirici di `Y` e `Z` ai livelli 0.1, 0.2, 0.5, 0.75 tramite
   l'inversa generalizzata `F^-1(u) = inf{x : u <= F(x)}`;
5. stima della parte discreta (`Z = 0`) e della parte continua (`Z = 1`) di `Y`
   sui punti `seq(0, 10, 0.25)`, con confronto rispetto alle densita' teoriche.

## `esercizio2.R`

Distribuzione **logistic-normal** su `[0, 1]`:

```
f(x | mu, sigma2) = 1/(x(1-x)) * 1/sqrt(2 pi sigma2) * exp(-(logit(x) - mu)^2 / (2 sigma2))
```

con `mu ~ U(-1.5, 1.5)` e `sigma2 ~ U(0.5, 1.5)`.

Punti svolti:

1. simulazione dei dati con il metodo **Accept-Reject** (proposta uniforme sul
   rettangolo `[0,1] x [0,M]`, con `M` massimo della densita');
2. derivazione della a-posteriori di `mu` con prior `mu ~ N(0, 100)` e `sigma2` noto.
   Il kernel risulta gaussiano:

   ```
   pi(mu | x) ~ N( sum(logit(x_i)) / (n + sigma2/100),  sigma2 / (n + sigma2/100) )
   ```

   quindi il campionamento e' diretto;
3. confronto tra prior e posterior al crescere delle osservazioni (10 vs 100);
4. stima di `E(X)` con **Importance Sampling** usando come densita' strumentale
   una `Beta(1.2, 1.2)`.

## `report/`

Relazione in LaTeX (`main.tex`) con le figure prodotte dagli script, in
`report/Immagini`. Riporta la derivazione della a-posteriori di `mu` e i
risultati numerici commentati.
