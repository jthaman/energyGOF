# EnergyGOF

Conduct one- and two-sample goodness-of-fit tests for univariate data. In the one-sample case, normal, uniform, exponential, Bernoulli, binomial, geometric, beta, Poisson, lognormal, Laplace, asymmetric Laplace, inverse Gaussian, half-normal, chi-squared, gamma, F, Weibull, Cauchy, and Pareto distributions are supported. egof.test() can also test goodness-of-fit to any distribution with a continuous distribution function. A subset of the available distributions can be tested for the composite goodness-of-fit hypothesis, that is, one can test for distribution fit with unknown parameters. P-values are calculated via parametric bootstrap.

# Examples

```r
x <- rnorm(10)
y <- rt(10, 4)

## Composite energy goodness-of-fit test (test for Normality with unknown
## parameters)

energyGOF.test(x, "normal", nsim = 10)

## Simple energy goodness-of-fit test (test for Normality with known
## parameters). egof.test is an alias for energyGOF.test.

egof.test(x, "normal", nsim = 10, mean = 0, sd = 1)

## Alternatively, use the energyGOFdist generic directly so that you do not need
## to pass parameter names into `...`

energyGOFdist(x, normal_dist(0, 1), nsim = 10)

## Conduct a two-sample test

egof.test(x, y, 0)

## Conduct a test against any continuous distribution function

egof.test(x, pcauchy, 0)

## Simple energy goodness-of-fit test for Weibull distribution

y <- rweibull(10, 1, 1)
energyGOF.test(y, "weibull", shape = 1, scale = 3, nsim = 10)

## Alternatively, use the energyGOFdist generic directly, which is slightly less
## verbose. egofd is an alias for energyGOFdist.

egofd(y, weibull_dist(1, 3), nsim = 10)

## Conduct a generalized GOF test. `pow` is the exponent *s* in the generalize ## energy statistic. Pow is only necessary when testing Cauchy, and
## Pareto distributions. If you don't set a pow, there is a default for each
## of the distributions, but the default isn't necessarily better than any
## other number.

egofd(rcauchy(100),
      cauchy_dist(location = 0, scale = 1, pow = 0.5),
      nsim = 10)

## energyGOF does not support tests with a mix of known and unknown
## parameters, so this will result in an error.

energyGOF.test(x, "normal", mean = 0, nsim = 10) # sd is missing
```
