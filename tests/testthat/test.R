set.seed(730)

##### Top Tests

test_that("formals same as switch", {
  b <- names(as.list(body(char_to_dist))[[2]][-c(1, 2)])
  b <- b[b != ""] # drop the ""
  f <- unlist(as.list(formals(ef.test)$dist)[-1])
  expect_setequal(b, f)
})


##### Check that distributions are correctly formed

test_that("Dists", {
  ## Uniform
  expect_s3_class(uniform_dist(0, 1), "EuclideanGOFDist")
  expect_s3_class(uniform_dist(-100, 100), "EuclideanGOFDist")
  expect_error(uniform_dist(100, -100))
  ## Normal
  expect_s3_class(normal_dist(0, 1), "GOFDist")
  expect_s3_class(normal_dist(), "GOFDist")
  expect_error(normal_dist(0, -1))
  expect_error(normal_dist(0, NULL))
  expect_error(normal_dist(NULL, 0))
  ## Exp
  expect_s3_class(exponential_dist(1), "GOFDist")
  expect_s3_class(exponential_dist(), "GOFDist")
  expect_error(exponential_dist(-1))
  ## Poisson
  expect_s3_class(poisson_dist(1), "GOFDist")
  expect_s3_class(poisson_dist(), "GOFDist")
  expect_error(poisson_dist(-1))
  ## Bernoulli
  expect_s3_class(bernoulli_dist(.5), "GOFDist")
  expect_s3_class(bernoulli_dist(), "GOFDist")
  expect_error(bernoulli_dist(-1))
  ## Binomial
  expect_s3_class(binomial_dist(5, .5), "GOFDist")
  expect_s3_class(binomial_dist(5), "GOFDist")
  expect_error(binomial_dist(NULL, 5))
  expect_error(binomial_dist(5, NULL))
  expect_error(binomial_dist(-1))
  ## Beta
  expect_s3_class(beta_dist(1, 1), "GOFDist")
  expect_s3_class(beta_dist(1e6, 1 / 1e5), "GOFDist")
  expect_error(beta_dist(1, NULL))
  expect_error(beta_dist(NULL, 1))
  expect_s3_class(beta_dist(), "GOFDist")
  ## Geometric
  expect_s3_class(geometric_dist(.5), "GOFDist")
  expect_s3_class(geometric_dist(), "GOFDist")
  expect_error(geometric_dist(-1))
  ## Half normal
  expect_s3_class(halfnormal_dist(.5), "GOFDist")
  expect_s3_class(halfnormal_dist(), "GOFDist")
  expect_error(halfnormal_dist(-1))
  ## Laplace
  expect_s3_class(laplace_dist(0, 1), "GOFDist")
  expect_s3_class(laplace_dist(), "GOFDist")
  expect_error(laplace_dist(0, -1))
  expect_error(laplace_dist(0, NULL))
  expect_error(laplace_dist(NULL, 0))
  ## Lognormal
  expect_s3_class(lognormal_dist(0, 1), "GOFDist")
  expect_s3_class(lognormal_dist(), "GOFDist")
  expect_error(lognormal_dist(0, -1))
  expect_error(lognormal_dist(0, NULL))
  expect_error(lognormal_dist(NULL, 0))
  ## Asym Laplace
  expect_s3_class(asymmetric_laplace_dist(0, 1, 1), "GOFDist")
  expect_s3_class(asymmetric_laplace_dist(), "GOFDist")
  expect_error(asymmetric_laplace_dist(1, 1, -1))
  expect_error(asymmetric_laplace_dist(1, -1, 1))
  expect_error(asymmetric_laplace_dist(NULL, -1, 1))
  expect_error(asymmetric_laplace_dist(1, NULL, 1))
  expect_error(asymmetric_laplace_dist(1, 1, NULL))
  ## Weibull
  expect_s3_class(weibull_dist(1, 1), "GOFDist")
  expect_s3_class(weibull_dist(), "GOFDist")
  expect_error(weibull_dist(1, NULL))
  expect_error(weibull_dist(NULL, 1))
  expect_error(weibull_dist(-1, 1))
  expect_error(weibull_dist(1, -1))
  ## Gamma
  expect_s3_class(gamma_dist(1, 1), "GOFDist")
  expect_s3_class(gamma_dist(), "GOFDist")
  expect_error(gamma_dist(1, NULL))
  expect_error(gamma_dist(NULL, 1))
  expect_error(gamma_dist(-1, 1))
  expect_error(gamma_dist(1, -1))
  ## Chi-sq
  expect_s3_class(chisq_dist(1), "GOFDist")
  expect_s3_class(chisq_dist(), "GOFDist")
  expect_error(chisq_dist(NULL))
  ## Inv Gaussian
  expect_s3_class(inverse_gaussian_dist(1, 1), "GOFDist")
  expect_s3_class(inverse_gaussian_dist(), "GOFDist")
  expect_error(inverse_gaussian_dist(1, NULL))
  expect_error(inverse_gaussian_dist(NULL, 1))
  expect_error(inverse_gaussian_dist(-1, 1))
  expect_error(inverse_gaussian_dist(1, -1))
  ## Pareto
  expect_s3_class(pareto_dist(1, 1), "GOFDist")
  expect_s3_class(pareto_dist(1e5, 1 / 1e5), "GOFDist")
  expect_error(pareto_dist(1 / 1e5, 1e5))
  expect_s3_class(pareto_dist(1e5, 1e5), "GOFDist") # Rounds to zero!
  expect_s3_class(pareto_dist(5, 5, 1), "GOFDist")
  expect_s3_class(pareto_dist(5, 5, 2), "GOFDist")
  expect_error(pareto_dist(5, 5, 3))
  expect_error(pareto_dist(5, -5, 2))
  expect_error(pareto_dist(-5, 5, 2))
  expect_error(pareto_dist(NULL, 4, 2))
  expect_error(pareto_dist(NULL, NULL, 2))
  expect_error(pareto_dist(4, NULL, 2))
  expect_error(pareto_dist(4, 4, NULL))
  ## Cauchy
  expect_s3_class(cauchy_dist(0, 1), "GOFDist")
  expect_s3_class(cauchy_dist(), "GOFDist")
  expect_error(cauchy_dist(0, -1))
  expect_error(cauchy_dist(0, NULL))
  expect_error(cauchy_dist(NULL, 0))
  ## Stable
  expect_s3_class(stable_dist(), "GOFDist")
  expect_error(stable_dist(NULL, 1, 1, 1))
  expect_error(stable_dist(0, NULL, 1, 1))
  expect_error(stable_dist(0, 1, NULL, 1))
  expect_error(stable_dist(0, 0, 1, NULL))
  expect_error(stable_dist(0, 0, 0, 0))
  expect_error(stable_dist(1, 1, 1, 1, 1))
})

##### Beta Test

## Fails
test_that("Beta Estat should be positive", {
  d <- beta_dist(20, 20)
  x <- rbeta(100, 20, 20)
  o <- ef(x, d)
  print(o)
  expect_gt(o$statistic, 0)})

test_that("Beta Estat should be positive", {
  d <- beta_dist(.5, 1.5)
  x <- rbeta(100, .5, 1.5)
  o <- ef(x, d)
  print(o)
  expect_gt(o$statistic, 0)})

test_that("Beta Estat should be sensitive to parameter change", {
  d <- beta_dist(.5, 1.5)
  x <- rbeta(100, 20, 20)
  o <- ef(x, d)
  print(o)
  expect_lt(o$p.value, .01)})


test_that("Beta Composite should work", {
  d <- beta_dist()
  x <- rbeta(100, 20, 20)
  o <- ef(x, d)
  print(o)
  expect_gt(o$statistic, 0)})

##### Pareto Test
test_that("Pareto: shape, scale >1", {
  # erratic
  d <- pareto_dist(3, 3)
  x <- d$sampler(100, d$par)
  o <- ef(x, d, nsim = 25)
  print(o)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})

test_that("Pareto: shape = scale = pow > 1", {
  # erratic
  expect_error(d <- pareto_dist(4, 4, 4))
})

test_that("Pareto: shape, scale >1", {
  d <- pareto_dist(5, 5)
  x <- d$sampler(100, d$par)
  o <- ef(x, d, nsim = 100)
  print(o)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})


test_that("Pareto: shape, scale <1", {
  d <- pareto_dist(.1, .1)
  x <- d$sampler(100, d$par)
  o <- ef(x, d, nsim = 0)
  print(o)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})

test_that("Pareto: shape, scale = 1", {
  d <- pareto_dist(1, 1)
  x <- d$sampler(100, d$par)
  o <- ef(x, d, nsim = 0)
  print(o)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})

test_that("Pareto: mixed", {
  d <- pareto_dist(.1, 10)
  x <- d$sampler(100, d$par)
  o <- ef(x, d, nsim = 0)
  print(o)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})

test_that("Pareto: mixed", {
  d <- pareto_dist(10, .1)
  x <- d$sampler(100, d$par)
  o <- ef(x, d, nsim = 0)
  print(o)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})


test_that("Pareto: pow > shape", {
  expect_error(pareto_dist(10, .1, pow = 5))
})


test_that("Pareto: pow > shape", {
  expect_error(pareto_dist(10, .1, 1))
})

##### Binomial tests

test_that("Binomial", {
  x <- rbinom(10, 10, .5)
  d <- binomial_dist(size = 10, prob = .5)
  o <- ef(x, d, 0)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
  x <- rbinom(10, 10, .9)
  o <- ef(x, d)
  expect_gt(o$statistic, 0)
  expect_lt(o$p.value, 0.5)
  x <- rexp(10)
  expect_error(binomial_dist(x))
})

##### Normal Tests
test_that("Normal should not be transformed", {
  x <- rnorm(10)
  o <- ef(x, normal_dist(0, 1))
  expect_identical(names(o$statistic), "E-statistic")
})

test_that("egf should return htest, even when nsim is missing", {
  x <- rnorm(10)
  o <- ef(x, normal_dist(0, 1))
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})

test_that("Normal p-vals should be uniform under Null hypothesis", {
  n <- 15
  save <- numeric(n)
  for (i in 1:n) {
    x <- rnorm(n, 0, 1)
    o <- ef(x, normal_dist(0, 1))
    save[i] <- unlist(o$p.value)
  }
  expect_gt(ef(save, uniform_dist(0, 1))$p.value, 0.05)
  expect_gt(o$statistic, 0)
})

test_that("Power to detect mean shift.", {
  n <- 10
  save <- numeric(n)
  for (i in 1:n) {
    x <- rnorm(n, 1, 1)
    o <- ef(x, normal_dist(0, 1))
    save[i] <- unlist(o$p.value)
  }
  expect_lt(mean(save[i]), 0.05)
  expect_gt(o$statistic, 0)
})

test_that("Power to detect sd shift.", {
  n <- 10
  save <- numeric(n)
  for (i in 1:n) {
    x <- rnorm(n, 0, 3)
    o <- ef(x, normal_dist(0, 1))
    save[i] <- unlist(o$p.value)
  }
  expect_lt(mean(save[i]), 0.05)
  expect_gt(o$statistic, 0)
})

test_that("Composite Test should work", {
  x <- rnorm (10)
  o <- ef(x, normal_dist(), nsim = 0)
  expect_gt(o$statistic, 0)
})
