set.seed(730)

##### Top Tests

test_that("formals same as switch", {
  b <- names(as.list(body(char_to_dist))[[2]][-c(1, 2)])
  b <- b[b != ""] # drop the ""
  f <- unlist(as.list(formals(ef.test)$dist)[-1])
  expect_setequal(b, f)
})


##### test that distributions are correctly formed

test_that("Dists", {
  ## Uniform
  expect_s3_class(uniform_dist(0, 1), "GOFDist")
  expect_s3_class(uniform_dist(-100, 100), "GOFDist")
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
  expect_s3_class(pareto_dist(NULL, NULL, 2), "GOFDist")
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

##### ef.test tests

test_that("ef.test", {
  x <- rnorm(10)
  expect_s3_class(ef.test(x, "normal"), "htest")
  expect_s3_class(ef.test(x, "laplace"), "htest")
  expect_s3_class(ef.test(x, "cauchy"), "htest")
})

##### Test that classes are correctly formed

##### Test that X is correctly validated

test_that ("validate support checks", {
  nd <- normal_dist()
  ed <- exponential_dist()
  bd <- bernoulli_dist()
  gd <- geometric_dist()
  bind <- binomial_dist()
  betad <- beta_dist()
  lnd <- lognormal_dist()
  ld <- laplace_dist()
  ald <- asymmetric_laplace_dist()
  igd <- inverse_gaussian_dist()
  hd <- halfnormal_dist()
  csd <- chisq_dist()
  gamd <- gamma_dist()
  wd <- weibull_dist()
  cd <- cauchy_dist()
  sd <- stable_dist(0, 1, -1, .5)
  pd <- pareto_dist(1, 1)

  normalx <- rnorm(10)
  posx <- rexp(10)
  intx <- c(0, rpois(10, 1))
  negx <- -posx
  x01 <- rbinom(10, 1, .5)
  unifx <- rbeta(10, 6, 6)

  ## Normal
  expect_no_error(validate_x(normalx, nd))
  expect_no_error(validate_x(posx, nd))
  expect_no_error(validate_x(intx, nd))
  expect_no_error(validate_x(negx, nd))
  expect_no_error(validate_x(x01, nd))
  expect_no_error(validate_x(unifx, nd))
  ## Exponential
  expect_error(validate_x(normalx, ed))
  expect_no_error(validate_x(posx, ed))
  expect_error(validate_x(intx, ed))
  expect_error(validate_x(negx, ed))
  expect_error(validate_x(x01, ed))
  expect_no_error(validate_x(unifx, ed))
  ## Bernoulli
  expect_error(validate_x(normalx, bd))
  expect_error(validate_x(posx, bd))
  expect_error(validate_x(intx, bd))
  expect_error(validate_x(negx, bd))
  expect_no_error(validate_x(x01, bd))
  expect_error(validate_x(unifx, bd))
  ##Geometric
  expect_error(validate_x(normalx, gd))
  expect_error(validate_x(posx, gd))
  expect_no_error(validate_x(intx, gd))
  expect_no_error(validate_x(intx[intx > 0], gd))
  expect_error(validate_x(negx, gd))
  expect_no_error(validate_x(x01, gd))
  expect_error(validate_x(unifx, gd))
  ## Binomial dist
  expect_error(validate_x(normalx, bind))
  expect_error(validate_x(posx, bind))
  expect_error(validate_x(intx, bind))
  expect_error(validate_x(negx, bind))
  expect_no_error(validate_x(x01, bind))
  expect_error(validate_x(unifx, bind))
  ## Beta Dist
  expect_error(validate_x(normalx, betad))
  expect_error(validate_x(posx, betad))
  expect_error(validate_x(intx, betad))
  expect_error(validate_x(negx, betad))
  expect_error(validate_x(x01, betad))
  expect_no_error(validate_x(unifx, betad))
  ## lognormal dist
  expect_error(validate_x(normalx, lnd))
  expect_no_error(validate_x(posx, lnd))
  expect_error(validate_x(intx, lnd))
  expect_error(validate_x(negx, lnd))
  expect_error(validate_x(x01, lnd))
  expect_no_error(validate_x(unifx, lnd))
  ## Laplace
  expect_no_error(validate_x(normalx, ld))
  expect_no_error(validate_x(posx, ld))
  expect_no_error(validate_x(intx, ld))
  expect_no_error(validate_x(negx, ld))
  expect_no_error(validate_x(x01, ld))
  expect_no_error(validate_x(unifx, ld))
  ## ALD
  expect_no_error(validate_x(normalx, ald))
  expect_no_error(validate_x(posx, ald))
  expect_no_error(validate_x(intx, ald))
  expect_no_error(validate_x(negx, ald))
  expect_no_error(validate_x(x01, ald))
  expect_no_error(validate_x(unifx, ald))
  ## IG
  expect_error(validate_x(normalx, igd))
  expect_no_error(validate_x(posx, igd))
  expect_error(validate_x(intx, igd))
  expect_error(validate_x(negx, igd))
  expect_error(validate_x(x01, igd))
  expect_no_error(validate_x(unifx, igd))
  ## Half Norm
  expect_error(validate_x(normalx, hd))
  expect_no_error(validate_x(posx, hd))
  expect_error(validate_x(intx, hd))
  expect_error(validate_x(negx, hd))
  expect_error(validate_x(x01, hd))
  expect_no_error(validate_x(unifx, hd))
  ## Chi Sq
  expect_error(validate_x(normalx, csd))
  expect_no_error(validate_x(posx, csd))
  expect_error(validate_x(intx, csd))
  expect_error(validate_x(negx, csd))
  expect_error(validate_x(x01, csd))
  expect_no_error(validate_x(unifx, csd))
  ## Gamma
  expect_error(validate_x(normalx, gamd))
  expect_no_error(validate_x(posx, gamd))
  expect_error(validate_x(intx, gamd))
  expect_error(validate_x(negx, gamd))
  expect_error(validate_x(x01, gamd))
  expect_no_error(validate_x(unifx, gamd))
  ## Weibull
  expect_error(validate_x(normalx, wd))
  expect_no_error(validate_x(posx, wd))
  expect_error(validate_x(intx, wd))
  expect_error(validate_x(negx, wd))
  expect_error(validate_x(x01, wd))
  expect_no_error(validate_x(unifx, wd))
  ## Cauchy
  expect_no_error(validate_x(normalx, cd))
  expect_no_error(validate_x(posx, cd))
  expect_no_error(validate_x(intx, cd))
  expect_no_error(validate_x(negx, cd))
  expect_no_error(validate_x(x01, cd))
  expect_no_error(validate_x(unifx, cd))
  ## Pareto
  expect_error(validate_x(normalx, pd))
  expect_error(validate_x(posx, pd))
  expect_error(validate_x(intx, pd))
  expect_error(validate_x(intx[intx > 0], pd))
  expect_error(validate_x(negx, pd))
  expect_error(validate_x(x01, pd))
  expect_error(validate_x(unifx, pd))
  ## Stable
  expect_error(validate_x(normalx, sd))
  expect_error(validate_x(posx, sd))
  expect_error(validate_x(intx, sd))
  expect_error(validate_x(intx[intx > 0], sd))
  expect_no_error(validate_x(negx, sd))
  expect_error(validate_x(x01, sd))
  expect_error(validate_x(unifx, sd))
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

##### Uniform tests

test_that("Test should work", {
  x <- runif(100)
  d <- uniform_dist(0, 1)
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_gt(o$p.value, 0.01)
})

test_that("Test should work", {
  x <- runif(100, -10, 10)
  d <- uniform_dist(-10, 10)
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_gt(o$statistic, 0.01)
})

test_that("Detect Beta", {
  x <- rbeta(100, 20, 20)
  d <- uniform_dist(0, 1)
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_lt(o$p.value, .1)
})

##### Exponential tests

test_that("Test should work", {
  x <- rexp(100)
  d <- exponential_dist(1)
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_gt(o$p.value, 0.01)
})

test_that("Test should work", {
  x <- rexp(1000, 50)
  d <- exponential_dist(50)
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_gt(o$p.value, 0.01)
})

test_that("Test should detect weibull", {
  x <- rweibull(1000, 4, 4)
  d <- exponential_dist(4)
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_lt(o$p.value, 0.01)
})

test_that("Composite should work", {
  x <- rexp(100, 4)
  d <- exponential_dist()
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_gt(o$p.value, 0.01)
})

test_that("Composite should work detect weibull", {
  x <- rweibull(100, 4, 4)
  d <- exponential_dist()
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_lt(o$p.value, 0.01)
})

##### Bernoulli tests
test_that("Test works", {
  x <- rbinom(100, 1, .5)
  d <- bernoulli_dist(.5)
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_gt(o$p.value, 0.01)
})

test_that("Detect p shift", {
  x <- rbinom(100, 1, .8)
  d <- bernoulli_dist(.5)
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_lt(o$p.value, 0.01)
})

##### Binomial tests

test_that("Test works", {
  x <- rbinom(100, 10, .5)
  d <- binomial_dist(10, .5)
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_gt(o$p.value, 0.01)
})

test_that("Detect p shift", {
  x <- rbinom(100, 10, .8)
  d <- binomial_dist(10, .5)
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_lt(o$p.value, 0.01)
})

test_that("Detect size shift", {
  x <- rbinom(100, 5, .5)
  d <- binomial_dist(10, .5)
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_lt(o$p.value, 0.01)
})

##### Geometric tests
test_that("Test works", {
  x <- rgeom(100, .5)
  d <- geometric_dist(.5)
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_gt(o$p.value, 0.01)
})

test_that("Detect p shift", {
  x <- rgeom(100, .2)
  d <- geometric_dist(.5)
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_lt(o$p.value, 0.01)
})

##### Poisson tests
test_that("Test works", {
  x <- rpois(100, 10)
  d <- poisson_dist(10)
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_gt(o$p.value, 0.01)
})

## Something is wrong
test_that("Detect lam shift", {
  x <- rpois(100,5)
  d <- poisson_dist(10)
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_lt(o$p.value, 0.01)
})

test_that("Composite Test works", {
  x <- rpois(100, 10)
  d <- poisson_dist()
  o <- ef(x, d, nsim = 60)
  print(o)
  expect_gt(o$statistic, 0)
  expect_gt(o$p.value, 0.01)
})


##### Asym Laplace tests
##### Inv Gaussian
##### Half Norm tests
##### Chi Sq tests
##### Gamma  tests
##### Weibull  tests
##### Cauchy tests
##### Stable tests

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

test_that("Beta Composite should work", {
  d <- beta_dist()
  x <- rbeta(100, 1 / 60, 70)
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


##### Lognormal Tests

test_that("lnorm test", {
  x <- rlnorm(1000)
  o <- ef(x, lognormal_dist(0, 1), nsim = 25)
  print(o)
  expect_gt(o$statistic, 0)
})

test_that("lnorm test", {
  x <- rlnorm(1000, 5, 5)
  d <- lognormal_dist(5, 5)
  o <- ef(x, d, nsim = 25)
  print(o)
  expect_gt(o$statistic, 0)
  expect_gt(o$p.value, 0.001)
})

test_that("lnorm test should detect meanlog shift", {
  x <- rlnorm(1000, 7, 5)
  d <- lognormal_dist(5, 5)
  o <- ef(x, d, nsim = 25)
  print(o)
  expect_gt(o$statistic, 0)
  expect_lt(o$p.value, 0.5)
})

test_that("lnorm composite test ", {
  x <- rlnorm(1000, 2, 2)
  d <- lognormal_dist()
  o <- ef(x, d, nsim = 25)
  print(o)
  expect_gt(o$statistic, 0)
  expect_gt(o$p.value, 0.01)
})

test_that("lnorm composite test can detect weibull ", {
  x <- rweibull(1000, 2, 2)
  d <- lognormal_dist()
  o <- ef(x, d, nsim = 25)
  print(o)
  expect_gt(o$statistic, 0)
  expect_lt(o$p.value, 0.05)
})

##### Laplace Tests

test_that("laplace test", {
  d <- laplace_dist(0, 1)
  x <- d$sampler(1000, d$par)
  o <- ef(x, d, nsim = 25)
  print(o)
  expect_gt(o$statistic, 0)
})

test_that("laplace test, different params. ", {
  d <- laplace_dist(5, 5)
  x <- d$sampler(1000, d$par)
  o <- ef(x, d, nsim = 25)
  print(o)
  expect_gt(o$statistic, 0)
  expect_gt(o$p.value, 0.001)
})

test_that("laplace test should detect location shift", {
  d <- laplace_dist(5, 5)
  x <- d$sampler(1000, list(location = 0, scale = 5))
  o <- ef(x, d, nsim = 25)
  print(o)
  expect_gt(o$statistic, 0)
  expect_lt(o$p.value, 0.05)
})

test_that("laplace composite test ", {
  d <- laplace_dist()
  x <- d$sampler(1000, list(location = 0, scale = 5))
  o <- ef(x, d, nsim = 25)
  print(o)
  expect_gt(o$statistic, 0)
  expect_gt(o$p.value, 0.01)
})

test_that("laplace composite test can detect normal ", {
  d <- laplace_dist()
  x <- rnorm(1000)
  o <- ef(x, d, nsim = 25)
  print(o)
  expect_gt(o$statistic, 0)
  expect_lt(o$p.value, 0.01)
})

#####
