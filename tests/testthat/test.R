set.seed(730)

##### Top Tests

test_that("formals same as switch", {
  b <- names(as.list(body(char_to_dist))[[2]][-c(1, 2)])
  b <- b[b != ""] # drop the ""
  f <- unlist(as.list(formals(ef.test)$dist)[-1])
  expect_setequal(b, f)
})


##### Check that distributions are correctly formed

test_that("Uniform", {
  expect_s3_class(uniform_dist(0, 1), "EuclideanGOFDist")
  expect_s3_class(uniform_dist(-100, 100), "EuclideanGOFDist")
  expect_error(uniform_dist(100, -100))
})

##### Pareto Test
test_that("Pareto: shape, scale >1", {
  # erratic
  d <- pareto_dist(3, 3, 2)
  x <- d$sampler(100, d$par)
  o <- ef(x, d, nsim = 0)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})

test_that("Pareto: shape = scale = pow > 1", {
  # erratic
  expect_error(d <- pareto_dist(4, 4, 4))
})

test_that("Pareto: shape, scale >1", {
  d <- pareto_dist(5, 5, 1)
  x <- d$sampler(100, d$par)
  o <- ef(x, d, nsim = 0)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})


test_that("Pareto: shape, scale <1", {
  d <- pareto_dist(.1, .1)
  x <- d$sampler(100, d$par)
  o <- ef(x, d, nsim = 0)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})

test_that("Pareto: shape, scale = 1", {
  d <- pareto_dist(1, 1)
  x <- d$sampler(100, d$par)
  o <- ef(x, d, nsim = 0)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})

test_that("Pareto: mixed", {
  d <- pareto_dist(.1, 10)
  x <- d$sampler(100, d$par)
  o <- ef(x, d, nsim = 0)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})

test_that("Pareto: mixed", {
  d <- pareto_dist(10, .1)
  x <- d$sampler(100, d$par)
  o <- ef(x, d, nsim = 0)
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
