test_that("formals same as switch", {
  b <- names(as.list(body(char_to_dist))[[2]][-c(1, 2)])
  b <- b[b != ""] # drop the ""
  f <- unlist(as.list(formals(ef.test)$dist)[-1])
  expect_setequal(b, f)
})

##### Pareto
test_that("Pareto: shape, scale >1", {
  d <- pareto_dist(5, 5, 4)
  x <- d$sampler(100, d$par)
  o <- ef(x, d)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})


test_that("Pareto: shape, scale >1", {
  d <- pareto_dist(5, 5, 1)
  x <- d$sampler(100, d$par)
  o <- ef(x, d)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})


test_that("Pareto: shape, scale <1", {
  d <- pareto_dist(.1, .1)
  x <- d$sampler(100, d$par)
  o <- ef(x, d)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})

test_that("Pareto: shape, scale = 1", {
  d <- pareto_dist(1, 1)
  x <- d$sampler(100, d$par)
  o <- ef(x, d)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})

test_that("Pareto: mixed", {
  d <- pareto_dist(.1, 10)
  x <- d$sampler(100, d$par)
  o <- ef(x, d)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})

test_that("Pareto: mixed", {
  d <- pareto_dist(10, .1)
  x <- d$sampler(100, d$par)
  o <- ef(x, d)
  expect_s3_class(o, "htest")
  expect_gt(o$statistic, 0)
})


test_that("Pareto: pow", {
  d <- pareto_dist(10, .1, pow = 5)
  x <- d$sampler(100, d$par)
  expect_error(ef(x, d))
})


test_that("Pareto: pow ", {
  d <- pareto_dist(10, .1, 1)
  x <- d$sampler(100, d$par)
  expect_error(ef(x, d))
})


##### Normal
test_that("egf should return htest, even when R is missing", {
  x <- rnorm(10)
  o <- ef(x, normal_dist(0, 1))
  expect_s3_class(o, "htest")
})
