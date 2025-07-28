test_that("formals same as switch", {
  b <- names(as.list(body(char_to_dist))[[2]][-c(1, 2)])
  b <- b[b != ""] # drop the ""
  f <- unlist(as.list(formals(egf.test)$dist)[-1])
  expect_equal(b, f)
})

test_that("egf should return htest, even when R is missing", {
  x <- rnorm(10)
  o <- egf(x, normal(0, 1))
  expect_s3_class(o, "htest")
})
