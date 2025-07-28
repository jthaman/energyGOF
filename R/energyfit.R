### energyfit: Goodness-of-fit tests via the energy of data

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' @title Energy goodness-of-fit tests for univariate distributions
#' @author John T. Haman
#'
#' @param x A numeric vector.
#' @param dist A string. The distribution to test \code{x} against.
#' @param R A positive integer. The number of parametric bootstrap replicates
#'   taken to calculate the p-value.
#' @param ... Parameters of the distribution \code{dist}. For distributions in
#'   the R `stats' library, parameter argument names are identical. To test the
#'   _composite_ goodness-of-fit hypothesis that \code{x} is distributed
#'   according to the _family of distributions_ \code{dist}, don't pass
#'   parameters in \code{...}. #'
#' @seealso \link[stats]{Distributions} for a list of distributions available
#'   in most R installations. \link[energy]{normal.test} for the energy
#'   goodness-of-fit test with unknown parameters. \link[energy]{normal.e} for
#'   the energy goodness-of-fit statistic. See the
#'   \link[energy]{poission.mtest} for a different poisson goodness-of-fit test
#'   based on mean distances. The tests for Normal and Poisson distribution in
#'   the \link[energy] package are implemented in C/C++ , and are faster than
#'   the ones available in the energyfit package.
#'
#' @return An object of class `htest' representing the result of the energy
#'   goodness-of-fit hypothesis test. The htest object has the elements:
#'
#' # TODO out of date
#' \describe{
#'   \item{\code{method}}{A test description.}
#'   \item{\code{data.name}}{Name of data passed in x.}
#'   \item{\code{parameters}}{What was passed in ..., or NULL for composite tests.}
#'   \item{\code{null-value}}{A description of the null hypothesis.}
#'   \item{\code{R}}{Number of simulation replicates.}
#'   \item{\code{composite_p}}{A logical, TRUE if composite test was performed.}
#'   \item{\code{statistic}}{NULL, or a list of MLEs calculated, if a composite test was performed.}
#'   \item{\code{p.value}}{A numeric p-value.}
#'   \item{\code{estimate}}{A numeric value of the energy statistic for testing \code{x} against {dist}.}
#' }
#'
#'
#' @details
#'
#' [TODO description of Energy GOF test here.]
#'
#' There are two types of goodness-of-fit tests covered by the energyfit function,
#' simple and composite. Simple GOF tests test the data `x` against a specific
#' distribution with _known parameters_ that you must pass to energyfit in the
#' ellipsis agrument (...). You should use a simple GOF test if you wish to
#' test questions like "my data is Normal with mean 1 and sd 2". energyfit can also
#' conduct _some_ composite GOF tests. A composite test is performed if no
#' parameters are passed in the ellpisis argument (...). You should conduct a
#' composite test if your research question is "my data is Normal."
#'
#' All the composite tests in energyfit assume that none of the parameters are
#' known. So while there is a statistical test of Normality with known mean and
#' unknown sd, this is not implemented in the energygof package. So, either
#' pass all the distribution parameters or none of them. (In the special case
#' of the Normal distribution, you can use the energy package to test the GOF
#' hypothesis with any combination of known and known parameters.)
#'
#' You should set R to be a very large number in practice. I recommend at least
#' 10,000. The default value is not a robust choice.
#'
#' @examples
#' x <- rnorm(10)
#'
#' ## Composite energy goodness-of-fit test (test for Normality with unknown
#' ## parameters)
#'
#' energyfit.test(x, "normal", R = 10)
#'
#' ## Simple energy goodness-of-fit test (test for Normality with known
#' ## parameters). ef is an alias for energyfit.test.
#'
#' ef.test(x, "normal", R = 10, mean = 0, sd = 1)
#'
#' ## Alternatively, use the energyfit generic directly so that you do not need
#' ## to pass parameter names into `...`
#'
#' energyfit(x, normal_dist(0, 1), R = 10)
#'
#' ## Simple energy goodness-of-fit test for Weibull distribution
#'
#' y <- rweibull(10, 1, 1)
#' energyfit.test(y, "weibull", shape = 1, scale = 3, R = 10)
#'
#' ## Alternatively, use the energyfit generic directly, which is slightly less
#' ## verbose. ef is an alias for energyfit.
#'
#' ef(y, weibull_dist(1, 3), R = 10
#'
#' ## energyfit does not support "partially composite" GOF tests, so this will
#' ## result in an error.
#'
#' ## Not run:
#' energyfit.test(x, "normal", mean = 0, R = 10) # sd is missing
#' ## End(Not run)
#'
#' @references
#'
#' Székely, G. J., & Rizzo, M. L. (2023). The energy of data and distance
#' correlation. Chapman and Hall/CRC.
#'
#' Székely, G. J., & Rizzo, M. L. (2013). Energy statistics: A class of
#' statistics based on distances. Journal of statistical planning and
#' inference, 143(8), 1249-1272.
#'
#' Li, Y. (2015). Goodness-of-fit tests for Dirichlet distributions with
#' applications. Bowling Green State University.
#'
#' Rizzo, M. L. (2002). A new rotation invariant goodness-of-fit test (PhD
#' thesis). Bowling Green State University
#'
#' Haman, J. T. (2018). The energy goodness-of-fit test and EM type estimator
#' for asymmetric Laplace distributions (Doctoral dissertation, Bowling Green
#' State University).
#'
#' Ofosuhene, P. (2020). The energy goodness-of-fit test for the inverse
#' Gaussian distribution (Doctoral dissertation, Bowling Green State
#' University).
#'
#' Rizzo, M. L. (2009). New goodness-of-fit tests for Pareto distributions.
#' ASTIN Bulletin: The Journal of the IAA, 39(2), 691-715.
#'
#'
#'
#' @export

### Code

#### energyfit.test (ef.test) user function
energyfit.test <- function(x, dist = c("uniform",
                                       "exponential",
                                       "bernoulli", "binomial",
                                       "geometric",
                                       "normal", "gaussian",
                                       "beta",
                                       "poisson",
                                       "lognormal", "lnorm",
                                       "laplace", "doubleexponential",
                                       "asymmetriclaplace", "alaplace",
                                       "inversegaussian",
                                       "halfnormal",
                                       "chisq", "chisquared",
                                       "gamma",
                                       "weibull",
                                       "cauchy", "stable",
                                       "pareto"),
                           R = 100,
                           ...) {
  valid_dists <- eval(formals(energyfit.test)$dist)
  distname <- match.arg(tolower(dist), choices = valid_dists)
  dots <- list(...)
  dist <- char_to_dist(distname, ...)
  energyfit(x, dist, R)
}

ef.test <- energyfit.test

#### Validation

##### Validate Parameters
validate_par <- function(dist) {
  if (!dist$parameter_domain(dist$parameter)) {
    stop(sprintf(
      "Parameters passed in ... failed domain check:  %s", paste0(deparse(body(dist$parameter_domain)), collapse = "")))
  }
}

##### Validate Distribution Obj
validate_dist <- function(dist) {
  # Not done!
  stopifnot(all(c("name", "parameter", "ref_parameter",
                  "support", "sampler",
                  "EYY", "EXYhat") %in% names(dist)))
  stopifnot(setequal(names(dist$parameter), formals(dist)))
  stopifnot(setequal(names(dist$ref_parameter), formals(dist)))
}

##### Validate Dots


##### Validate x
validate_x <- function(x, dist) {
  if (any(is.na(x)) || any(is.null(x)) || any(is.infinite(x))) {
    stop ("Missing data are not supported.")
  }
  if (!dist$support(x, dist$parameter)) {
    stop(sprintf("Not all elements of x lie in the support of distribution: %s
Support test:  %s",
dist$name, paste0(deparse(body(dist$support)),
                  collapse = "")))
  }
}

##### Validate R
validate_R <- function(R) {
  if (!is.numeric(R))
    stop("R must be numeric.")
  if (!(R >= 0))
    stop("R must be non-negative.")
  floor(R)
}

### Switchers
#### Distribution Switcher
char_to_dist <- function(name, ...) {
  ## ... should have the params and pow
  switch(name,
         ## Euclidean GOF Dists
         "normal" = normal_dist(...),
         "gaussian" = normal_dist(...),
         "uniform" = uniform_dist(...),
         "exponential" = exponential_dist(...),
         "beta" = beta_dist(...),
         "gamma" = gamma_dist(...),
         "weibull" = weibull_dist(...),
         "lognormal" = lognormal_dist(...),
         "lnorm" = lognormal_dist(...),
         "laplace" = laplace_dist(...),
         "doubleexponential" = laplace_dist(...),
         "asymmetriclaplace" = asymmetric_laplace_dist(...),
         "alaplace" = asymmetric_laplace_dist(...),
         "inversegaussian" = inverse_gaussian_dist(...),
         "halfnormal" = halfnormal_dist(...),
         "chisq" = chisq_dist(...),
         "chisquared" = chisq_dist(...),
         "binomial" = binomial_dist(...),
         "bernoulli" = bernoulli_dist(...),
         "geometric" = geometric_dist(...),
         "poisson" = poisson_dist(...),
         ## Generalized GOF Dists
         "cauchy" = cauchy_dist(...),
         "stable" = stable_dist(...),
         "pareto" = pareto_dist(...),
         stop("Unsupported distribution: ", name)
         )
}




#### energyfit_test Generic & Methods
energyfit <- function(x, dist, R = 100, ...) {
  validate_par(dist)
  validate_x(x, dist)
  R <- validate_R(R)
  UseMethod("energyfit", dist)
}

ef <- energyfit

energyfit.function <- function (x, dist, R = 100) {
  # TODO, for supplying a quantile function.
}

energyfit.GOFDist <- function(x, dist, R = 100) {
  ## Setup
  EYYpar <- if (dist$composite_p) dist$ref_parameter else dist$parameter
  ## Run functions
  EYY <- dist$EYY(EYYpar)
  E_stat <- Qhat(x, dist, EYY)
  sim <- simulate_pval(x, dist, R, E_stat)
  output_htest(x, dist, R, E_stat, sim)
}

#### Compute Energy GOF statistic
Qhat <- function(x, dist, EYY, ...) {
  UseMethod("Qhat", dist)
}

Qhat.CauchyDist <- function(x, dist, EYY) {
  x <- dist$xform(x, dist$parameter)
  NextMethod(object = dist)
}

Qhat.ParetoDist <- function(x, dist, EYY) {
  initpar <- dist$parameter
  initshape <- initpar$shape
  initscale <- initpar$scale
  initpow <- initpar$pow
  r <- initpar$r
  if (initshape > 1 && initpow != 1) {
    xshape <- initshape / r
    xpow <- if (initpow >= xshape) initpow / r else initpow
    xpar <- list(scale = initscale^r,
                 shape = xshape,
                 pow = xpow,
                 r = 1)
    # New ingredients
    x <- dist$xform(x, initpar)
    dist <- do.call("pareto_dist", xpar)
    validate_par(dist) # My fault if this is broken
    EYY <- dist$EYY(xpar)
  }
  NextMethod(object = dist)
}

Qhat.GOFDist <- function(x, dist, EYY) {
  if (dist$composite_p) {
    mle <- lapply(dist$statistic, function(f) f(x))
    x <- dist$xform(x, mle)
    EXYpar <- dist$ref_parameter
  } else {
    EXYpar <- dist$parameter
  }
  n <- length(x)
  EXY <- dist$EXYhat(x, EXYpar)
  EXX <- EXXhat(x, dist)
  E_stat <- n * (2 * EXY - EYY - EXX)
  names(E_stat) <- paste0("E-statistic",
                          if (dist$composite_p) " (standardized data)" else "")
  E_stat
}


#### EXXhat
EXXhat <- function(x, dist, ...) {
  UseMethod("EXXhat", dist)
}

EXXhat.EuclideanGOFDist <- function(x, dist) {
  n <- length(x)
  xs <- sort(x)
  prefix <- 2 * seq_len(n) - 1 - n
  2 * mean(prefix * xs) / n
}

EXXhat.GeneralizedGOFDist <- function(x, dist) {
  pow <- dist$pow
  mean(as.matrix(dist(x, "minkowski", p = pow))^pow)
}

#### Simulate P-values
simulate_pval <- function(x, dist, R, E_stat) {
  if (R == 0) return(list(sim_reps = 0, p_value = NA))
  ran.gen.args <- if (dist$composite_p) dist$ref_parameter else dist$parameter
  bootobj <- boot::boot(x,
                        statistic = Qhat,
                        R = R,
                        sim = "parametric",
                        ran.gen = dist$sampler,
                        mle = ran.gen.args,
                        dist = dist,
                        EYY = dist$EYY(ran.gen.args))
  list(
    sim_reps = bootobj$t,
    p_value = mean(bootobj$t > bootobj$t0)
  )
}


#### Output Htest
output_htest <- function(x, dist, R, E_stat, sim) {
  cp <- dist$composite_p
  if (cp) mle <- unlist(lapply(dist$statistic, function(f) f(x)))
  structure(list(
    method = paste0((if (cp) "Composite" else "Simple"),
                    " energy goodness-of-fit test"),
    data.name = deparse(substitute(x)),
    distribution = dist,
    parameter = c("distribution" = dist$name, (if (cp) NULL else dist$parameter)),
    R = R,
    pow = if (inherits(dist, "GeneralizedGOFTest")) dist$pow else NULL,
    composite_p = cp,
    statistic = E_stat,
    p.value = sim$p_value,
    sim_reps = sim$sim_reps,
    estimate = if (cp) mle  else NULL
  ), class = "htest")
}

#### Distributions
is_composite <- function(...) {
  nulls <- sapply(list(...), is.null)
  n_null <- sum(nulls)

  if (n_null == 0) {
    FALSE  # simple test: all parameters supplied
  } else if (n_null == length(nulls)) {
    TRUE   # composite test: all parameters NULL
  } else {
    stop("Partially composite tests not implemented.")
  }
}


##### Normal

normal_dist <- function(mean = NULL, sd = NULL) {
  structure(
    list(
      name = "Normal",
      composite_p = is_composite(mean, sd),
      parameter = list(mean = mean, sd = sd),
      ref_parameter = list(mean = 0, sd = 1),
      parameter_domain = function (par) {
        all(
          par$sd > 0 || is.null(par$sd),
          is.finite(par$mean) || is.null(par$mean))
      },
      support = function(x, par) all(is.finite(x)),
      sampler = function(n, par) rnorm(n, par$mean, par$sd),
      EYY = function(par) 2 * par$sd / sqrt(pi),
      EXYhat = function(x, par) {
        mean(2 * (x - par$mean) * pnorm(x, par$mean, par$sd) +
               2 * par$sd^2 * dnorm(x, par$mean, par$sd) - (x - par$mean))
      },
      xform = function(x, par) (x - par$mean) / par$sd,
      statistic = list(mean = function(x) mean(x),
                       sd = function(x) sd(x))
    ), class = c("NormalDist", "EuclideanGOFDist", "GOFDist")
  )
}


##### Uniform

uniform_dist<- function(min = 0, max = 1) {
  structure(
    list(
      name = "Uniform",
      composite_p = FALSE,
      parameter = list(min = min, max = max),
      ref_parameter = list(min = 0, max = 1),
      parameter_domain = function (par) {
        par$max - par$min > 0
      },
      support = function(x, par) is.numeric(x),
      sampler = function(n, par) runif(n, par$min, par$max),
      EYY = function(par) (par$max - par$min) / 3,
      EXYhat = function(x, par) {
        mean((x - par$min)^2 / (par$max - par$min) - x +
               (par$max - par$min) / 2)},
      xform = function(x, par) (x - par$min) / (par$max - par$min),
      statistic = list(min = function(x) min(x),
                       max = function(x) max(x))
    ), class = c("UniformDist", "EuclideanGOFDist", "GOFDist")
  )
}
##### Exponential
exponential_dist <- function(rate = NULL) {
  structure(
    list(
      name = "Exponential",
      composite_p = is.null(rate),
      parameter = list(rate = rate),
      ref_parameter = list(rate = 1),
      parameter_domain = function (par) {
        par$rate > 0
      },
      support = function(x, par) all(x > 0),
      sampler = function(n, par) rexp(n, par$rate),
      EYY = function(par) 1 / par$rate,
      EXYhat = function(x, par) {
        mean(x + par$rate * (1 - 2 * pexp(x, par$rate)))
      },
      xform = function(x, par) x / par$rate,
      statistic = list(rate = function(x) mean(x))
    ), class = c("ExponentialDist", "EuclideanGOFDist", "GOFDist")
  )
}

##### Poisson
poisson_dist <- function(lambda = NULL) {
  structure(
    list(
      name = "Poisson",
      composite_p = is.null(lambda),
      parameter = list(lambda = lambda),
      ref_parameter = list(lambda = mean(x)),
      parameter_domain = function (par) {
        par$lambda > 0 || is.null(par$lambda)
      },
      support = function (x) {
        all(x >= 0) && all(is.integer(x))},
      sampler = function(n, par) {
        rpois(n, par$lambda)},
      EYY = function(par) {
        2 * par$lambda * exp(-2 * par$lambda) * (besselI(2 * par$lambda, 0) +
                                                   besselI(2 * par$lambda, 1))
      },
      EXYhat = function(x, par) {
        n <- length(x)
        mean(2 * x * ppois(x, par$lambda) -
               2 * par$lambda * ppois(x - 1, par$lambda) + par$lambda - x)
      },
      xform = function(x) x,
      statistic = list(lambda = function(x) mean(x))
    ), class = c("PoissonDist", "EuclideanGOFDist", "GOFDist")
  )
}


##### Skew-Normal?

##### Bernoulli

bernoulli_dist <- function(prob = 0.5) {
  structure(
    list(
      name = "Bernoulli",
      composite_p = FALSE,
      parameter = list(prob = prob),
      ref_parameter = list(prob = NULL),
      parameter_domain = function (par) {
        par$prob > 0 && par$prob < 1
      },
      support = function(x, par) all(x %in% c(0L, 1L)),
      sampler = function(n, par) {
        rbinom(n, size = 1, prob = par$prob)},
      EYY = function(par) {
        2 * par$prob * (1 - par$prob)},
      EXYhat = function(x, par) {
        h <- sum(x)
        n <- length(x)
        (h * (1 - par$prob) + (n - h) * par$prob) / n
      },
      statistic = list(prob = function(x) mean(x))
    ), class = c("BernoulliDist", "EuclideanGOFDist", "GOFDist")
  )
}

##### Binomial
####### Seems to have a bug.
## EXYhat.binomial <- function(x, n, size, prob, ...) {
##   stopifnot(all(x >= 0), all(x == floor(x)))
##   k <- 0:size
##   pmf <- dbinom(k, n, prob)
##   mean(sapply(x, function(t) sum(abs(t - k) * pmf)))
## }
##
## EYY.binomial <- function(size, prob, ...) {
##   k <- 0:size
##   pmf <- dbinom(k, size, prob)
##   outer_diff <- outer(k, k, function(i, j) abs(i - j))
##   sum(outer_diff * outer(pmf, pmf))
## }

##### Beta
beta_dist <- function(shape1 = 1, shape2 = 1) {
  structure(
    list(
      name = "Beta",
      composite_p = FALSE,
      parameter = list(shape1 = shape1, shape2 = shape2),
      parameter_domain = function (par) {
        par$shape1 > 0 && par$shape2 > 0
      },
      sampler = function(n, par) {
        rbeta(n, shape1 = par$shape1, shape2 = par$shape2)},
      support = function(x, par) all(x < 1) && all(x > 0),
      EYY = function(par)  {
        integrand <- function(x, par) {
          shape1 <- par$shape1
          shape2 <- par$shape2
          2 * x * pbeta(x, shape1, shape2) - x + (shape1 / (shape1 + shape2)) -
            2 * (beta(shape1 + 1, shape2) / beta(shape1, shape2)) *
              pbeta(x, shape1 + 1, shape2)
        }
        integrate(integrand, 0, 1, par)$value
      },
      EXYhat = function(x, par) {
        mean(2 * x * pbeta(x, par$shape1, par$shape2) - x + (par$shape1 / (par$shape1 + par$shape2)) -
               2 * (beta(par$shape1 + 1, par$shape2) / beta(par$shape1, par$shape2)) *
                 pbeta(x, par$shape1 + 1, par$shape2))
      }
    ), class = c("BetaDist", "EuclideanGOFDist", "GOFDist")
  )
}

##### Dirchlet?

##### Geometric
geometric_dist  <- function(prob = 0.5) {
  structure(
    list(
      name = "Geometric",
      composite_p = FALSE,
      parameter = list(prob = prob),
      parameter_domain = function (par) {
        par$prob > 0 && par$prob < 1
      },
      support = function(x, par) all(x == floor(x)) && all(x > 0),
      sampler = function(n, par) rgeom(n, par$prob),
      EYY = function(par) {
        q <- 1 - par$prob
        (2 * q) / (1 - q^2)
      },
      EXYhat = function(x, par) {
        mean(x + 1 + (1 - 2 * pgeom(x, par$prob)) / par$prob)
      }
    ), class = c("GeometricDist", "EuclideanGOFDist", "GOFDist")
  )
}


##### Negative Binomial?

##### Half-Normal
## TODO, this seems to be bugged
halfnormal_dist <- function(scale = NULL) {
  structure(
    list(
      name = "Half-Normal",
      composite_p = is.null(scale),
      parameter = list(scale = scale),
      parameter_domain = function (par) {
        par$scale > 0 || is.null(par$scale)
      },
      support = function(x, par) all(x > 0),
      sampler = function(n, par) {
        abs(rnorm(n, 0, sd = par$scale))},
      EXYhat = function(x, par) {
        mean(2 * x * (2 * pnorm(x, 0, scale) - 1)
             - x + par$scale * sqrt(2 / pi) -
               2 * sqrt(2 / pi) * par$scale *
                 (1 - exp(-x^2 / (2 * par$scale^2))))
      },
      EYY = function(par) {
        par$scale * 2 * (2 - sqrt(2)) / sqrt(pi)
      },
      xform = function (x) x / sd(x),
      statistic = list(scale = function(x) sd(x))
    ), class = c("HalfNormalDist", "EuclideanGOFDist", "GOFDist")
  )
}

##### Laplace
laplace_dist <- function(location = NULL, scale = NULL) {
  structure(
    list(
      name = "Laplace",
      composite_p = is_composite(location, scale),
      parameter = list(location = location, scale = scale),
      ref_parameter = list(location = 0, scale = 1),
      parameter_domain = function (par) {
        par$scale > 0 || is.null(par$scale)
      },
      support = function(x, par) all(is.finite(x)),
      sampler =  function(n, par) {
        par$location + sign(runif(n) - 0.5) * rexp(n, 1 / par$scale)
      },
      EXYhat = function(x, par) {
        mean(abs(x - par$location) + par$scale
             * exp(-abs(x - par$location) / par$scale))
      },
      EYY = function(par) {
        par$scale * (3 / 2)
      },
      statistic = list(location = function(x) median(x),
                       scale = function(x) mean(abs(x - median(x)))),
      xform = function(x) {
        (x - median(x)) / mean(abs(x - median(x)))
      }
    ), class = c("LaplaceDist", "EuclideanGOFDist", "GOFDist")
  )
}

##### Log-Normal
lognormal_dist <- function(meanlog = NULL, sdlog = NULL) {
  structure(
    list(
      name = "Log-Normal",
      composite_p = is_composite(meanlog, sdlog),
      parameter = list(meanlog = meanlog, sdlog = sdlog),
      ref_parameter = list(meanlog = 0, sdlog = 1),
      support = function(x, par) {
        all(x > 0) && all(is.finite(x))
      },
      parameter_domain = function (par) {
        par$sdlog > 0 || is.null(par$sdlog)
      },
      sampler = function(n, par) {
        rlnorm(n, par$meanlog, par$sdlog)
      },
      EXYhat = function(x, par) {
        m <- par$meanlog
        s <- par$sdlog
        A <- exp(m + s^2 / 2)
        z <- (log(x) - m) / s
        w <- (m + s^2 - log(x)) / s
        mean(x * (2 * pnorm(z) - 1) + A * (2 * pnorm(w) - 1))
      },
      EYY = function(par) {
        integrand <- function(t, par) {
          m <- par$meanlog
          s <- par$sdlog
          A <- exp(m + s^2 / 2)
          z <- (log(t) - m) / s
          w <- (m + s^2 - log(t)) / s
          ExYhat <- t * (2 * pnorm(z) - 1) + A * (2 * pnorm(w) - 1)
          ExYhat * dlnorm(t, m, s)
        }
        integrate(integrand, lower = 0, upper = Inf, par = par)$value
      },
      xform = function(x) x, #todo
      statistic = list(meanlog = x, sdlog = x)
    ), class = c("LogNormalDist", "EuclideanGOFDist", "GOFDist")
  )
}



##### Asymmetric Laplace
asymmetric_laplace_dist <- function(location = NULL, scale = NULL,
                               skew = NULL) {
  structure(
    list(
      name = "Asymmetric Laplace",
      composite_p = is_composite(location, scale, skew),
      parameter = list(location = location, scale = scale, skew = skew),
      ref_parameter = list(location = 0, scale = 1, skew = 1), # yes?
      parameter_domain = function (par) {
        all(par$scale > 0 || is.null(par$scale),
            par$skew > 0 || is.null(par$skew))
      },
      support = function(x, par) all(is.finite(x)),
      sampler = function(n, par) {
        loc <- par$location
        scale <- par$scale
        k <- par$skew
        u1 <- runif(n)
        u2 <- runif(n)
        loc + scale / sqrt(2) * log(u1^k / (u2^(1 / k)))
      },
      EXYhat = function(x, par) {
        loc <- par$location
        scale <- par$scale
        k <- par$skew
        mu <- (1 / k - k) / sqrt(2)
        lam <- sqrt(2) * k / scale
        beta <- sqrt(2) / (k * scale)
        pk <- 1 / ( 1 + k^2)
        qk <- 1 - pk
        mean(ifelse(x >= loc,
                    x - loc - mu + (2 * pk / lam) * exp(-lam * abs(x - loc)),
                    -x + loc + mu + (2 * qk / beta) * exp(-beta * abs(x - loc))))
      },
      EYY = function(par){
        loc <- par$location
        scale <- par$scale
        k <- par$skew
        mu <- (1 / k - k) / sqrt(2)
        lam <- sqrt(2) * k / scale
        beta <- sqrt(2) / (k * scale)
        pk <- 1 / (1 + k^2)
        qk <- 1 - pk
        pk / beta + qk / lam + pk^2 / lam + qk^2 / beta
      },
      notes = if (composite_p)
        message("Composite Test conditional on estimation of skewness parameter.")
    ), class = c("AsymmetricLaplaceDist", "EuclideanGOFDist", "GOFDist")
  )
}

##### F??

##### Weibull
weibull_dist <- function(shape = NULL, scale = NULL) {
  structure(
    list(
      name = "Weibull",
      composite_p = is_composite(shape, scale),
      parameter = list(shape = shape, scale = scale),
      ref_parameter = list(shape = 1, scale = 1),
      support = function(x, par) {
        all(x > 0)
      },
      parameter_domain = function (par) {
        all(par$shape > 0 || is.null(par$shape),
            par$scale > 0 || is.null(par$shape))
      },
      sampler = function(n, par) {
        rweibull(n, shape = par$shape, scale = par$scale)},
      EXYhat = function(x, par) {
        z = (x / par$scale)^par$shape
        mean(2 * x * pweibull(x, par$shape, par$scale) - x +
               par$scale * gamma(1 + 1 / par$shape) *
               (1 - 2 * pgamma(z, 1 + 1 / par$shape, 1)))
      },
      EYY = function(par) {
        # par$shape = k
        # scale = lambda
        (2 * par$scale / par$shape) * gamma(1 / par$shape) *
          (1 - 2^(-1 / par$shape))
      },
      xform = function(x, par) {
        (x / par$shape)^par$scale
      }
    ), class = c("WeibullDist", "EuclideanGOFDist", "GOFDist")
  )
}

##### Gamma
gamma_dist <- function(shape = NULL, rate = NULL) {
  structure(
    list(
      name = "Gamma",
      composite_p = is_composite(shape, rate),
      parameter = list(shape = shape, rate = rate),
      ref_parameter = list(shape = 1, rate = 1),
      support = function(x, par) {
        all(x > 0) && all(is.finite(x))
      },
      parameter_domain = function (par) {
        all(
          par$shape > 0 || is.null(par$shape),
          par$rate > 0 || is.null(par$rate)
        )
      },
      sampler = function(n, par) {
        rgamma(n, shape = par$shape, rate = par$rate)},
      EXYhat = function(x, par) {
        a <- par$shape
        b <- par$rate
        mean(2 * x * pgamma(x, a, b) - x + a / b -
               2 * a / b * pgamma(x, a + 1, b))
      },
      EYY = function(par) {
        a <- par$shape
        b <- par$rate
        2 * gamma(a + 1 / 2) / (b * gamma(a) * sqrt(pi))
      }
    ), class = c("GammaDist", "EuclideanGOFDist", "GOFDist")
  )
}


##### Chi-Square
chisq_dist <- function(df = 2) {
  structure(
    list(
      name = "Chi-Squared",
      composite_p = FALSE,
      parameter = list(df = df),
      ref_parameter = list(df = NULL),
      support = function(x, par) {
        all(x > 0) && all(is.finite(x))
      },
      parameter_domain = function (par) {
        par$df > 0
      },
      sampler = function(n, par) {
        rchisq(n, df = par$df, ncp = 0)},
      EXYhat = function(x, par) {
        v <- par$df
        mean(2 * x * pchisq(x, v) - x + v -
               2 * v * pchisq(x, v + 2))
      },
      EYY = function(par) {
        v <- par$df
        4 * gamma((v + 1) / 2) / gamma(v / 2) / sqrt(pi)
      }
    ), class = c("ChiSquaredDist", "EuclideanGOFDist", "GOFDist")
  )
}


##### Inverse Gaussian
inverse_gaussian_dist <- function(mu = NULL, lambda = NULL) {
  structure(
    list(
      name = "Inverse Gaussion",
      composite_p = is_composite(mu, lambda),
      parameter = list(mu = mu, lambda = lambda),
      support = function(x, par) {
        all(x > 0) && all(is.finite(x))
      },
      parameter_domain = function (par) {
        all(par$mu > 0 || is.null(par$mu),
            par$lambda > 0 || is.null(par$mu))
      },
      sampler = function(n, par) {
        mu <- par$mu
        lam <- par$lambda
        v <- rnorm(n)^2
        x <- mu + mu^2 * y / 2 / lam - mu / 2 / lam *
          sqrt(4 * mu * lam * y + mu^2 * y^2)
        ifelse(runif(n) < mu / (mu + x), x, mu^2 / x)
      },
      EXYhat = function(x, par) {
        mu <- par$mu
        lam <- par$lambda
        A <- sqrt(lam / x) * (x / mu - 1)
        B <- exp(2 * lam / mu)
        C <- sqrt(lam / x) * (x / mu + 1)
        pinvg <- function(x, mu, lam, A, B, C) {
          pnorm(A) + B * pnorm(-C)
        }
        2 * x * pinvg(x, mu, lam, A, B, C) + mu - x - 2 *
          (mu * pnorm(A) - mu * B * pnorm(-C))
      },
      EYY = function(par) {
        mu <- par$mu
        lam <- par$lambda
        integrand <- function(t, mu, lam) {
          phi <- sqrt(mu / lam)
          erf <- function(w) 2 * pnorm(w * sqrt(2)) - 1
          8 * exp(-t^2) * erf(t) / sqrt(pi) / sqrt(t^2 + 2 * phi^2)
        }
        mu * integrate(integrand, 0, Inf, mu = mu, lam = lam)$value
      }
    ), class = c("InverseGaussianDist", "EuclideanGOFDist", "GOFDist")
  )
}


##### Inverse Gamma?

##### Gumbel?

#### Generalized Goodness-of-fit Tests

##### Pareto
pareto_dist <- function(scale = NULL, shape = NULL,
                        pow = shape / 2,
                        r = shape){
  # Use equations 8.15
  # r is only needed if shape > 1 and pow != 1.
  ## X^r ~ P(scale^r, shape / r)
  structure(
    list(
      name = "Pareto (Type I)",
      composite_p = is_composite(scale, shape),
      parameter = list(scale = scale, shape = shape,
                       pow = pow, r = r),
      ref_parameter = list(scale = 1, shape = 1,
                           pow = .5, r = 1), #??
      pow = pow,
      support = function (x, par) {
        all(x > par$scale)
      },
      parameter_domain = function (par) {
        all(
          par$scale > 0 || is.null(par$scale),
          par$shape > 0 || is.null(par$shape),
          par$pow < par$shape,
          par$shape < 1 || (par$pow == 1 || par$r >= par$shape)
        )
      },
      sampler = function(n, par) {
        u <- runif(n)
        par$scale / u^(1/par$shape)
      },
      EXYhat = function(x, par) {
        shape <- par$shape
        scale < par$scale
        pow <- par$pow
        x0 <- (x - scale) / x
        if (shape == 1) {
          A <- (x - scale)^pow
          B <- scale * pow * x^(pow - 1)
          C <- x0^pow / pow + x0^(pow + 1) / (pow + 1) *
            gsl::hyperg_2F1(1, pow + 1, pow + 2, x0)
          D <- scale * x^(pow - 1) * beta(pow + 1, 1 - pow)
          mean(A - B * C + D)
        } else if (shape > 1){
          mean(x + (2 * scale^shape * x^(1 - shape) - shape * scale) /
                 (shape - 1))
        } else {
          mean((x - scale)^pow - scale^shape *
                 (pow * beta(pow, 1 - shape) * pbeta(x0, pow, 1 - shape) -
                    shape * beta(shape - pow, pow + 1)) / x^(shape - pow))
        }
      },
      EYY = function(par) {
        shape <- par$shape
        scale <- par$scale
        pow <- par$pow
        if (shape == 1) {
          2 * scale^pow / (2 - pow) * beta(1 - pow, pow + 1)
        } else if (shape > 1) {
          2 * shape * scale / (shape - 1) / (2 * shape - 1)
        } else {
          2 * shape^2 * scale^pow * beta(shape - pow, pow + 1) / (2 * shape - pow)
        }
      },
      statistic = list(scale = function(x) min(x),
                       shape = function(x) {
                         n <- length(x)
                         n / (sum(log(x / min(x))))}),
      xform = function(x, par) {x^par$r},
      notes = if(!is.null(shape) && shape > 1 && pow != 1)
        message("Shape > 1 and pow != 1. Transforming data by data^r to conduct energy GOF test.")
    ), class = c("ParetoDist", "GeneralizedGOFDist", "GOFDist")
  )
}


##### Cauchy
cauchy_dist <- function(location = NULL, scale = NULL,
                        pow = 0.5) {
  structure(
    list(
      name = "Cauchy",
      composite_p = is_composite(location, scale),
      parameter = list(location = location, scale = scale, pow = pow),
      ref_parameter = list(location = 0, scale = 1, pow = 0.5),
      pow = pow,
      support = function(x, par) {
        all(is.finite(x))
      },
      parameter_domain = function (par) {
        all(par$scale > 0 || is.null(par$scale),
            par$pow < 1 && par$pow > 0)
      },
      sampler = function(n, par) {
        rcauchy(n, location = par$location, scale = par$scale)},
      EXYhat = function(x, par) {
        pow <- par$pow
        mean((1 + x^2)^(pow / 2) * cos(pow * atan(x)) / cospi(pow / 2))
      },
      EYY = function(par) {
        pow <- par$pow
        2^pow / cospi(pow / 2)
      },
      xform = function(x, par) {
        (x - par$location) / par$scale
      }
    ), class = c("CauchyDist", "GeneralizedGOFDist", "GOFDist")
  )
}

##### Stable
stable_dist <- function(location = NULL, scale = NULL,
                        skew = NULL, stability = NULL,
                        pow = stability / 4) {
  structure(
    list(
      name = "Stable",
      composite_p = FALSE,
      parameter = list(location = location, scale = scale, skew = skew,
                       stability = stability, pow = pow),
      ref_parameter = list(location = 0, scale = 1, skew = skew,
                           stability = stability),
      pow = pow,
      support = function(x, par) {
        if (par$skew < 1 && par$scale == 1) {
          all(x > par$location) && is.finite(x)
        } else if (par$skew < 1 && par$scale == -1) {
          all(x < par$location) && is.finite(x)
        } else
          is.finite(x)
      },
      parameter_domain = function(par) {
        all(
          par$stability > 0 && par$stability <= 2,
          is.finite(par$location),
          par$scale > 0,
          par$skew > -1 && par$stability < 1,
          par$pow < par$stability / 2
        )
      },
      sampler = function(n, par) {
        d <- par$location
        s <- par$scale
        b <- par$skew
        a <- par$stability
        u <- runif(n, -pi / 2, pi / 2)
        w <- rexp(n)
        zeta <- -b * tan(pi * a / 2)
        xi <- if (a == 1) pi / 2 else atan(-zeta) / a
        X <- if (a == 1) {
          A <- (pi / 2 + b * u) * tan(u)
          B <- b * log(pi * w / 2 * cos(u) / (pi / 2 + b * u))
          (A - B) / zeta
        } else {
          A <- (1 + zeta^2)^(1 / 2 / a)
          B <- sin(a * (u + xi)) / cos(u)^(1 / a)
          C <- (cos(u - a * (u + xi)) / w)^(1 / a - 1)
          A * B * C
        }
        if (a != 1)
          s * X + d
        else
          s * X + d + 2 * b * s * log(s) / pi
      },
      EXYhat = function(x, par) {
        n <- length(x)
        a <- par$stability
        b <- par$skew
        pow <- par$pow
        if (a == 1 && b != 0) {
          A <- 2 / pi * gamma(pow + 1)
          B <- sin(pi * pow / 2)
          integrand <- function(t, x, par, pow) {
            a <- par$stability
            b <- par$skew
            (1 - exp(-t^a) * cos(b * t^a * log(t) + x * t)) / t^(pow + 1)
          }
          I <- 0
          for (i in 1:n) {
            I <- if (x[i] > 2000) {
              I + abs(x[i])
            } else {
              I + integrate(integrand, x = x[i], par = par, pow = pow)$value
            }
            return (A * B * I / n)
          }
        } else if (a != 1 && b != 0) {
          # General case
          A <- 2 / pi * gamma(pow + 1)
          B <- sin(pi * pow / 2)
          integrand <- function(t, x, par, pow) {
            a <- par$stability
            b <- par$skew
            (1 - exp(-t^a) * cos(b * t^a * tan(pi * a / 2) - x * t)) / t^(pow + 1)
          }
          I <- 0
          for (i in 1:n) {
            I <- if (x[i] > 2000) {
              I + abs(x[i])
            } else {
              I + integrate(integrand, x = x[i], par = par, pow = pow)$value
            }
          }
          return(A * B * I / n)
        } else if (a != 1 && b == 0){
          #Symmetric Case
          integrand <- function(t, x, par, pow) {
            a <- par$stability
            (1 - exp(-t^a) * cos(x * t)) / t^(pow + 1)
          }
          I <- 0
          for (i in 1:n) {
            I <- if (x[i] > 2000) {
              I + abs(x[i])
            } else {
              I + integrate(integrand, x = x[i], par = par, pow = pow)$value
            }
          }
          return(I / n)
        } else {
          # Cauchy Case
          mean((1 + x^2)^(pow / 2) * cos(pow * atan(x)) / cos(pi * pow / 2))
        }
      },
      EYY = function(par) {
        a <- par$stability
        b <- par$skew
        pow <- par$pow
        2^(pow / a + 1) * gamma(1 - pow / a) * gamma(pow) * sin(pi * pow / 2) / pi
      },
      xform = function(x, par) {
        d <- par$location
        s <- par$scale
        a <- par$stability
        b <- par$skew
        A <- (1 / s)^(1 / a)
        B <- if (a != 1) {
          -d * (1 / s)^(1 / a - 1)
        } else {
          -d + 2 * b * log(1 / s) / pi
        }
        A * x + s * B
      }
    ), class = c("StableDist", "GeneralizedGOFDist", "GOFDist")
  )
}

#### Extras
tabular <- function(df, ...) {
  stopifnot(is.data.frame(df))

  align <- function(x) if (is.numeric(x)) "r" else "l"
  col_align <- vapply(df, align, character(1))

  cols <- lapply(df, format, ...)
  contents <- do.call("paste",
                      c(cols, list(sep = " \\tab ", collapse = "\\cr\n#'   ")))

  paste("#' \\tabular{", paste(col_align, collapse = ""), "}{\n#'   ",
        paste0("\\strong{", names(df), "}", sep = "", collapse = " \\tab "), " \\cr\n#'   ",
        contents, "\n#' }\n", sep = "")
}

deats <- data.frame(
  Distribution = character(1),
  Paramater = character(1),
  CompositeAllowed = character(1)
)

deats <- rbind(deats, list("Normal", "mean, sd", "Yes"))
deats <- rbind(deats, list("Uniform", "min, max", "No"))
deats <- rbind(deats, list("Exponential", "rate", "Yes"))
deats <- rbind(deats, list("Poisson", "lambda", "Yes"))
deats <- rbind(deats, list("Bernoulli", "prob", "No"))
deats <- rbind(deats, list("Binomial", "prob", "Yes"))
deats <- rbind(deats, list("Beta", "shape1, shape2", "No"))
deats <- rbind(deats, list("Half-Normal", "theta", "No"))
