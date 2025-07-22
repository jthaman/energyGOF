### egof: Energy goodness-of-fit tests

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
#' @param dist A string. The distribution to test.
#' @param R A positive integer. The number of parametric bootstrap replicates
#'   taken to calculate the p-value.
#' @param htest A logical. If TRUE, return an htest object, otherwise return R6
#'   EGOFTest object. The EGOFTest object has some more information about the
#'   test and a different print method.
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
#'   the ones available in the egof package.
#'
#' @return An object of class `htest' representing the result of the energy
#'   goodness-of-fit hypothesis test. The htest object has the elements:
#'
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
#' @export
#
#'
#' @details
#'
#' [TODO description of Energy GOF test here.]
#'
#' There are two types of goodness-of-fit tests covered by the egof function,
#' simple and composite. Simple GOF tests test the data `x` against a specific
#' distribution with _known parameters_ that you must pass to egof in the
#' ellipsis agrument (...). You should use a simple GOF test if you wish to
#' test questions like "my data is Normal with mean 1 and sd 2". egof can also
#' conduct _some_ composite GOF tests. A composite test is performed if no
#' parameters are passed in the ellpisis argument (...). You should conduct a
#' composite test if your research question is "my data is Normal."
#'
#' All the composite tests in egof assume that none of the parameters are
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
#' egof(x, "normal", R = 10)
#'
#' ## Simple energy goodness-of-fit test (test for Normality with known
#' ## parameters)
#'
#' egof(x, "normal", R = 10, mean = 0, sd = 1)
#'
#' ## Simple energy goodness-of-fit test for Weibull distribution
#'
#' y <- rweibull(10, 1, 1)
#' egof(y, "weibull", shape = 1, scale = 3)
#'
#' ## Error, egof does not support "partially composite" GOF tests
#'
#' ## Not run:
#' egof(x, "normal", R = 10, mean = 0)
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

#### EGOF
egof <- function(x, dist = c("uniform",
                             "exponential",
                             "bernoulli", "binomial",
                             "geometric",
                             "normal", "gaussian",
                             "beta",
                             "poisson",
                             "lognormal", "lnorm",
                             "laplace", "doubleexponential",
                             "asymmetriclaplace",
                             "inversegaussian",
                             "halfnormal",
                             "chisq", "chisquared",
                             "gamma",
                             "weibull",
                             "cauchy",
                             "pareto"),
                 R = 100, ...) {
  valid_dists <- eval(formals(egof)$dist)
  distname <- match.arg(tolower(dist), choices = valid_dists)
  validate_R(R)
  dots <- list(...)
  dist <- distribution_factory(distname, ...)
  validate_dots(dots, distname)
  validate_x(x, dist)
  egof_test(x, dist, R, ...)
}

egof_test <- function(x, dist, R, ...) {
  UseMethod("egof_test", dist)
}

#### Validation
validate_dist <- function(dist) {
  stopifnot(all(c("name", "parameter", "ref_parameter", "support", "sampler",
                  "EYY", "EXYhat") %in% names(dist)))
  stopifnot(is.logical(dist$composite_allow))
  stopifnot(setequal(names(dist$parameter), formals(dist)))
  stopifnot(setequal(names(dist$ref_parameter), formals(dist)))
}

validate_dots <- function(dots, distname) {
  dist <- distribution_factory(distname)
  required_params <- names(dist$parameter)
  supplied_params <- names(dots)
  missing_params <- setdiff(required_params, supplied_params)
  extra_params <- setdiff(supplied_params, required_params)
  no_required_params_in_dots_p <- setequal(missing_params, required_params)

  ## composite case
  if (no_required_params_in_dots_p) {
    if (dist$composite_allowed) {
      # OK
    } else {
      ## Stop if the distribution does not permit composite test.
      stop(sprintf("Cannot conduct a composite test of distribution %s.",
                   dist$name))
    }
  } else if (length(missing_params) > 0){
    ## Error if partially composite test
    stop(sprintf("Missing required parameter(s) needed for *simple* test of '%s' distribution: %s",
                 dist$name, paste(missing_params, collapse = ", ")))
  }
  ## Warning if extra stuff in ...
  if (length(extra_params) > 0) {
    warning(sprintf("Ignoring unexpected parameter(s) passed to %s distribution: %s.",
                    dist$name,
                    paste(extra_params, collapse = ", ")))
  }
}

validate_x <- function(x, dist) {
  if (any(is.na(x)) || any(is.null(x)) || any(is.infinite(x))) {
    stop ("Missing data are not supported.")
  }
  if (!dist$support(x)) {
    stop(sprintf("Not all elements of x lie in the support of distribution: %s
Support test:  %s",
dist$name, paste0(deparse(body(dist$support)),
                  collapse = "")))
  }
}

validate_R <- function(R) {
  if (!is.numeric(R))
    stop("R must be numeric.")
  if (!(R >= 0))
    stop("R must be non-negative.")
}

### Distribution Factory

distribution_factory <- function(name, ...) {
  switch(name,
         "normal" = make_normal_dist(...),
         "gaussian" = make_normal_dist(...),
         "uniform" = make_uniform_dist(...),
         "exponential" = make_exponential_dist(...),
         "beta" = make_beta_dist(...),
         "gamma" = make_gamma_dist(...),
         "weibull" = make_weibull_dist(...),
         "cauchy" = make_cauchy_dist(...),
         "pareto" = make_pareto_dist(...),
         "lognormal" = make_lognormal_dist(...),
         "lnorm" = make_lognormal_dist(...),
         "laplace" = make_laplace_dist(...),
         "doubleexponential" = make_laplace_dist(...),
         "asymmetriclaplace" = make_asymmetric_laplace_dist(...),
         "inversegaussian" = make_inverse_gaussian_dist(...),
         "standardhalfnormal" = make_halfnormal_dist(...),
         "halfnormal" = make_halfnormal_dist(...),
         "chisq" = make_chisq_dist(...),
         "chisquared" = make_chisq_dist(...),
         "binomial" = make_binomial_dist(...),
         "bernoulli" = make_bernoulli_dist(...),
         "geometric" = make_geometric_dist(...),
         "poisson" = make_poisson_dist(...),
         stop("Unsupported distribution: ", name)
         )
}


egof_test.GOFDist <- function(x, dist, R, ...) {
  composite_p <- all(sapply(dist$parameter, is.null))
  EYY <- dist$EYY(if (composite_p) dist$ref_parameter else dist$parameter)
  E_stat <- compute_E_stat(x, dist, EYY, composite_p)
  sim <- simulate_pval(x, dist, R, E_stat, composite_p)
  structure(list(
    method = paste0((if (composite_p) "Composite" else "Simple"),
                    " energy goodness-of-fit test"),
    data.name = deparse(substitute(x)),
    distribution = dist,
    parameter = c("distribution" = dist$name, (if (composite_p) NULL else dist$parameter)),
    R = R,
    composite_p = composite_p,
    statistic = E_stat,
    p.value = sim$p_value,
    sim_reps = sim$sim_reps,
    estimate = if (composite_p) lapply(dist$statistic, function(f) f(x)) else NULL
  ), class = "htest")
}

compute_E_stat <- function(x, dist, EYY, composite_p) {
  if (composite_p) x <- dist$xform(x)
  n <- length(x)
  EXYpar <- if (composite_p) dist$ref_parameter else dist$parameter
  EXY <- dist$EXYhat(x, EXYpar)
  EXX <- EXXhat(x)
  stat <- n * (2 * EXY - EYY - EXX)
  names(stat) <- paste0("E-statistic", if (composite_p) " (standardized data)" else "")
  stat
}

EXXhat <- function(x) {
  n <- length(x)
  xs <- sort(x)
  prefix <- 2 * seq_len(n) - 1 - n
  2 * mean(prefix * xs) / n
}

simulate_pval <- function(x, dist, R, E_stat, composite_p) {
  if (R == 0) return(NA)
  ran.gen.args <- if (composite_p) dist$ref_parameter else dist$parameter
  bootobj <- boot::boot(x,
                        statistic = compute_E_stat,
                        R = R,
                        sim = "parametric",
                        ran.gen = dist$sampler,
                        mle = ran.gen.args,
                        dist = dist,
                        EYY = dist$EYY(ran.gen.args),
                        composite_p = composite_p)
  list(
    sim_reps = bootobj$t,
    p_value = mean(bootobj$t > bootobj$t0)
  )
}

#### Distributions


##### Normal

make_normal_dist <- function(mean = NULL, sd = NULL) {
  structure(
    list(
      name = "Normal",
      composite_allowed = TRUE,
      parameter = list(mean = mean, sd = sd),
      ref_parameter = list(mean = 0, sd = 1),
      support = function(x) all(is.finite(x)),
      sampler = function(n, par) rnorm(n, par$mean, par$sd),
      EYY = function(par) 2 * par$sd / sqrt(pi),
      EXYhat = function(x, par) {
        mean(2 * (x - par$mean) * pnorm(x, par$mean, par$sd) +
               2 * par$sd^2 * dnorm(x, par$mean, par$sd) - (x - par$mean))
      },
      xform = function(x) (x - mean(x)) / sd(x),
      statistic = list(mean = function(x) mean(x),
                       sd = function(x) sd(x))
    ), class = c("NormalDist", "GOFDist")
  )
}


##### Uniform

make_uniform_dist <- function(min = NULL, max = NULL) {
  structure(
    list(
      name = "Uniform",
      composite_allowed = FALSE,
      parameter = list(min = min, max = max),
      ref_parameter = list(min = 0, max = 1),
      support = function(x) is.numeric(x),
      sampler = function(n, par) runif(n, par$min, par$max),
      EYY = function(par) (par$max - par$min) / 3,
      EXYhat = function(x, par = self$parameter) {
        mean((x - par$min)^2 / (par$max - par$min) - x +
               (par$max - par$min) / 2)},
      xform = function(x) (x - min(x)) / (max(x) - min(x)),
      statistic = list(min = function(x) min(x),
                       max = function(x) max(x))
    ), class = c("UniformDist", "GOFDist")
  )
}
##### Exponential
make_exponential_dist <- function(rate = NULL) {
  structure(
    list(
      name = "Exponential",
      composite_allowed = TRUE,
      parameter = list(rate = rate),
      ref_parameter = list(rate = 1),
      support = function(x) all(x > 0),
      sampler = function(n, par) rexp(n, par$rate),
      EYY = function(par) 1 / par$rate,
      EXYhat = function(x, par = self$parameter) {
        mean(x + par$rate * (1 - 2 * pexp(x, par$rate)))
      },
      xform = function(x) x / mean(x),
      statistic = list(rate = function(x) mean(x))
    ), class = c("ExponentialDist", "GOFDist")
  )
}

##### Poisson
make_poisson_dist <- function(lambda = NULL) {
  structure(
    list(
      name = "Poisson",
      composite_allowed = TRUE,
      parameter = list(lambda = lambda),
      ref_parameter = list(lambda = mean(x)),
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
    ), class = c("PoissonDist", "GOFDist")
  )
}


##### Skew-Normal?

##### Bernoulli

make_bernoulli_dist <- function(prob = NULL) {
  structure(
    list(
      name = "Bernoulli",
      composite_allowed = FALSE,
      parameter = list(prob = prob),
      ref_parameter = list(prob = NULL),
      support = function(x) all(x %in% c(0L, 1L)),
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
    ), class = c("BernoulliDist", "GOFDist")
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
make_beta_dist <- function(shape1 = NULL, shape2 = NULL) {
  structure(
    list(
      name = "Beta",
      composite_allowed = FALSE,
      parameter = list(shape1 = shape1, shape2 = shape2),
      sampler = function(n, par) {
        rbeta(n, shape1 = par$shape1, shape2 = par$shape2)},
      support = function(x) all(x < 1) && all(x > 0),
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
    ), class = c("BetaDist", "GOFDist")
  )
}

##### Dirchlet?

##### Geometric
make_geometric_dist  <- function(prob = NULL) {
  structure(
    list(
      name = "Geometric",
      composite_allowed = FALSE,
      parameter = list(prob = prob),
      support = function(x) all(x == floor(x)) && all(x > 0),
      sampler = function(n, par) rgeom(n, par$prob),
      EYY = function(p = self$parameter$prob) {
        q <- 1 - p
        (2 * q) / (1 - q^2)
      },
      EXYhat = function(x, par) {
        mean(x + 1 + (1 - 2 * pgeom(x)) / par$prob)
      }
    ), class = c("GeometricDist", "GOFDist")
  )
}


##### Negative Binomial?

##### Half-Normal
make_halfnormal_dist <- function(theta = NULL) {
  structure(
    list(
      name = "Half-Normal",
      composite_allowed = FALSE,
      parameter = list(theta = theta),
      support = function(x) all(x > 0),
      sampler = function(n, par) {
        abs(rnorm(n, 0, sd = par$theta))},
      EXYhat = function(x, par) {
        mean(2 * par$theta * (dnorm(x / par$theta) + (x / par$theta) * (pnorm(x / par$theta) - 1)))
      },
      EYY = function(par) {
        par$theta * sqrt(2) * (1 - 2 / pi)
      }
    ), class = c("HalfNormalDist", "GOFDist")
  )
}

##### Laplace
make_laplace_dist <- function(theta = NULL, sigma = NULL) {
  structure(
    list(
      name = "Laplace",
      composite_allowed = TRUE,
      parameter = list(mu = mu, sigma = sigma),
      parameter = list(mu = 0, sigma = 1),
      support = function(x) is.finite(x),
      sampler =  function(n, par) {
        u <- runif(n, -0.5, 0.5)
        par$mu - par$sigma * sign(u) * log(1 - 2 * abs(u))
      },
      EXYhat = function(x, par) {
        mean(par$sigma * exp(-abs(x - par$mu) / par$sigma) + abs(x - par$mu))
      },
      EYY = function(par) {
        2 * par$sigma
      }
    ), class = c("LaplaceDist", "GOFDist")
  )
}

##### Log-Normal
make_lognormal_dist <- function(meanlog = NULL, sdlog = NULL) {
  structure(
    list(
      name = "Log-Normal",
      composite_allowed = TRUE,
      parameter = list(meanlog = meanlog, sdlog = sdlog),
      ref_parameter = list(meanlog = 0, sdlog = 1),
      support = function(x) {
        all(x > 0)
      },
      sampler = function(n, par) {
        rlnorm(n, par$meanlog, par$sdlog)
      },
      EXYhat = function(x, par) {
        m <- par$meanlog
        s <- par$sdlog
        A <- exp(m + s^2 / 2)
        z <- (log(x) - m) / s
        z_prime <- (log(x) - m - s^2) / s
        mean(x * (2 * pnorm(z) - 1) + A * (2 * pnorm(z_prime) - 1))
      },
      EYY = function(par) {
        integrand <- function(w) {
          abs(exp(w) - 1) * dnorm(w, mean = 0, sd = sqrt(2) * par$sdlog)
        }
        scaling <- exp(par$meanlog + par$sdlog^2 / 2)
        scaling * integrate(integrand, lower = -Inf, upper = Inf)$value
      }
    ), class = c("LogNormalDist", "GOFDist")
  )
}



##### Asymmetric Laplace
make_asymmetric_laplace_dist <- function(mu = NULL, sigma = NULL, kappa = NULL) {
  structure(
    list(
      name = "Asymmetric Laplace",
      composite_allowed = TRUE,
      parameter = list(mu = mu, sigma = sigma, kappa = kappa),
      ref_parameter = list(mu = 0, sigma = 1, kappa = 1), # yes?
      support = function(x) is.finite(x),
      sampler = function(n, par = self$parameter) {
        #stuff
        #TODO
      },
      EXYhat = function(x, par) {
        mu <- (1 / par$kappa - par$kappa) / sqrt(2)
        lambda <- sqrt(2) * par$kappa / par$sigma
        beta <- sqrt(2) / (par$kappa * par$sigma)
        pk <- 1 / ( 1 + par$kappa^2)
        qk <- 1 - pk
        mean(ifelse(x >= par$theta,
                    x - par$theta - mu + (2 * pk / lambda) * exp(-lambda * abs(x - par$theta)),
                    -x + par$theta + mu + (2 * qk / beta) * exp(-beta * abs(x - par$theta))))
      },
      EYY = function(par){
        mu <- (1 / par$kappa - par$kappa) / sqrt(2)
        lambda <- sqrt(2) * par$kappa / par$sigma
        beta <- sqrt(2) / (par$kappa * par$sigma)
        pk <- 1 / (1 + par$kappa^2)
        qk <- 1 - pk
        pk / beta + qk / lambda + pk^2 / lambda + qk^2 / beta
      }
    ), class = c("AsymmetricLaplaceDist", "GOFDist")
  )
}


##### Weibull
make_weibull_dist <- function(shape = NULL, scale = NULL) {
  structure(
    list(
      name = "Weibull",
      composite_allowed = TRUE,
      parameters = list(shape = shape, scale = scale),
      parameters = list(shape = 1, scale = 1),
      support = function(x) {
        all(x > 0)
      },
      sampler = function(n, par) {
        rweibull(n, shape = par$shape, scale = par$scale)},
      EXYhat = function(x, par) {
        z = (x / par$scale)^par$shape
        mean(2 * x * pweibull(x, par$shape, par$scale) - x +
               par$scale * gamma(1 + 1 / par$shape) * (1 - 2 * pgamma(z, 1 + 1 / par$shape, 1)))
      },
      EYY = function(par) {
        # par$shape = k
        # scale = lambda
        (2 * par$scale / par$shape) * gamma(1 / par$shape) * (1 - 2^(-1 / par$shape))
      },
      xform = function(x, par) {
        (x / par$shape)^par$scale
      }
    ), class = c("WeibullDist", "GOFDist")
  )
}




##### Gamma
make_gamma_dist <- function(shape = NULL, rate = NULL) {
  structure(
    list(
      name = "gamma",
      composite_allowed = TRUE,
      parameter = list(shape = shape, rate = rate),
      ref_parameter = list(shape = 1, rate = 1),
      support = function(x) {
        all(x > 0)
      },
      sampler = function(n, par) {
        rgamma(n, shape = par$shape, rate = par$rate)},
      EXYhat = function(x, par) {
        mean(x * (2 * pgamma(x, par$shape, par$rate) - 1) +
               gamma(par$shape + 1) / (par$rate * gamma(par$shape)))
      },
      EYY = function(par) {
        2 * gamma(par$shape + 1 / 2) / (par$rate * gamma(par$shape) * sqrt(pi))
      }
    ), class = c("GammaDist", "GOFDist")
  )
}


##### Chi-Square
ChiSquaredGOFGen <- function(df = NULL) {
  structure(
    list(
      name = "Chi-Squared",
      composite_allowed = FALSE,
      parameter = list(df = df),
      parameter = list(df = NULL),
      support = function(x) {
        all(x > 0) && all(is.finite(x))
      },
      sampler = function(n, par) {
        rchisq(n, df = par$df, ncp = 0)},
      EXYhat = function(x, par) {
        mean((par$df - x) + 2 * x * pchisq(x, par$df, 0) - 2 * par$df * pchisq(x, par$df + 2, 0))
      },
      EYY = function(par) {
        4 * gamma((par$df + 1) / 2) / (sqrt(pi) * gamma(par$df / 2))
      }
    ), class = c("ChiSquaredDist", "GOFDist")
  )
}


##### Inverse Gaussian


##### Inverse Gamma?

##### Gumbel?

#### Generalized Goodness-of-fit Tests

##### Standard Cauchy
StandardCauchyGOFGen <- R6::R6Class(
  "StandardCauchyGOF",
  inherit = DistGOFGen,
  public = list(
    initialize = function(exponent = 0.5) {
      super$initialize("standardcauchy", composite_allowed = TRUE)
      self$expontent <- exponent
    },
    support = function(x) {
      is.numeric(x)
    },
    sampler = function(n) {
      rcauchy(n, 0, 1)
    },
    EXYhat = function(x, s = self$exponent) {
      mean((1 + x^2)^(s / 2) * (cos(s * arctan(x)) / cos(pi * s / 2)))
    },
    EYY = function(s = self$exponent) {
      2^s / cos(pi * s / 2)
    }
  )
)

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

## ##### Cauchy
##
## ##### Pareto
## ###### Case 1: alpha > 1, s = 1
## EXYhat.pareto.alpha.greaterthan1 <- function(x, alpha, sigma, ...) {
##   mean(y + (2 * sigma^alpha * x^(1 - alpha) - alpha * sigma) / (alpha - 1))
## }
##
## EYY.pareto.alphage.greaterthan1 <- function(alpha, sigma, ...) {
##   stopifnot(alpha > 1)
##   EY <- alpha * sigma / (alpha - 1)
##   EY / (alpha - 0.5)
## }
##
## ###### Case 2: alpha > 1, s = alpha - 1
##
## ###### Case 3: 0 < s < alpha < 1
## EXYhat.pareto.alphalessthan1<- function(x, alpha, sigma, ...) {
##
## }
##
##
## EYY.pareto.alphalessthan1 <- function(alpha, sigma, ...) {
##
## }
##
## ###### Case 4: alpha = 1 and 0 < s < 1
