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
#'
#'

### Code

egof <- function(x, dist =  c("uniform",
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
                              "standardhalfnormal", "halfnormal",
                              "chisq", "chisquared",
                              "gamma",
                              "weibull",
                              "cauchy",
                              "pareto"),
                 R = 100,
                 ...) {
  valid_dists <- eval(formals(egof)$dist)
  distname <- match.arg(tolower(dist), choices = valid_dists)
  validate_R(R)
  dots <- list(...)
  dots <- validate_dots(dots, distname)
  dist_obj <- distribution_factory(distname, ...)
  validate_x(x, dist_obj)
  test <- EGOFTest$new(x, dist = dist_obj, R = R)
  test$as_htest()
}

energygof <- egof

validate_dots <- function(dots, distname) {
  dist <- distribution_factory(distname)
  required_params <- names(dist$parameter)
  supplied_params <- names(dots)
  missing_params <- setdiff(required_params, supplied_params)
  extra_params <- setdiff(supplied_params, required_params)
  no_required_params_p <- length(setdiff(required_params, missing_params)) == 0

  ## composite case
  if (no_required_params_p) {
    # do nothing if it seems to be the composite case.
  } else if (length(missing_params) > 0){
    ## Error if partially composite test
    stop(sprintf("Missing required parameters needed for *simple* test of '%s' distribution: %s", dist$name, paste(missing_params, collapse = ", ")))
  }
  ## Warning if extra stuff in ...
  if (length(extra_params) > 0) {
    warning(sprintf("Unexpected parameters passed to '%s': %s.",
                    dist$name,
                    paste(extra_params, collapse = ", ")))
  }
  dots[names(dots) %in% required_params]
}

validate_x <- function(x, dist) {
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

distribution_factory <- function(name, ...) {
  switch(name,
         "normal" = NormalGOF$new(...),
         "gaussian" = NormalGOF$new(...),
         "uniform" = UniformGOF$new(...),
         "exponential" = ExponentialGOF$new(...),
         "bernoulli" = BinomialGOF$new(...),
         "binomial" = BinomialGOF$new(...),
         "geometric" = GeometricGOF$new(...),
         "beta" = BetaGOF$new(...),
         "poisson" = PoissonGOF$new(...),
         "lognormal" = LognormalGOF$new(...),
         "lnorm" = LognormalGOF$new(...),
         "laplace" = LaplaceGOF$new(...),
         "doubleexponential" = LaplaceGOF$new(...),
         "asymmetriclaplace" = LaplaceGOF$new(...),
         "inversegaussian" = InverseGaussianGOF$new(...),
         "halfnormal" = HalfNormalGOF$new(...),
         "chisq" = ChiSquaredGOF$new(...),
         "chisquared" = ChiSquaredGOF$new(...),
         "gamma" = GammaGOF$new(...),
         "weibull" = WeibullGOF$new(...),
         "cauchy" = CauchyGOF$new(...),
         "pareto" = ParetoGOF$new(...),
         stop("Unsupported distribution"))
}


EGOFTest <- R6::R6Class(
  "EGOFTest",
  public = list(
    dist = NULL,
    R = 0,
    x = NULL,
    composite_p = FALSE,
    E_stat = NULL,
    p_value = NULL,
    EYY = 0,

    initialize = function(x, dist, R) {
      self$x <- x
      self$R <- R
      self$dist <- dist
      self$composite_p <- all(sapply(self$dist$parameter, is.null))
      self$EYY <- self$dist$EYY(
      (if (self$composite_p)
        self$dist$ref_parameter
        else
          self$dist$parameter))
      self$E_stat <-self$compute_E_stat(x)
      self$p_value <- self$simulate_pval(x)
    },

    compute_E_stat = function(x = self$x,
                              d = self$dist,
                              EYY = self$EYY) {
      if (self$composite_p) x <- self$dist$xform(x)
      EXYpar <- (if (self$composite_p)
        self$dist$ref_parameter
        else
          self$dist$parameter)
      n <- length(x)
      EXY <- d$EXYhat(x, EXYpar)
      EXX <- self$EXXhat(x)
      out <- n * (2 * EXY - EYY - EXX)
      names(out) <- paste0("E-statistic",
      (if (self$composite_p) " (standardized data)" else ""))
      out
},

    simulate_pval = function(x = self$x, R = self$R) {
      if (self$R == 0) return(NA)
      bootobj <- boot::boot(x, statistic = self$compute_E_stat,
                            R = R, sim = "parametric",
                            ran.gen = self$dist$sampler,
                            mle =
                              (if (self$composite_p)
                                self$dist$ref_parameter
                                else
                                  self$dist$parameter),
                            EYY = self$EYY)
      mean(bootobj$t > bootobj$t0)
    },

    EXXhat = function(x = self$x) {
      n <- length(x)
      xs <- sort(x)
      prefix <- 2 * seq_len(n) - 1 - n
      2 * mean(prefix * xs) / n
    },

    as_htest = function() {
      structure(list(
        method = paste((if (self$composite_p) "Composite" else "Simple"),
                       " Energy goodness-of-fit test for",
                       self$dist$name, " distribution"),
        data.name = deparse(substitute(self$x)),
        parameters = self$dist$parameters,
        null.value = paste0(self$dist$name,
                            "Distribution with Parameters: ",
                            self$dist$parameter),
        R = self$R,
        composite_p = self$composite_p,
        statistic = self$E_stat,
        p.value = self$p_value,
        estimate = if (self$composite_p) self$dist$statistics else NULL
      ), class = "htest")
    }
  )
)

DistributionGOF <- R6::R6Class(
  "DistributionGOF",
  public = list(
    name = NULL,
    composite_allowed = FALSE,
    parameter = NULL,
    ref_parameter = NULL,
    statistic = NULL,
    initialize = function(name = NULL,
                          composite_allowed = FALSE) {
      self$name <- name
      self$composite_allowed <- composite_allowed
      self$parameter <- list()
      self$ref_parameter <- list()
      self$statistic <- list()
    },
    support = function(x) stop("Not implemented."),
    sampler = function(n, ...) stop("Not implemented."),
    EYY = function(...) stop("Not implemented."),
    EXYhat = function(x, ...) stop("Not implemented.")
  )
)


#### Distributions

##### Normal
NormalGOF <- R6::R6Class(
  "NormalGOF",
  inherit = DistributionGOF,
  public = list(
    initialize = function(mean = NULL, sd = NULL) {
      super$initialize("normal", composite_allowed = TRUE)
      self$parameter <- list(mean = mean, sd = sd)
      self$ref_parameter <- list(mean = 0, sd = 1)
      if (is.null(mean) || is.null(sd)) self$estimator(x)
    },
    estimator = function(x) {
      self$statistic <- list(mean = mean(x), sd = sd(x))
    },
    support = function(x) {
      is.numeric(x)
    },
    sampler = function(n, par = self$parameter) {
      rnorm(n, par$mean, par$sd)
    },
    EYY = function(par = self$parameter) {
      2 * par$sd / sqrt(pi)
    },
    EXYhat = function(x, par = self$parameter) {

      mean(2 * (x - par$mean) * pnorm(x, par$mean, par$sd) +
             2 * par$sd^2 * dnorm(x, par$mean, par$sd) - (x - par$mean))
    },
    xform = function(x, stat = self$statistic) {
      (x - stat$mean) / stat$sd
    }
  )
)


##### Uniform
UniformGOF <- R6::R6Class(
  "UniformGOF", inherit = DistributionGOF,
  public = list(
    initialize = function(min = NULL, max = NULL) {
      super$initialize("uniform",
                       composite_allowed = FALSE)
      self$parameter <- list(min = min, max = max)
      self$ref_parameter <- list(min = 0, max = 1)
    },
    support = function (x) {
      all(x > self$parameter$min) && all(x < self$parameter$max)},
    sampler = function(n, par = self$parameter) {
      runif(n, par$min, par$max)},
    EYY =  function(par = self$parameter) {
      (par$max - par$min) / 3},
    EXYhat = function(x, par = self$parameter) {
      mean((x - par$min)^2 / (par$max - par$min) - x +
             (par$max - par$min) / 2)
    }
  )
)

##### Exponential
ExponentialGOF <- R6::R6Class(
  "ExponentialGOF", inherit = DistributionGOF,
  public = list(
    initialize = function(rate = NULL) {
      super$initialize("exponential",
                       composite_allowed = TRUE)
      self$parameter <- list(rate = rate)
      self$ref_parameter <- list(rate = 1)
      if (is.null(rate)) self$estimator(x)
    },
    estimator = function(x) {
      self$statistic <- list(rate = 1 / mean(x))
    },
    support = function (x) all(x > 0),
    sampler = function(n, par = self$parameter) {
      rexp(n, par$rate)},
    EYY = function(par = self$parameter) {
      1 / par$rate},
    EXYhat = function(x, par = self$parameter) {
      mean(x + par$rate * (1 - 2 * pexp(x, par$rate)))
    }
  )
)

##### Poisson
PoissonGOF <- R6::R6Class(
  "PoissonGOF", inherit = DistributionGOF,
  public = list(
    initialize = function(lambda = NULL) {
      super$initialize("poisson",
                       composite_allowed = FALSE)
      self$parameter <- list(lambda = lambda)
      self$ref_parameter <- list(lambda = 1)
    },
    statistic = function(x) {
      self$statistic <- list(lambda = mean(x))
    },
    support = function (x) {
      all(x >= 0) && all(x == floor(x))},
    sampler = function(n, par = self$parameter) {
      rpois(n, par$lambda)},
    EYY = function(par = self$parameter) {
      2 * par$lambda * exp(-2 * par$lambda) * (besselI(2 * par$lambda, 0) -
                                                 besselI(2 * par$lambda, 1))
    },
    EXYhat = function(x, par = self$parameter) {
      n <- length(x)
      mean(2 * n * ppois(x, par$lambda) -
             2 * lambda * ppois(x - 1, par$lambda) + par$lambda - x)
    }
  )
)

##### Skew-Normal?

##### Bernoulli
BernoulliGOF <- R6::R6Class(
  "BernoulliGOF", inherit = DistributionGOF,
  public = list(
    initialize = function(prob = NULL) {
      super$initialize("bernoulli", composite_allowed = FALSE)
      # Set parameter values
      self$parameter$prob <- prob
    },
    support = function(x) all(x %in% c(0L, 1L)),
    sampler = function(n, par = self$parameter) rbinom(n, size = 1, prob = par$prob),
    EYY = function(prob = self$prob) 2 * prob * (1 - prob),
    EXYhat = function(x, prob = self$prob) {
      h <- sum(x)
      (h * (1 - prob) + (n - h) * prob) / n
    }
  )
)

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
BetaGOF <- R6::R6Class(
  "BetaGOF", inherit = DistributionGOF,
  public = list(
    initialize = function(shape1 = NULL, shape2 = NULL) {
      super$initialize("beta", composite_allowed = FALSE)
      # Set parameter values
      self$parameter$shape1 <- shape1
      self$parameter$shape2 <- shape2
    },
    sampler = function(n, par = self$parameter) rbeta(n, shape1 = par$shape1,
                                                      shape2 = par$shape2),
    support = function(x) all(x < 1) && all(x > 0),
    ExY = function(x, shape1 = self$parameter$shape1,
                   shape2 = self$parameter$shape2) {
      2 * x * pbeta(x, shape1, shape2) - x + (shape1 / (shape1 + shape2)) -
        2 * (beta(shape1 + 1, shape2) / beta(shape1, shape2)) *
          pbeta(x, shape1 + 1, shape2)
    },
    EYY = function(shape1 = self$parameter$shape1,
                   shape2 = self$parameter$shape2)  {
      integrate(ExY.beta, 0, 1, shape1 = shape1, shape2 = shape2)$value
    },
    EXYhat = function(x, shape1 = self$parameter$shape1,
                      shape2 = self$parameter$shape2) {
      mean(2 * x * pbeta(x, shape1, shape2) - x + (shape1 / (shape1 + shape2)) -
             2 * (beta(shape1 + 1, shape2) / beta(shape1, shape2)) *
               pbeta(x, shape1 + 1, shape2))
    }
  )
)

##### Dirchlet?

##### Geometric
GeometricGOF <- R6::R6Class(
  "GeometricGOF", inherit = DistributionGOF,
  public = list(
    initialize = function(prob = NULL) {
      super$initialize("geometric", composite_allowed = FALSE)
      # Set parameter values
      self$parameter$prob <- prob
    },
    support = function(x) all(x == floor(x)) && all(x > 0),
    sampler = function(n, par) rgeom(n, par$prob),
    EYY = function(p = self$parameter$prob) {
      q <- 1 - p
      (2 * q) / (1 - q^2)
    },
    EXYhat = function(x, p = self$parameter$prob) {
      mean(x + 1 + (1 - 2 * pgeom(x)) / p)
    }
  )
)


##### Negative Binomial?

##### Standard Half-Normal
StandardHalfNormalGOF <- R6::R6Class(
  "StandardHalfNormalGOF", inherit = DistributionGOF,
  public = list(
    initialize = function() {
      super$initialize("standardhalfnormal", composite_allowed = FALSE)
      ## No Parameters
    },
    support = function(x) all(x > 0),
    sampler = function(n, par) abs(rnorm(n)),
    EXYhat = function(x) {
      mean(4 * x * pnorm(x) + 4 * dnorm(x) - 3 * x + sqrt(2 / pi))
    },
    EYY = function(x) {
      2 * (2 - sqrt(2)) / sqrt(pi)
    }
  )
)

##### Half-Normal
HalfNormalGOF <- R6::R6Class(
  "HalfNormalGOF", inherit = DistributionGOF,
  public = list(
    initialize = function(theta = NULL) {
      super$initialize("halfnormal", composite_allowed = FALSE)
      self$parameter$theta <- theta
    },
    support = function(x) all(x > 0),
    sampler = function(n, par = self$parameter) abs(rnorm(n, 0, sd = par$theta)),
    EXYhat = function(x, theta) {
      mean(2 * theta * (dnorm(x / theta) + (x / theta) * (pnorm(x / theta) - 1)))
    },
    EYY = function(theta) {
      theta * sqrt(2) * (1 - 2 / pi)
    }
  )
)

##### Laplace
LaplaceGOF <- R6::R6Class(
  "LaplaceGOF", inherit = DistributionGOF,
  public = list(
    initialize = function(mu = NULL, sigma = NULL) {
      super$initialize("laplace", composite_allowed = TRUE)
      self$parameter$theta <- theta
      self$parameter$sigma <- sigma
    },
    support = function(x) is.numeric(x),
    sampler =  function(n, par = self$parameter) {
      mu <- par$mu; sigma <- par$sigma
      u <- runif(n, -0.5, 0.5)
      mu - sigma * sign(u) * log(1 - 2 * abs(u))
    },
    EXYhat = function(x, mu, sigma) {
      mean(sigma * exp(-abs(x - mu) / sigma) + abs(x - mu))
    },
    EYY = function(sigma) {
      2 * sigma
    }
  )
)

##### Log-Normal
LogNormalGOF <- R6::R6Class(
  "LogNormalGOF", inherit = DistributionGOF,
  public = list(
    initialize = function(meanlog = NULL, sdlog = NULL) {
      super$initialize("lognormal", composite_allowed = TRUE)
      self$parameter$meanlog <- meanlog
      self$parameter$sdlog <- sdlog
    },
    support = function(x) {
      all(x > 0)
    },
    sampler = function(n, par = self$parameter) {
      rlnorm(n, par$meanlog, par$sdlog)
    },
    EXYhat = function(x, meanlog = self$parameter$meanlog,
                      sdlog = self$parameter$sdlog) {
      A <- exp(meanlog + sdlog^2 / 2)
      z <- (log(x) - meanlog) / sdlog
      z_prime <- (log(x) - meanlog - sdlog^2) / sdlog
      mean(x * (2 * pnorm(z) - 1) + A * (2 * pnorm(z_prime) - 1))
    },
    EYY = function(meanlog = self$parameter$meanlog,
                   sdlog = self$parameter$sdlog) {
      integrand <- function(w) {
        abs(exp(w) - 1) * dnorm(w, mean = 0, sd = sqrt(2) * sigma)
      }
      scaling <- exp(meanlog + sdlog^2 / 2)
      scaling * integrate(integrand, lower = -Inf, upper = Inf)$value
    }
  )
)


##### Asymmetric Laplace
AsymmetricLaplaceGOF <- R6::R6Class(
  "AsymmetricLaplaceGOF", inherit = DistributionGOF,
  public = list(
    initialize = function(mu = NULL, sigma = NULL) {
      super$initialize("asymmetric laplace", composite_allowed = TRUE)
      self$parameter$theta <- theta
      self$parameter$sigma <- sigma
      self$parameter$kappa <- kappa
    },
    support = function(x) is.numeric(x),
    sampler = function(n, par) {
      #stuff
      #TODO
    },
    EXYhat = function(x, theta, sigma, kappa) {
      mu <- (1 / kappa - kappa) / sqrt(2)
      lambda <- sqrt(2) * kappa / sigma
      beta <- sqrt(2) / (kappa * sigma)
      pk <- 1 / ( 1 + kappa^2)
      qk <- 1 - pk
      mean(ifelse(x >= theta,
                  x - theta - mu + (2 * pk / lambda) * exp(-lambda * abs(x - theta)),
                  -x + theta + mu + (2 * qk / beta) * exp(-beta * abs(x - theta))))
    },
    EYY = function (theta, sigma, kappa){
      mu <- (1 / kappa - kappa) / sqrt(2)
      lambda <- sqrt(2) * kappa / sigma
      beta <- sqrt(2) / (kappa * sigma)
      pk <- 1 / (1 + kappa^2)
      qk <- 1 - pk
      pk / beta + qk / lambda + pk^2 / lambda + qk^2 / beta
    }
  )
)


##### Weibull
WeibullGOF <- R6::R6Class(
  "WeibullGOF", inherit = DistributionGOF,
  public = list(
    initialize = function(shape = NULL, scale = NULL) {
      super$initialize("weibull", composite_allowed = TRUE)
      self$parameter$shape <- shape
      self$parameter$scale <- scale
    },
    support = function(x) {
      all(x > 0)
    },
    sampler = function(n, par = self$parameter) rweibull(n, shape = par$shape,
                                                         scale = par$scale),
    EXYhat = function(x, shape = self$parameter$shape,
                      scale = self$parameter$scale) {
      z = (x / scale)^shape
      mean(2 * x * pweibull(x, shape, scale) - x +
             scale * gamma(1 + 1 / shape) * (1 - 2 * pgamma(z, 1 + 1 / shape, 1)))
    },
    EYY = function(shape, scale) {
      # shape = k
      # scale = lambda
      (2 * scale / shape) * gamma(1 / shape) * (1 - 2^(-1 / shape))
    }
  )
)



##### Gamma
GammaGOF <- R6::R6Class(
  "GammaGOF", inherit = DistributionGOF,
  public = list(
    initialize = function(shape = NULL, rate = NULL) {
      super$initialize("gamma", composite_allowed = TRUE)
      self$parameter$shape <- shape
      self$parameter$rate <- rate
    },
    support = function(x) {
      all(x > 0)
    },
    sampler = function(n, par = self$parameter)
      rgamma(n, shape = par$shape, rate = par$rate),
    EXYhat = function(x, shape = self$parameter$shape,
                      rate = self$parameter$rate) {
      mean(x * (2 * pgamma(x, shape, rate) - 1) +
             gamma(shape + 1) / (b * gamma(shape)))
    },
    EYY = function(shape = self$parameter$shape,
                   rate = self$parameter$rate) {
      2 * gamma(shape + 1 / 2) / (b * gamma(shape) * sqrt(pi))
    }
  )
)

##### Chi-Square
ChiSquaredGOF <- R6::R6Class(
  "ChiSquaredGOF", inherit = DistributionGOF,
  public = list(
    initialize = function(df = NULL) {
      super$initialize("chi-squared", composite_allowed = FALSE)
      self$parameter$df <- df
    },
    support = function(x) {
      all(x > 0)
    },
    sampler = function(n, par = self$parameter) rchisq(n, df = par$df,
                                                       ncp = 0),
    EXYhat = function(x, df = self$parameter$df) {
      mean((df - x) + 2 * x * pchisq(x, df, 0) - 2 * df * pchisq(x, df + 2, 0))
    },
    EYY.chisq = function(df = self$parameter$df) {
      4 * gamma((df + 1) / 2) / (sqrt(pi) * gamma(df / 2))
    }
  )
)


##### Inverse Gaussian
EXYhat.inversegaussian <- function(x, ...) {
  r <- abs(sqrt(lambda / x) * (x - mu) / mu)
  EXYhat.standardhalfnormal(r)
}


##### Inverse Gamma?

#### Generalized Goodness-of-fit Tests

##### Standard Cauchy
StandardCauchyGOF <- R6::R6Class(
  "StandardCauchyGOF",
  inherit = DistributionGOF,
  public = list(
    initialize = function(exponent = 0.5) {
      super$initialize("Standard Cauchy", composite_allowed = TRUE)
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

##### Cauchy

##### Pareto
###### Case 1: alpha > 1, s = 1
EXYhat.pareto.alpha.greaterthan1 <- function(x, alpha, sigma, ...) {
  mean(y + (2 * sigma^alpha * x^(1 - alpha) - alpha * sigma) / (alpha - 1))
}

EYY.pareto.alphage.greaterthan1 <- function(alpha, sigma, ...) {
  stopifnot(alpha > 1)
  EY <- alpha * sigma / (alpha - 1)
  EY / (alpha - 0.5)
}

###### Case 2: alpha > 1, s = alpha - 1

###### Case 3: 0 < s < alpha < 1
EXYhat.pareto.alphalessthan1<- function(x, alpha, sigma, ...) {

}


EYY.pareto.alphalessthan1 <- function(alpha, sigma, ...) {

}

###### Case 4: alpha = 1 and 0 < s < 1


## egof <- function(x,
##                  dist =  c("uniform",
##                            "exponential",
##                            "bernoulli", "binomial",
##                            "geometric",
##                            "normal", "gaussian",
##                            "beta",
##                            "poisson",
##                            "lognormal", "lnorm",
##                            "laplace", "doubleexponential",
##                            "asymmetriclaplace",
##                            "inversegaussian",
##                            "standardhalfnormal", "halfnormal",
##                            "chisq",
##                            "gamma",
##                            "weibull",
##                            "cauchy",
##                            "pareto"),
##                  ...,
##                  R = 0,
##                  s = NULL) {
##   stopifnot(is.numeric(x))
##   stopifnot(is.character(dist), length(dist) == 1)
##   stopifnot(R >= 0, length(R) == 1)
##   distname <- match.arg(dist)
##   distname <- switch(distname,
##                      "doubleexponential" = "laplace",
##                      "gaussian" = "normal",
##                      "lnorm" = "lognormal",
##                      distname)
##   n <- length(x)
##   dots <- list(...)
##   known.param <- length(dots) != 0
##   if (known.param) {
##     stopifnot(dist.data[[distname]]$param %in% names(dots))
##     EYY <- do.call(paste0("EYY.", distname), dots)
##     E.stat <- E.stat(x, distname, n, EYY, ...)
##     p.value <- simulation.pval(x, distname, n, EYY, E.stat, R, known.param, ...)
##   } else if (!known.param){
##     stopifnot(dist.data[[distname]]$composite)
##     estimates <- lapply(dist.data[[distname]]$estimator,
##                         function(estimator) estimator(x))
##     x <- dist.data[[distname]]$parameterless.dist(x)
##     x.est <- lapply(dist.data[[distname]]$estimator,
##                     function(estimator) estimator(x))
##     EYY <- do.call(paste0("EYY.", distname), x.est)
##     E.stat <- E.stat(x, distname, n, EYY, unlist(x.est))
##     p.value <- simulation.pval(x, distname, n, EYY, E.stat, R, known.param, ...)
##   }
##   out <- list(
##     method = paste0("Energy goodness-of-fit test for ",
##                     distname, " distribution",
##                     ifelse(known.param, " (Simple hypothesis test)",
##                            " (Composite hypothesis test)")),
##     data.name = paste0(deparse(substitute(x)),
##                        ", sample size ", n,
##                        ", replicates ", R),
##     statistic = E.stat,
##     #    estimate = if (!known.param) estimates else "Simple test: No parameter estimate",
##     if (!known.param) statistic = estimates,
##     parameter = dots[dist.data[[distname]]$param],
##     p.value = if (R == 0) NA else p.value
##   )
##   class(out) <- "htest"
##   out
## }


## dist.data <- list(
##   "uniform" = list(
##     composite = FALSE,
##     sampler = "runif",
##     param = c("min", "max")),
##   "exponential" = list(
##     composite = TRUE,
##     sampler = "rexp",
##     param = "rate",
##     estimator = c(rate = function(x) 1 / mean(x)),
##     parameterless.dist= function (x) x / mean(x)),
##   "bernoulli" = list(
##     composite = FALSE,
##     sampler = "rbinom",
##     sampler.arg = c(size = 1),
##     param = "prob"),
##   "binomial" = list(
##     composite = FALSE,
##     sampler = "rbinom",
##     param = c("size", "prob")),
##   "beta" = list(
##     composite = FALSE,
##     sampler = "rbeta",
##     param = c("shape1", "shape2")),
##   "chisq" = list(
##     composite = TRUE,
##     sampler = "rchisq",
##     param = c("df"),
##     sampler.arg = c(ncp = 0)),
##   "normal" = list(
##     composite = TRUE,
##     sampler = "rnorm",
##     param = c("mean", "sd"),
##     estimator = c(mean = function(x) mean(x),
##                   sd = function(x) sd(x)),
##     parameterless.dist = function (x) (x - mean(x)) / sd(x)),
##   "lognormal" = list(
##     composite = TRUE,
##     sampler = "rlnorm",
##     param = c("meanlog", "sdlog")),
##   "gamma" = list(
##     composite = FALSE,
##     sampler = "rgamma",
##     param = c("shape", "rate")),
##   "poisson" = list(
##     composite = TRUE,
##     sampler = "rpois",
##     param = "lambda"),
##   "asymmetriclaplace" = list(
##     composite = TRUE,
##     sampler = "rasymmetriclaplace",
##     param = c("theta", "sigma", "kappa")),
##   "laplace" = list(
##     composite = TRUE,
##     sampler = "rlaplace",
##     param = c("theta", "sigma")),
##   "geometric" = list(
##     composite = FALSE,
##     sampler = "rgeom",
##     param = "prob"),
##   "standardhalfnormal" = list(
##     composite = FALSE,
##     sampler = "rhnorm",
##     sampler.arg = c(theta = 1)),
##   "inversegaussian" = list(
##     composite = FALSE,
##     sampler = "rshnorm",
##     param = "theta"),
##   "weibull" = list(
##     composite = TRUE,
##     sampler = "rweibull",
##     param = c("shape", "scale"),
##     parameterless.dist = function(x) (x / scale)^shape),
##   "cauchy" = list(
##     composite = TRUE,
##     sampler = "rcauchy",
##     param = c("location", "scale")
##   )
## ) # add pareto


## E.stat <- function(x, distname, n, EYY, ...) {
##   dots <- list(...)
##   all.args <- c(as.list(environment()), dots)
##   EXYhat <-  do.call(paste0("EXYhat.", distname), all.args)
##   EXXhat <- EXXhat(x, n)
##   out <- n * (2 * EXYhat - EYY - EXXhat)
##   names(out) <- "E-statistic"
##   out
## }
##
## simulation.pval <- function(x, distname, n, EYY, E.stat, R, known.param, ...) {
##   dots <- list(...)
##   all.args <- c(as.list(environment()), dots,
##                 dist.data[[distname]]$sampler.arg)
##   sampler.name <- dist.data[[distname]]$sampler
##   Estat.nulldist <- rep(NA, R)
##   for (i in 1:R) {
##     if (known.param) {
##       samp <- do.call(sampler.name, all.args[names(formals(sampler.name))])
##       Estat.nulldist[i] <- E.stat(samp, distname, n, EYY, ...)
##     } else {
##       samp <- do.call(sampler.name, all.args[names(formals(sampler.name))])
##       samp.free <- dist.data[[distname]]$distribution.free(samp)
##       samp.free.param <- lapply(dist.data[[distname]]$estimator,
##                                 function(estimator) estimator(samp.free))
##       EYY.samp <- do.call(paste0("EYY.", distname), samp.free.param)
##       Estat.nulldist[i] <- E.stat(samp.free, distname, n,
##                                   EYY.samp, samp.free.param)
##     }
##   }
##   mean(Estat.nulldist > E.stat)
## }
##
## simulation.pval <- memoise::memoise(simulation.pval)
