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
#'
#'

### Code

#### EGOF
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
                 htest = TRUE,
                 ...) {
  valid_dists <- eval(formals(egof)$dist)
  distname <- match.arg(tolower(dist), choices = valid_dists)
  validate_R(R)
  dots <- list(...)
  validate_dots(dots, distname)
  dist_obj <- distribution_factory(distname, ...)
  validate_x(x, dist_obj)
  test <- EGOFTestGen$new(x, dist = dist_obj, R = R)
  if (htest) test$as_htest() else test
}

#### Validation
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
         "normal" = NormalGOFGen$new(...),
         "gaussian" = NormalGOFGen$new(...),
         "uniform" = UniformGOFGen$new(...),
         "exponential" = ExponentialGOFGen$new(...),
         "bernoulli" = BinomialGOFGen$new(...),
         "binomial" = BinomialGOFGen$new(...),
         "geometric" = GeometricGOFGen$new(...),
         "beta" = BetaGOFGen$new(...),
         "poisson" = PoissonGOFGen$new(...),
         "lognormal" = LognormalGOFGen$new(...),
         "lnorm" = LognormalGOFGen$new(...),
         "laplace" = LaplaceGOFGen$new(...),
         "doubleexponential" = LaplaceGOFGen$new(...),
         "asymmetriclaplace" = LaplaceGOFGen$new(...),
         "inversegaussian" = InverseGaussianGOFGen$new(...),
         "halfnormal" = HalfNormalGOFGen$new(...),
         "chisq" = ChiSquaredGOFGen$new(...),
         "chisquared" = ChiSquaredGOFGen$new(...),
         "gamma" = GammaGOFGen$new(...),
         "weibull" = WeibullGOFGen$new(...),
         "cauchy" = CauchyGOFGen$new(...),
         "pareto" = ParetoGOFGen$new(...),
         stop("Unsupported distribution"))
}

#### EGOFTestGen Class
EGOFTestGen <- R6::R6Class(
  "EGOFTest",
  public = list(
    dist = NULL,
    R = 0,
    sim_reps = NULL,
    x = NULL,
    sim_reps = NULL,
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
      self$E_stat <-self$compute_E_stat(self$x, self$dist, self$EYY)
      self$p_value <- self$simulate_pval(self$x, self$R)
    },

    #### Compute Energy Statistic
    compute_E_stat = function(x = self$x,
                              d = self$dist,
                              EYY = self$EYY) {
      if (self$composite_p)
        x <- self$dist$xform(x)
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

    #### Simulate Pvalue
    simulate_pval = function(x = self$x, R = self$R) {
      if (R == 0)
        return(NA)
      ran.gen.args <-
        (if (self$composite_p)
          self$dist$ref_parameter
          else
            self$dist$parameter)
      bootobj <- boot::boot(x, statistic = self$compute_E_stat,
                            R = R, sim = "parametric",
                            ran.gen = self$dist$sampler,
                            mle = ran.gen.args)
      self$sim_reps <- bootobj$t
      mean(bootobj$t > bootobj$t0)
    },

    #### EXXhat
    EXXhat = function(x = self$x) {
      n <- length(x)
      xs <- sort(x)
      prefix <- 2 * seq_len(n) - 1 - n
      2 * mean(prefix * xs) / n
    },

    #### As htest objective
    as_htest = function() {
      structure(list(
        method = paste0((if (self$composite_p) "Composite" else "Simple"),
                        " energy goodness-of-fit test"),
        data.name = deparse(substitute(x)),
        parameter = c("Distribution" = self$dist$name,
                      self$dist$parameter),
        R = self$R,
        composite_p = self$composite_p,
        statistic = self$E_stat,
        p.value = self$p_value,
        estimate = if (self$composite_p) self$dist$statistic else NULL
      ), class = "htest")
    }
  )
)

#### DistributionGOF Class
DistributionGOFGen <- R6::R6Class(
  "DistributionGOF",
  public = list(
    name = NULL,
    exponent = NULL,
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
NormalGOFGen <- R6::R6Class(
  "NormalGOF",
  inherit = DistributionGOFGen,
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
    xform = function(x) {
      (x - mean(x)) / sd(x)
    }
  )
)


##### Uniform
UniformGOFGen <- R6::R6Class(
  "UniformGOF", inherit = DistributionGOFGen,
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
ExponentialGOFGen <- R6::R6Class(
  "ExponentialGOF", inherit = DistributionGOFGen,
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
PoissonGOFGen <- R6::R6Class(
  "PoissonGOF", inherit = DistributionGOFGen,
  public = list(
    initialize = function(lambda = NULL) {
      super$initialize("poisson",
                       composite_allowed = FALSE)
      self$parameter <- list(lambda = lambda)
      self$ref_parameter <- list(lambda = 1)
      if (is.null(lambda)) self$estimator(x)
    },
    estimator = function(x) {
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
             2 * par$lambda * ppois(x - 1, par$lambda) + par$lambda - x)
    },
    xform = function(x) x
  )
)

##### Skew-Normal?

##### Bernoulli
BernoulliGOFGenGen <- R6::R6Class(
  "BernoulliGOFGen", inherit = DistributionGOF,
  public = list(
    initialize = function(prob = NULL) {
      super$initialize("bernoulli",
                       composite_allowed = FALSE)
      self$parameter$prob <- prob
    },
    support = function(x) all(x %in% c(0L, 1L)),
    sampler = function(n, par = self$parameter) {
      rbinom(n, size = 1, prob = par$prob)},
    EYY = function(prob = self$prob) {
      2 * par$prob * (1 - par$prob)},
    EXYhat = function(x, par = self$parameter) {
      h <- sum(x)
      (h * (1 - par$prob) + (n - h) * par$prob) / n
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
BetaGOFGen <- R6::R6Class(
  "BetaGOF", inherit = DistributionGOFGen,
  public = list(
    initialize = function(shape1 = NULL,
                          shape2 = NULL) {
      super$initialize("beta",
                       composite_allowed = FALSE)
      self$parameter <- list(shape1 = shape1, shape2 = shape2)
    },
    sampler = function(n, par = self$parameter) {
      rbeta(n, shape1 = par$shape1, shape2 = par$shape2)},
    support = function(x) all(x < 1) && all(x > 0),
    EYY = function(par = self$parameter)  {
      integrand <- function(x, shape1, shape2) {
        2 * x * pbeta(x, shape1, shape2) - x + (shape1 / (shape1 + shape2)) -
          2 * (beta(shape1 + 1, shape2) / beta(shape1, shape2)) *
            pbeta(x, shape1 + 1, shape2)
      }
      integrate(ExY.beta, 0, 1, shape1 = par$shape1, shape2 = par$shape2)$value
    },
    EXYhat = function(x, par = self$parameter) {
      mean(2 * x * pbeta(x, par$shape1, par$shape2) - x + (par$shape1 / (par$shape1 + par$shape2)) -
             2 * (beta(par$shape1 + 1, par$shape2) / beta(par$shape1, par$shape2)) *
               pbeta(x, par$shape1 + 1, par$shape2))
    }
  )
)

##### Dirchlet?

##### Geometric
GeometricGOFGen <- R6::R6Class(
  "GeometricGOF", inherit = DistributionGOFGen,
  public = list(
    initialize = function(prob = NULL) {
      super$initialize("geometric",
                       composite_allowed = FALSE)
      # Set parameter values
      self$parameter$prob <- prob
    },
    support = function(x) all(x == floor(x)) && all(x > 0),
    sampler = function(n, par = self$parameter) rgeom(n, par$prob),
    EYY = function(p = self$parameter$prob) {
      q <- 1 - p
      (2 * q) / (1 - q^2)
    },
    EXYhat = function(x, par = self$parameter) {
      mean(x + 1 + (1 - 2 * pgeom(x)) / par$prob)
    }
  )
)


##### Negative Binomial?

##### Standard Half-Normal
StandardHalfNormalGOFGen <- R6::R6Class(
  "StandardHalfNormalGOF", inherit = DistributionGOFGen,
  public = list(
    initialize = function() {
      super$initialize("standardhalfnormal",
                       composite_allowed = FALSE)
      ## No Parameters
    },
    support = function(x) all(x > 0),
    sampler = function(n, par) abs(rnorm(n)),
    EXYhat = function(x, par) {
      mean(4 * x * pnorm(x) + 4 * dnorm(x) - 3 * x + sqrt(2 / pi))
    },
    EYY = function(par) {
      2 * (2 - sqrt(2)) / sqrt(pi)
    }
  )
)

##### Half-Normal
HalfNormalGOFGen <- R6::R6Class(
  "HalfNormalGOF", inherit = DistributionGOFGen,
  public = list(
    initialize = function(theta = NULL) {
      super$initialize("halfnormal", composite_allowed = FALSE)
      self$parameter$theta <- theta
    },
    support = function(x) all(x > 0),
    sampler = function(n, par = self$parameter) {
      abs(rnorm(n, 0, sd = par$theta))},
    EXYhat = function(x, par = self$parameter) {
      mean(2 * par$theta * (dnorm(x / par$theta) + (x / par$theta) * (pnorm(x / par$theta) - 1)))
    },
    EYY = function(par = self$parameter) {
      par$theta * sqrt(2) * (1 - 2 / pi)
    }
  )
)

##### Laplace
LaplaceGOFGen <- R6::R6Class(
  "LaplaceGOF", inherit = DistributionGOFGen,
  public = list(
    initialize = function(mu = NULL, sigma = NULL) {
      super$initialize("laplace", composite_allowed = TRUE)
      self$parameter$theta <- theta
      self$parameter$sigma <- sigma
    },
    support = function(x) is.numeric(x),
    sampler =  function(n, par = self$parameter) {
      u <- runif(n, -0.5, 0.5)
      par$mu - par$sigma * sign(u) * log(1 - 2 * abs(u))
    },
    EXYhat = function(x, par = self$parameter) {
      mean(par$sigma * exp(-abs(x - par$mu) / par$sigma) + abs(x - par$mu))
    },
    EYY = function(par = self$parameter) {
      2 * par$sigma
    }
  )
)

##### Log-Normal
LogNormalGOFGen <- R6::R6Class(
  "LogNormalGOF", inherit = DistributionGOFGen,
  public = list(
    initialize = function(meanlog = NULL, sdlog = NULL) {
      super$initialize("lognormal",
                       composite_allowed = TRUE) # TODO
      self$parameter$meanlog <- meanlog
      self$parameter$sdlog <- sdlog
    },
    support = function(x) {
      all(x > 0)
    },
    sampler = function(n, par = self$parameter) {
      rlnorm(n, par$meanlog, par$sdlog)
    },
    EXYhat = function(x, par = self$parameter) {
      m <- par$meanlog
      s <- par$sdlog
      A <- exp(m + s^2 / 2)
      z <- (log(x) - m) / s
      z_prime <- (log(x) - m - s^2) / s
      mean(x * (2 * pnorm(z) - 1) + A * (2 * pnorm(z_prime) - 1))
    },
    EYY = function(par = self$parameter) {
      integrand <- function(w) {
        abs(exp(w) - 1) * dnorm(w, mean = 0, sd = sqrt(2) * par$sdlog)
      }
      scaling <- exp(par$meanlog + par$sdlog^2 / 2)
      scaling * integrate(integrand, lower = -Inf, upper = Inf)$value
    }
  )
)


##### Asymmetric Laplace
AsymmetricLaplaceGOFGen <- R6::R6Class(
  "AsymmetricLaplaceGOF", inherit = DistributionGOFGen,
  public = list(
    initialize = function(mu = NULL, sigma = NULL) {
      super$initialize("asymmetriclaplace",
                       composite_allowed = TRUE)
      self$parameter <- list(theta = theta, sigma = sigma, kappa = kappa)
    },
    support = function(x) is.numeric(x),
    sampler = function(n, par = self$parameter) {
      #stuff
      #TODO
    },
    EXYhat = function(x, par = self$parameter) {
      mu <- (1 / par$kappa - par$kappa) / sqrt(2)
      lambda <- sqrt(2) * par$kappa / par$sigma
      beta <- sqrt(2) / (par$kappa * par$sigma)
      pk <- 1 / ( 1 + par$kappa^2)
      qk <- 1 - pk
      mean(ifelse(x >= par$theta,
                  x - par$theta - mu + (2 * pk / lambda) * exp(-lambda * abs(x - par$theta)),
                  -x + par$theta + mu + (2 * qk / beta) * exp(-beta * abs(x - par$theta))))
    },
    EYY = function(par = self$parameter){
      mu <- (1 / par$kappa - par$kappa) / sqrt(2)
      lambda <- sqrt(2) * par$kappa / par$sigma
      beta <- sqrt(2) / (par$kappa * par$sigma)
      pk <- 1 / (1 + par$kappa^2)
      qk <- 1 - pk
      pk / beta + qk / lambda + pk^2 / lambda + qk^2 / beta
    }
  )
)


##### Weibull
WeibullGOFGen <- R6::R6Class(
  "WeibullGOF", inherit = DistributionGOFGen,
  public = list(
    initialize = function(shape = NULL, scale = NULL) {
      super$initialize("weibull",
                       composite_allowed = TRUE)
      self$parameter <- list(shape = shape, scale = scale)
    },
    support = function(x) {
      all(x > 0)
    },
    sampler = function(n, par = self$parameter) {
      rweibull(n, shape = par$shape, scale = par$scale)},
    EXYhat = function(x, par = parameter$scale) {
      z = (x / par$scale)^par$shape
      mean(2 * x * pweibull(x, par$shape, par$scale) - x +
             par$scale * gamma(1 + 1 / par$shape) * (1 - 2 * pgamma(z, 1 + 1 / par$shape, 1)))
    },
    EYY = function(par = self$parameter) {
      # par$shape = k
      # scale = lambda
      (2 * par$scale / par$shape) * gamma(1 / par$shape) * (1 - 2^(-1 / par$shape))
    }
  )
)



##### Gamma
GammaGOFGen <- R6::R6Class(
  "GammaGOF", inherit = DistributionGOFGen,
  public = list(
    initialize = function(shape = NULL, rate = NULL) {
      super$initialize("gamma",
                       composite_allowed = TRUE)
      self$parameter <- list(shape = shape, rate = rate)
    },
    support = function(x) {
      all(x > 0)
    },
    sampler = function(n, par = self$parameter) {
      rgamma(n, shape = par$shape, rate = par$rate)},
    EXYhat = function(x, par = self$parameter) {
      mean(x * (2 * pgamma(x, par$shape, par$rate) - 1) +
             gamma(par$shape + 1) / (par$rate * gamma(par$shape)))
    },
    EYY = function(par = self$parameter) {
      2 * gamma(par$shape + 1 / 2) / (par$rate * gamma(par$shape) * sqrt(pi))
    }
  )
)

##### Chi-Square
ChiSquaredGOFGen <- R6::R6Class(
  "ChiSquaredGOF", inherit = DistributionGOFGen,
  public = list(
    initialize = function(df = NULL) {
      super$initialize("chisquared",
                       composite_allowed = FALSE)
      self$parameter$df <- df
    },
    support = function(x) {
      all(x > 0)
    },
    sampler = function(n, par = self$parameter) {
      rchisq(n, df = par$df, ncp = 0)},
    EXYhat = function(x, par = self$parameter) {
      mean((par$df - x) + 2 * x * pchisq(x, par$df, 0) - 2 * par$df * pchisq(x, par$df + 2, 0))
    },
    EYY = function(par = self$parameter) {
      4 * gamma((par$df + 1) / 2) / (sqrt(pi) * gamma(par$df / 2))
    }
  )
)


##### Inverse Gaussian


##### Inverse Gamma?

#### Generalized Goodness-of-fit Tests

##### Standard Cauchy
StandardCauchyGOFGen <- R6::R6Class(
  "StandardCauchyGOF",
  inherit = DistributionGOFGen,
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
