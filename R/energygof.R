### egof: Energy goodness-of-fit tests

#' @title Energy goodness-of-fit tests for univariate distributions
#' @author John T. Haman
#'
#' @param x A numeric vector.
#' @param dist A string. The distribution to test.
#' @param R A positive integer. The number of parametric bootstrap replicates taken to calculate the p-value.
#' @param ... Assumed parameters of the distribution `dist'. For distributions in the R `stats' library, parameter argument names are identical. To test the composite goodness-of-fit hypothesis that x is distributed according to the family of distributions dist, don't pass parameters here.
#'#'
#' @seealso \link[stats]{Distributions} for a list of distributions available in most R installations. \link[energy]{normal.test} for the energy goodness-of-fit test with unknown parameters. \link[energy]{normal.e} for the energy goodness-of-fit statistic. See the \link[energy]{poission.mtest} for a different poisson goodness-of-fit test based on mean distances. The tests for Normal and Poisson distribution in the \link[energy] package are implemented in C/C++ , and are faster than the ones available in the egof package.
#'
#' @return An object of class `htest' representing the result of the energy goodness-of-fit hypothesis test.
#'
#' @export
#
#'
#' @details
#'
#'
#'
#' @examples
#' x <- rnorm(10)
#'
#' ## Composite energy goodness-of-fit test (Test for Normality with unknown parameters)
#' egof(x, "normal", R = 100)
#'
#' ## Simple energy goodness-of-fit test (Test for Normality with known parameters)
#' egof(x, "normal", mean = 0, sd = 1)
#'
#' @references
#'
#' Székely, G. J., & Rizzo, M. L. (2023). The energy of data and distance correlation. Chapman and Hall/CRC.
#'
#' Székely, G. J., & Rizzo, M. L. (2013). Energy statistics: A class of statistics based on distances. Journal of statistical planning and inference, 143(8), 1249-1272.
#''
#' Li, Y. (2015). Goodness-of-fit tests for Dirichlet distributions with applications. Bowling Green State University.
#'
#' Rizzo, M. L. (2002). A new rotation invariant goodness-of-fit test (PhD thesis). Bowling Green State University
#'
#' Haman, J. T. (2018). The energy goodness-of-fit test and EM type estimator for asymmetric Laplace distributions (Doctoral dissertation, Bowling Green State University).
#'
#' Ofosuhene, P. (2020). The energy goodness-of-fit test for the inverse Gaussian distribution (Doctoral dissertation, Bowling Green State University).
#'
#' Rizzo, M. L. (2009). New goodness-of-fit tests for Pareto distributions. ASTIN Bulletin: The Journal of the IAA, 39(2), 691-715.
#'
#'
#'

# User-facing wrapper
energygof <- function(x, dist =  c("uniform",
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
                      R = 100, ...) {
  dist <- match.arg(dist)
  dist_obj <- distribution_factory(dist, ...)
  test <- EGOFTest$new(x, dist = dist_obj, R = R, ...)
  test$as_htest()
}

check_known_parameters <- function(dist, dots) {
  required_params <- names(dist$parameter)
  supplied_params <- names(dots)

  # 1. Check for missing required parameters
  missing_params <- setdiff(required_params, supplied_params)
  if (length(missing_params) > 0) {
    warning(sprintf("Missing required parameters for '%s': %s",
                    dist$distribution_name,
                    paste(missing_params, collapse = ", ")))
    return(FALSE)
  }

  # 2. Check for unexpected extra parameters
  extra_params <- setdiff(supplied_params, required_params)
  if (length(extra_params) > 0) {
    warning(sprintf("Unexpected parameters passed to '%s': %s",
                    dist$distribution_name,
                    paste(extra_params, collapse = ", ")))
    # Optional: decide if you want to allow extra parameters or not.
    # If not, return FALSE here.
  }

  # 3. Check that all required parameters are not NULL
  no_nulls <- all(!vapply(dots[required_params], is.null, logical(1)))

  return(no_nulls)
}

validate_x <- function(x, dist) {
  if (!dist$support(x)) {
    stop(sprintf("Not all elements of x lie in the support of distribution: %s
Support test:  %s",
dist$distribution_name, paste0(deparse(body(dist$support)),
                               collapse = "")))
  }
}

EGOFTest <- R6::R6Class(
  "EGOFTest",
  public = list(
    dist = NULL,
    R = 0,
    known_param = FALSE,
    x = NULL,
    estimates = NULL,
    E_stat = NULL,
    p_value = NULL,

    initialize = function(x, dist, R = 0, ...) {
      validate_x(x, dist)
      dots <- list(...)
      self$x <- x
      self$R <- R
      self$dist <- dist
      self$known_param <- check_known_parameters(dist, dots)
      self$E_stat <- self$compute_E_stat(x)
      self$p_value <- self$simulate_pval(x)
    },

    compute_E_stat = function(x = self$x) {
      n <- length(self$x)
      EYY <- self$dist$EYY()
      EXY <- self$dist$EXYhat(x)
      EXX <- EXXhat(self$x)
      out <- n * (2 * EXY - EYY - EXX)
      names(out) <- "E-statistic"
      out
    },

    simulate_pval = function(x = self$x) {
      if (self$R == 0) return (NA)
      bootobj <- boot::boot(x, statistic = self$compute_E_stat,
                            R = R, sim = "parametric",
                            ran.gen = self$dist$sampler)
      1 - mean(bootobj$t < bootobj$t0)
    },

    as_htest = function() {
      structure(list(
        method = paste("Energy goodness-of-fit test for",
                       self$dist$distribution_name, " distribution"),
        data.name = deparse(substitute(self$x)),
        parameters = self$parameters,
        null.value = paste0(self$dist$distribution_name, "Distribution with Parameters: ", self$dist$parameter),
        R = self$R,
        statistic = self$E_stat,
        p.value = self$p_value,
        estimate = if (!self$known_param) self$estimates else NULL
      ), class = "htest")
    }
  )
)

EXXhat <- function(x) {
  n <- length(x)
  xs <- sort(x)
  prefix <- 2 * seq_len(n) - 1 - n
  2 * mean(prefix * xs) / n
}

# Factory function
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

DistributionGOF <- R6::R6Class(
  "DistributionGOF",
  public = list(
    distribution_name = NULL,
    composite_allowed = FALSE,
    parameter = NULL,
    statistic = NULL,
    initialize = function(distribution_name = NULL,
                          composite_allowed = FALSE) {
      self$distribution_name <- distribution_name
      self$composite_allowed <- composite_allowed
      self$parameter <- list()
      self$statistic <- list()
    },
    support = function(x) stop("Not implemented."),
    sampler = function(n, ...) stop("Not implemented."),
    EYY = function(...) stop("Not implemented."),
    EXYhat = function(x, ...) stop("Not implemented.")
  )
)


NormalGOF <- R6::R6Class(
  "NormalGOF",
  inherit = DistributionGOF,
  public = list(
    initialize = function(mean = NULL, sd = NULL) {
      super$initialize("normal", composite_allowed = TRUE)
      self$parameter$mean <- mean
      self$parameter$sd <- sd
    },
    # Statistic estimator
    statistic = function(x) {
      self$statistic$mean <- mean(x)
      self$statistic$sd <- sd(x)
    },
    support = function(x) {
      is.numeric(x)
    },
    sampler = function(n) {
      rnorm(n, self$parameter$mean, self$parameter$sd)
    },
    EYY = function() {
      2 * self$parameter$sd / sqrt(pi)
    },
    EXYhat = function(x) {
      mu <- self$parameter$mean
      s <- self$parameter$sd
      mean(2 * (x - mu) * pnorm(x, mu, s) +
             2 * s^2 * dnorm(x, mu, s) -
               (x - mu))
    },
    x_std = function(x) {
      (x - self$statistic$mean) / self$statistic$sd
    }
  )
)


##### Uniform
UniformGOF <- R6::R6Class(
  "UniformGOF", inherit = DistributionGOF,
  public = list(
    initialize = function(min = NULL, max = NULL) {
      super$initialize("uniform", composite_allowed = FALSE)
      # Set parameter values
      self$parameter$min <- min
      self$parameter$max <- max
    },
    support = function (x) all(x > self$parameter$min) && all(x < self$parameter$max),
    sampler = function(n) runif(n, self$parameter$min, self$parameter$max),
    EYY =  function() (self$parameter$max - self$parameter$min) / 3,
    EXYhat = function(x) {
      mean((x - self$parameter$min)^2 / (self$parameter$max - self$parameter$min) - x +
             (self$parameter$max - self$parameter$min) / 2)
    }
  )
)

##### Exponential
ExponentialGOF <- R6::R6Class(
  "ExponentialGOF", inherit = DistributionGOF,
  public = list(
    initialize = function(rate = NULL) {
      super$initialize("exponential", composite_allowed = TRUE)
      # Set parameter values
      self$parameter$rate <- rate
    },
    support = function (x) all(x > 0),
    sampler = function(n) rexp(n, self$parameter$rate),
    EYY = function() 1 / self$parameter$rate,
    EXYhat = function(x) {
      mean(x + self$parameter$rate * (1 - 2 * pexp(x, self$parameter$rate)))
    }
  )
)

##### Poisson
PoissonGOF <- R6::R6Class(
  "PoissonGOF", inherit = DistributionGOF,
  public = list(
    initialize = function(lambda = NULL) {
      super$initialize("poisson", composite_allowed = FALSE)
      # Set parameter values
      self$parameter$lambda <- lambda
    },
    statistic = function(x) {
      self$statistic$lambda <- mean(x)
    },
    support = function (x) all(x >= 0) && all(x == floor(x)),
    sampler = function(n) rpois(n, self$lambda),
    EYY = function(lambda = self$lambda) {
      2 * lambda * exp(-2 * lambda) * (besselI(2 * lambda, 0) - besselI(2 * lambda, 1))
    },
    EXYhat = function(x, lambda = self$lambda) {
      n <- length(x)
      mean(2 * n * ppois(x, lambda) - 2 * lambda * ppois(x - 1, lambda) + lambda - x)
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
    sampler = function(n) rbinom(n, size = 1, prob = self$parameter$prob)
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
    sampler = function(n) rbeta(n, shape1 = self$parameter$shape1,
                                shape2 = self$parameter$shape2)
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
    initialize = function(p = NULL) {
      super$initialize("beta", composite_allowed = FALSE)
      # Set parameter values
      self$parameter$p
    },
    support = function(x) all(x == floor(x)) && all(x > 0),
    EYY = function(p = self$parameter$p) {
      q <- 1 - p
      (2 * q) / (1 - q^2)
    },
    EXYhat = function(x, p = self$parameter$p) {
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
    sampler = function(n) abs(rnorm(n)),
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
    sampler = function(n) abs(rnorm(n, 0 sd = self$parameter$theta))
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
    support = function(x) is.numeric(x)
    sampler =  function(n, mu = self$parameter$mu, sigma = self$parameter$sigma) {
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
    sampler = function(n, meanlog = self$parameter$meanlog,
                       sdlog = self$parameter$sdlog) {
      rlnorm(n, mean, sd)
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
    support = function(x) is.numeric(x)
    sampler =  function(n, mu = self$parameter$mu, sigma = self$parameter$sigma) {
      u <- runif(n, -0.5, 0.5)
      mu - sigma * sign(u) * log(1 - 2 * abs(u))
    },
    sampler = function(n, theta, sigma, kappa) {
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
    sampler = function(n) rweibull(n, shape = self$parameter$shape,
                                   scale = self$parameter$scale),
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
    sampler = function(n) rgamma(n, shape = self$parameter$shape,
                                 rate = self$parameter$rate),
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
    sampler = function(n) rchisq(n, df = self$parameter$df,
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

EYY.inversegaussian <- EYY.standardhalfnormal

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
