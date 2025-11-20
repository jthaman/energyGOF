### energyGOF: Goodness-of-fit tests for univariate data via energy

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

#' @title Goodness-of-fit tests for univariate data via energy
#' @author John T. Haman
#' @description Perform a goodness-of-fit test of univariate data `x` against a
#'   target `y`. `y` may be one of the following:
#'
#'   * A string naming a distribution. For example, "normal". Both simple
#'   (known parameter) and composite (unknown parameter) tests are supported,
#'   but not all distributions allow for a composite test. See
#'   [energyGOF-package] for the table of supported distributions.
#'
#'     * Result: A parametric goodness-of-fit test is performed.
#'     * Allowable values: "uniform", "exponential", "bernoulli", "binomial", "geometric",
#'     "normal" ("gaussian"), "beta", "poisson", "lognormal" ("lnorm"),
#'     "laplace" ("doubleexponential"), "asymmetriclaplace" ("alaplace"),
#'     "inversegaussian", ("invgaussian"), "halfnormal", "chisq" ( "chisquared"),
#'     "F", "gamma", "weibull", "cauchy", "pareto".
#'
#'   * A numeric vector of data.
#'
#'     * Result: A two-sample, non-parametric goodness-of-fit test is performed
#'       to test if x and y are equally distributed.
#'
#'   * A **continuous** cumulative distribution function. For example, `pt`.
#'   Only simple tests are supported.
#'
#'     * Result: \eqn{y(x)} is tested for uniformity.
#'
#'
#'   *P*-values are determined via parametric bootstrap. For distributions
#'   where \eqn{E|Y|} is not finite (Cauchy, Pareto), a *generalized* energy
#'   goodness-of-fit test is performed, and an additional tuning parameter
#'   `pow` is required.
#' @param x A numeric vector.
#' @param y A string, distribution function, or numeric vector. The
#'   distribution to test `x` against.
#' @param nsim A non-negative integer. The number of parametric bootstrap
#'   replicates taken to calculate the *p*-value. If 0, no simulation.
#' @param ... If `y` is a string or distribution function, parameters of the
#'   distribution `y`. Required for a simple test. For distributions in the
#'   [stats] library, parameter argument names are identical. If `y` is a
#'   string, to test the *composite* goodness-of-fit hypothesis that `x` is
#'   distributed according to the *family of distributions* `y`, don't pass
#'   parameters in `...`. For *generalized* energy tests, you can also
#'   optionally pass the generalized energy exponent `pow` here. Composite
#'   testing is not supported if `y` is a function. (As you can see, there is a
#'   lot going on in `...` and if you don't like that, you may want to check
#'   out [energyGOFdist] for a structured interface.)
#'
#' @seealso
#'
#'  * [energyGOF-package] for specifics on the distributions available to test.
#'
#'  * [energyGOFdist] for the alternate S3 interface for parametric testing.
#'
#'  * \link[stats]{Distributions} for a list of distributions available
#'   in most R installations.
#'
#'  * [energy::eqdist.etest()] for information on the two-sample test.
#'
#'  * [energy::normal.test()] for the energy goodness-of-fit test with unknown
#' parameters. The tests for (multivariate) Normal in the energy package are
#' implemented with compiled code, and are faster than the one available in the
#' energyGOF package.
#'
#'  * [energy::poisson.mtest()] for a different Poisson goodness-of-fit test
#'   based on mean distances.
#'
#' @return If `y` is a string or function, return an object of class `htest'
#'   representing the result of the energy goodness-of-fit hypothesis test. The
#'   htest object has the elements:
#'
#' * `method`: Simple or Composite
#' * `data.name`
#' * `distribution`: The distribution object created to test
#' * `parameter`: List of parameters if the test is simple
#' * `nsim`: Number of bootstrap replicates
#' * `composite_p`: TRUE/FALSE composite predicate
#' * `statistic`: The value of the energy statistic (\eqn{Q=nE^*})
#' * `p.value`
#' * `sim_reps`: bootstrap replicates of energy statistic
#' * `estimate`: Any parameter estimates, if the test is composite
#'
#' If `y` is numeric, return the same htest object as [energy::eqdist.etest].
#'
#' @aliases egof.test
#'
#' @examples
#' x <- rnorm(10)
#' y <- rt(10, 4)
#'
#' ## Composite energy goodness-of-fit test (test for Normality with unknown
#' ## parameters)
#'
#' energyGOF.test(x, "normal", nsim = 10)
#'
#' ## Simple energy goodness-of-fit test (test for Normality with known
#' ## parameters). egof.test is an alias for energyGOF.test.
#'
#' egof.test(x, "normal", nsim = 10, mean = 0, sd = 1)
#'
#' ## Alternatively, use the energyGOFdist generic directly so that you do not need
#' ## to pass parameter names into `...`
#'
#' energyGOFdist(x, normal_dist(0, 1), nsim = 10)
#'
#' ## Conduct a two-sample test
#'
#' egof.test(x, y, 0)
#'
#' ## Conduct a test against any continuous distribution function
#'
#' egof.test(x, pcauchy, 0)
#'
#' ## Simple energy goodness-of-fit test for Weibull distribution
#'
#' y <- rweibull(10, 1, 1)
#' energyGOF.test(y, "weibull", shape = 1, scale = 3, nsim = 10)
#'
#' ## Alternatively, use the energyGOFdist generic directly, which is slightly less
#' ## verbose. egofd is an alias for energyGOFdist.
#'
#' egofd(y, weibull_dist(1, 3), nsim = 10)
#'
#' ## Conduct a generalized GOF test. `pow` is the exponent *s* in the generalized
#' ## energy statistic. Pow is only necessary when testing Cauchy, and
#' ## Pareto distributions. If you don't set a pow, there is a default for each
#' ## of the distributions, but the default isn't necessarily better than any
#' ## other number.
#'
#' egofd(rcauchy(100),
#'    cauchy_dist(location = 0, scale = 1, pow = 0.5),
#'    nsim = 10)
#'
#' ## energyGOF does not support tests with a mix of known and unknown
#' ## parameters, so this will result in an error.
#'
#' \dontrun{
#'   energyGOF.test(x, "normal", mean = 0, nsim = 10) # sd is missing
#' }
#'
#'
#' @importFrom stats dlnorm dnorm integrate median pbeta pchisq pexp pgamma
#' pgeom pnorm ppois pweibull rbeta rbinom rcauchy rchisq rexp rgamma rgeom
#' rlnorm rnorm rpois runif rweibull sd dbeta dbinom mahalanobis pf rf setNames
#'
#' @importFrom utils lsf.str
#'
#' @export energyGOF.test

energyGOF.test <- function(x, y, nsim, ...) {
  UseMethod("energyGOF.test", y)
}

#' @rdname energyGOF.test
#' @export egof.test
egof.test <- energyGOF.test

#' @export
energyGOF.test.function <- function(x, y, nsim, ...) {
  nsim <- validate_nsim(nsim)
  validate_cdf(y, x, ...)
  args <- c(list(x), list(...))
  first_arg <- names(formals(y))[1]
  args <- c(setNames(list(x), first_arg), list(...))
  xu <- do.call(y, args)
  d <- uniform_dist(0, 1)
  egofd(xu, d, nsim = nsim)
}

#' @export
energyGOF.test.numeric <- function(x, y, nsim, ...) {
  nsim <- validate_nsim(nsim)
  pooled <- matrix(c(x, y))
  sizes <- c(NROW(x), NROW(y))
  energy::eqdist.etest(pooled, sizes = sizes, method = "original", R = nsim)
}

#' @export
energyGOF.test.character <- function(
  x,
  y = c(
    "uniform",
    "exponential",
    "bernoulli",
    "binomial",
    "geometric",
    "normal",
    "gaussian",
    "beta",
    "poisson",
    "lognormal",
    "lnorm",
    "laplace",
    "doubleexponential",
    "asymmetriclaplace",
    "alaplace",
    "inversegaussian",
    "invgaussian",
    "halfnormal",
    "chisq",
    "chisquared",
    "f",
    "gamma",
    "weibull",
    "cauchy",
    "pareto"
  ),
  nsim,
  ...
) {
  dist <- y
  valid_dists <- eval(formals(energyGOF.test.character)$y)
  distname <- match.arg(tolower(dist), choices = valid_dists)
  dots <- list(...)
  dist <- char_to_dist(distname, ...)
  energyGOFdist(x, dist, nsim)
}


#### Validation

##### Validate distribution function
validate_cdf <- function(y, x, n = 1e5, tol1 = 1e-10, tol2 = 1e-2, ...) {
  # Grid of points over support subset
  xseq <- seq(min(x), max(x), length.out = n)
  first_arg <- names(formals(y))[1]
  args <- c(setNames(list(xseq), first_arg), list(...))
  vals <- do.call(y, args)
  if (any(vals < -tol1 | vals > 1 + tol1)) {
    warning("Distribution function [0, 1] range violation.")
  }
  diffs <- diff(vals)
  if (any(diffs < -tol1)) {
    warning("Distribution function seems to not be monotonic.")
  }
  if (any(diffs > tol2)) {
    warning("Distribution function may not be continuous.")
  }
  d <- data.frame(x = xseq, y = vals)
  ## Seems like a bad method for checking
  if (NROW(d[with(d, x < 0 & y == 0), ]) > 0) {
    warning(
      "You may be testing negative data against a distribution with positive support."
    )
  }
}

##### Validate Parameters
validate_par <- function(dist) {
  if (!dist$par_domain(dist$par)) {
    stop(sprintf(
      "Parameters passed in ... failed domain check:  %s",
      paste0(deparse(body(dist$par_domain)), collapse = "")
    ))
  }
  if (!dist$par_domain(dist$sampler_par)) {
    stop(sprintf(
      "Parameters passed in ... failed domain check:  %s",
      paste0(deparse(body(dist$par_domain)), collapse = "")
    ))
  }
}


##### Validate x
validate_x <- function(x, dist) {
  if (any(is.na(x), is.null(x), is.infinite(x))) {
    stop("Missing data are not supported.")
  }
  if (!dist$support(x, dist$par)) {
    stop(sprintf(
      "Not all elements of x lie in the support of distribution: %s
Support test:  %s",
      dist$name,
      paste0(deparse(body(dist$support)), collapse = "")
    ))
  }
}

##### Validate nsim
validate_nsim <- function(nsim) {
  if (!is.numeric(nsim)) {
    stop("nsim must be numeric.")
  }
  if (!(nsim >= 0)) {
    stop("nsim must be non-negative.")
  }
  floor(nsim)
}

### Switchers
#### Distribution Switcher
char_to_dist <- function(name, ...) {
  ## ... should have the params and pow
  switch(
    name,
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
    "invgaussian" = inverse_gaussian_dist(...),
    "halfnormal" = halfnormal_dist(...),
    "chisq" = chisq_dist(...),
    "chisquared" = chisq_dist(...),
    "f" = f_dist(...),
    "binomial" = binomial_dist(...),
    "bernoulli" = bernoulli_dist(...),
    "geometric" = geometric_dist(...),
    "poisson" = poisson_dist(...),
    ## Generalized GOF Dists
    "cauchy" = cauchy_dist(...),
    ##"stable" = stable_dist(...),
    "pareto" = pareto_dist(...),
    stop("Unsupported distribution: ", name)
  )
}

#### Transform x for some tests
## Separate function for outside bootstrap loop

xform_x <- function(x, dist) {
  UseMethod("xform_x", dist)
}

#' @export
xform_x.CauchyDist <- function(x, dist) {
  # Must transform in Simple case.
  if (!dist$composite) {
    x <- dist$xform(x, dist$par)
  }
  x
}

#' @export
xform_x.StableDist <- function(x, dist) {
  # Must transform in Simple case.
  ## I don't think there will be a composite test, but i can leave the if statement.
  if (!dist$composite) {
    x <- dist$xform(x, dist$par)
  }
  x
}


#' @export
xform_x.ParetoDist <- function(x, dist) {
  if (dist$par$shape > 1 && dist$par$pow != 1) {
    # New ingredients
    x <- dist$xform(x, dist$par)
  }
  x
}

#' @export
xform_x.GOFDist <- function(x, dist) {
  x
}

#### Transform dist for some tests
xform_dist <- function(x, dist) {
  UseMethod("xform_dist", dist)
}

#' @export
xform_dist.PoissonDist <- function(x, dist) {
  # Must transform in Simple case.
  if (dist$composite) {
    dist$sampler_par <- dist$statistic(x)
  }
  dist
}

#' @export
xform_dist.InverseGaussianDist <- function(x, dist) {
  # Must transform in composite case.
  if (dist$composite) {
    dist$sampler_par <- dist$statistic(x)
  }
  dist
}

#' @export
xform_dist.GammaDist <- function(x, dist) {
  # Must transform in Composite case.
  if (dist$composite) {
    dist$sampler_par <- dist$statistic(x)
  }
  dist
}

#' @export
xform_dist.AsymmetricLaplaceDist <- function(x, dist) {
  # Must transform in Composite case.
  if (dist$composite) {
    mle <- dist$statistic(x)
    dist$sampler_par <- list(location = 0, scale = 1, skew = mle$skew)
  }
  dist
}

#' @export
xform_dist.WeibullDist <- function(x, dist) {
  # Must transform in Composite case.
  if (dist$composite) {
    dist$sampler_par <- dist$statistic(x)
  }
  dist
}

#' @export
xform_dist.GOFDist <- function(x, dist) {
  # Done
  dist
}
#### EXXhat

EXXhat <- function(x, dist) {
  UseMethod("EXXhat", dist)
}

#' @export
EXXhat.EuclideanGOFDist <- function(x, dist) {
  n <- length(x)
  xs <- sort(x)
  prefix <- 2 * seq_len(n) - 1 - n
  2 * mean(prefix * xs) / n
}

#' @export
EXXhat.GeneralizedGOFDist <- function(x, dist) {
  pow <- dist$sampler_par$pow
  n <- length(x)
  2 * sum(dist(x)^pow) / n^2
}


#### Compute Energy GOF statistic

Qhat <- function(x, dist, EYY) {
  UseMethod("Qhat", dist)
}

#' @export
Qhat.CompositeGOFDist <- function(x, dist, EYY) {
  mle <- dist$statistic(x)
  x <- dist$xform(x, mle)
  NextMethod(object = dist)
}

#' @export
Qhat.GOFDist <- function(x, dist, EYY) {
  n <- length(x)
  EXY <- dist$EXYhat(x, dist$sampler_par)
  EXX <- EXXhat(x, dist)
  n * (2 * EXY - EYY - EXX)
}

#### Simulate P-values
simulate_pval <- function(x, dist, nsim, Qhat, EYY) {
  if (nsim == 0) {
    return(list(sim_reps = 0, p_value = NA))
  }
  bootobj <- boot::boot(
    x,
    statistic = Qhat,
    R = nsim,
    sim = "parametric",
    ran.gen = dist$sampler,
    mle = dist$sampler_par,
    dist = dist,
    EYY = EYY
  )
  list(
    sim_reps = bootobj$t,
    p_value = mean(bootobj$t > bootobj$t0)
  )
}


#### Output Htest

output_htest <- function(x, dist, nsim, E_stat, sim) {
  cp <- inherits(dist, "CompositeGOFDist")
  gen <- inherits(dist, "GeneralizedGOFDist")
  names(E_stat) <- paste0("E-statistic")
  if (cp) {
    mle <- unlist(dist$statistic(x))
  }
  structure(
    list(
      method = paste0(
        (if (cp) "Composite" else "Simple"),
        " energy goodness-of-fit test",
        (if (cp) {
          paste0(" (conditional on ", deparse(substitute(x)), ")")
        })
      ),
      data.name = deparse(substitute(x)),
      distribution = dist,
      parameter = c("Distribution" = dist$name, if (cp) NULL else dist$par),
      nsim = nsim,
      composite_p = cp,
      statistic = E_stat,
      p.value = sim$p_value,
      sim_reps = sim$sim_reps,
      estimate = if (cp) mle else NULL
    ),
    class = "htest"
  )
}

#### energyGOFdist Generic & Methods

#' @title S3 Interface to Parametric Goodness-of-Fit Tests via Energy
#' @description This is an alternative interface that provides the same
#'   parametric tests as [energyGOF.test()], but allows the user to directly
#'   pass a distribution object like [normal_dist()] (Distribution objects are
#'   specific to the implementation of this R package). The advantage is that
#'   you do not need to pass distribution parameters into a `...` argument as
#'   in `energyGOF.test`. `energyGOF.test` uses this function under the hood,
#'   but it's perfectly suitable for the user to use as well.
#' @param x A numeric vector.
#' @param dist An object of class GOFDist. The distribution to test `x`
#'   against. GOFDist objects are created with the various "`name_dist()`"
#'   functions in this package. See, for example, [normal_dist()] for details
#'   on these class objects.
#' @param nsim A non-negative integer. The number of parametric bootstrap
#'   replicates taken to calculate the *p*-value. If 0, no simulation.
#' @inherit energyGOF.test return author
#' @return  Return an object of class `htest'
#'   representing the result of the energy goodness-of-fit hypothesis test. The
#'   htest object has the elements:
#'
#'  * `method`: Simple or Composite
#'  * `data.name`
#'  * `distribution`: The distribution object created to test
#'  * `parameter`: List of parameters if the test is simple
#'  * `nsim`: Number of bootstrap replicates
#'  * `composite_p`: TRUE/FALSE composite predicate
#'  * `statistic`: The value of the energy statistic (\eqn{Q=nE^*})
#'  * `p.value`
#'  * `sim_reps`: bootstrap replicates of energy statistic
#'  * `estimate`: Any parameter estimates, if the test is composite
#'
#' @aliases egofd
#' @examples
#' ## Simple normal test
#' energyGOFdist(rnorm(10), normal_dist(0, 1), nsim = 10)
#'
#' ## Simple Poisson test
#' egofd(rpois(10,1), poisson_dist(1), nsim = 0) # No p-value
#'
#' ## Composite Normal test
#' egofd(rnorm(10), normal_dist(), nsim = 10)
#'
#'
#' @export energyGOFdist
energyGOFdist <- function(x, dist, nsim) {
  validate_x(x, dist)
  nsim <- validate_nsim(nsim)
  UseMethod("energyGOFdist", dist)
}

#' @rdname energyGOFdist
#' @export egofd
egofd <- energyGOFdist

#' @export
energyGOFdist.GOFDist <- function(x, dist, nsim) {
  ## Setup
  cp <- inherits(dist, "CompositeGOFDist")
  ## Run functions
  x <- xform_x(x, dist)
  dist <- xform_dist(x, dist)
  EYY <- dist$EYY(dist$sampler_par)
  E_stat <- Qhat(x, dist, EYY)
  sim <- simulate_pval(x, dist, nsim, Qhat, EYY)
  output_htest(x, dist, nsim, E_stat, sim)
}

#### Distributions
is_composite <- function(...) {
  nulls <- sapply(list(...), is.null)
  n_null <- sum(nulls)
  if (n_null == 0) {
    FALSE # simple test: all pars supplied
  } else if (n_null == length(nulls)) {
    TRUE # composite test: all pars NULL
  } else {
    stop("Partially composite tests not implemented.")
  }
}

composite_not_allowed <- function(...) {
  nulls <- sapply(list(...), is.null)
  n_null <- sum(nulls)
  if (n_null == 0) {
    FALSE # simple test: all pars supplied
  } else {
    stop("Composite test not implemented for this distribution.")
  }
}


set_composite_class <- function(dist) {
  if (dist$composite_p) {
    class(dist) <- append(class(dist), "CompositeGOFDist", 2L)
  } else {
    class(dist) <- append(class(dist), "SimpleGOFDist", 2L)
  }
  dist
}

#' @export
print.GOFDist <- function(x, ...) {
  cat(" Energy goodness-of-fit test for:\n")
  cat("   *", x$name, "Distribution\n")
  if (!x$composite_p) {
    cat(
      "   * Test Parameters: ",
      paste(names(x$par), unlist(x$par), sep = "=", collapse = ", "),
      "\n"
    )
  }
  cat(
    "   * Sampler Parameters: ",
    paste(
      names(x$sampler_par),
      unlist(x$sampler_par),
      sep = "=",
      collapse = ", "
    ),
    "\n"
  )
  cat(
    "   * Test type:",
    if (x$composite_p) {
      "Composite (parameters unknown)"
    } else {
      "Simple (parameters known)"
    },
    "\n"
  )
  cat("   * S3 Classes:", class(x), "\n")
}

##### Normal

#' @title Create a Normal distribution object for energy testing
#' @author John T. Haman
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against a normal
#'   distribution. If `mean` and `sd` are both NULL, perform a composite test.
#' @param mean NULL, or if specified, same as [rnorm()], but must be length 1.
#' @param sd NULL, or if specified, Same as [rnorm()], but must be length 1
#'
#' @return S3 data object containing the following fields.
#' * `name`: String
#' * `composite_p`: Composite predicate. TRUE if test is composite.
#' * `par`: Distribution parameters, list of the formals.
#' * `sampler_par`: Distribution parameters used for the calculation of energy
#' statistic. These may be different than `par`.
#' * `par_domain`: Function used to ensure `par` and `sampler_par` are valid for
#' this distribution
#' * `support`: Function to check that data `x` can be tested against `y`
#' * `sampler`: Function used for rng by [boot::boot()]
#' * `EYY`: Function to compute \eqn{E|Y-Y'|} (or \eqn{E|Y-Y'|^{pow}}, for the
#' generalized test.)
#' * `EXYhat`: Function to compute \eqn{\frac{1}{n} \sum_i E|x_i - Y|} (or
#' \eqn{\frac{1}{n} \sum_i E|x_i - Y|^{pow}}), where Y is distributed according
#' to `y` and x is the data under test (which is passed in `egof.test` or `egofd`).
#' * `xform`: Function that may be used to transform x. Only available in certain
#' distribution objects.
#' * `statistic`: Function that returns a list of maximum likelihood estimates.
#' Only available in certain distribution objects.
#' * `notes`: Distribution specific messages. Only used in certain distribution
#' objects.
#'
#' *Note*: Some distributions do not have notes, xform, and statistic fields.
#' This is because either a composite test is not implemented, or because a
#' data transformation is not needed.
#'
#' @examples
#' d <- normal_dist(0, 1)
#'
#' # Composite test
#' dc <- normal_dist()
#' egofd(rnorm(10), dc, 0)
#'
#'
#' ### Expected distances:
#'
#' d$EYY(d$par)
#'
#' ## should be about the same as mean(abs(rnorm(1e5) - rnorm(1e5)))
#'
#' x <- 3
#'
#' d$EXYhat(3, d$par)
#'
#' ## should be about the same as mean(abs(x - rnorm(1e5)))
#'
#' @export
normal_dist <- function(mean = NULL, sd = NULL) {
  cp <- is_composite(mean, sd)
  dist <- structure(
    list(
      name = "Normal",
      composite_p = cp,
      par = list(mean = mean, sd = sd),
      sampler_par = if (cp) {
        list(mean = 0, sd = 1)
      } else {
        list(mean = mean, sd = sd)
      },
      par_domain = function(par) {
        all(
          length(mean) == 1 || is.null(mean),
          length(sd) == 1 || is.null(sd),
          par$sd > 0 || is.null(par$sd),
          is.finite(par$mean) || is.null(par$mean)
        )
      },
      support = function(x, par) all(is.finite(x)),
      sampler = function(n, par) rnorm(n, par$mean, par$sd),
      EYY = function(par) 2 * par$sd / sqrt(pi),
      EXYhat = function(x, par) {
        mean(
          2 *
            (x - par$mean) *
            pnorm(x, par$mean, par$sd) +
            2 * par$sd^2 * dnorm(x, par$mean, par$sd) -
            (x - par$mean)
        )
      },
      xform = function(x, par) (x - par$mean) / par$sd,
      statistic = function(x) {
        list(mean = mean(x), sd = sd(x))
      }
    ),
    class = c("NormalDist", "EuclideanGOFDist", "GOFDist")
  )
  validate_par(dist)
  set_composite_class(dist)
}


##### Uniform

#' @title Create a Uniform distribution object for energy testing
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against a uniform
#'   distribution. Only simple tests are implemented.
#' @inherit normal_dist return author
#' @param min Same as in [runif()], but must be length 1
#' @param max Same as in [runif()], but must be length 1
#'
#' @examples
#'
#' d <- uniform_dist(0, 1)
#'
#' egofd(runif(10), d, 0)
#'
#'
#' @export
uniform_dist <- function(min = 0, max = 1) {
  dist <- structure(
    list(
      name = "Uniform",
      composite_p = composite_not_allowed(min, max),
      par = list(min = min, max = max),
      sampler_par = list(min = min, max = max),
      par_domain = function(par) {
        all(
          length(par$min) == 1,
          length(par$max) == 1,
          par$max - par$min > 0
        )
      },
      support = function(x, par) is.numeric(x),
      sampler = function(n, par) runif(n, par$min, par$max),
      EYY = function(par) (par$max - par$min) / 3,
      EXYhat = function(x, par) {
        mean(
          (x - par$min)^2 / (par$max - par$min) - x + (par$max - par$min) / 2
        )
      },
      xform = function(x, par) (x - par$min) / (par$max - par$min),
      statistic = function(x) {
        list(min = min(x), max = max(x))
      }
    ),
    class = c("UniformDist", "EuclideanGOFDist", "SimpleGOFDist", "GOFDist")
  )
  validate_par(dist)
  dist
}
##### Exponential

#' @title Create an Exponential distribution object for energy testing
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against an
#'   exponential distribution. If rate is NULL, a composite test is performed.
#' @inherit normal_dist return author
#' @param rate NULL, or a positive rate parameter as in [rexp()], but must be
#'   length 1.
#' @aliases exp_dist
#' @examples
#'
#' d <- exponential_dist(1)
#' egofd(rexp(10, 1), d, 0)
#'
#' @export
exponential_dist <- function(rate = NULL) {
  cp <- is.null(rate)
  dist <- structure(
    list(
      name = "Exponential",
      composite_p = cp,
      par = list(rate = rate),
      sampler_par = if (cp) list(rate = 1) else (list(rate = rate)),
      par_domain = function(par) {
        (par$rate > 0 && length(par$rate) == 1) || is.null(par$rate)
      },
      support = function(x, par) all(x > 0),
      sampler = function(n, par) rexp(n, par$rate),
      EYY = function(par) 1 / par$rate,
      EXYhat = function(x, par) {
        mean(x + 1 / par$rate * (1 - 2 * pexp(x, par$rate)))
      },
      xform = function(x, par) x / par$rate,
      statistic = function(x) list(rate = mean(x))
    ),
    class = c("ExponentialDist", "EuclideanGOFDist", "GOFDist")
  )
  validate_par(dist)
  set_composite_class(dist)
}

#' @export
exp_dist <- exponential_dist

##### Poisson

#' @title Create a Poisson distribution object for energy testing
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against a Poisson
#'   distribution. If lambda is NULL, a composite test is performed.
#' @inherit normal_dist return author
#' @param lambda NULL, or if specified, same as the lambda in [rpois()], but
#'   must be length 1.
#'
#' @examples
#' d <- poisson_dist(1)
#'
#' egofd(rpois(10, 1), d, 0)
#'
#' @export
poisson_dist <- function(lambda = NULL) {
  cp <- is.null(lambda)
  dist <- structure(
    list(
      name = "Poisson",
      composite_p = cp,
      par = list(lambda = lambda),
      sampler_par = if (cp) {
        list(lambda = NULL)
      } else {
        list(lambda = lambda)
      },
      par_domain = function(par) {
        (par$lambda > 0 && length(par$lambda) == 1) || is.null(par$lambda)
      },
      support = function(x, par) {
        all(x >= 0) && all(is.integer(x))
      },
      sampler = function(n, par) {
        rpois(n, par$lambda)
      },
      EYY = function(par) {
        2 *
          par$lambda *
          exp(-2 * par$lambda) *
          (besselI(2 * par$lambda, 0) +
            besselI(2 * par$lambda, 1))
      },
      EXYhat = function(x, par) {
        n <- length(x)
        mean(
          2 *
            x *
            ppois(x, par$lambda) -
            2 * par$lambda * ppois(x - 1, par$lambda) +
            par$lambda -
            x
        )
      },
      xform = function(x, par) x,
      statistic = function(x) list(lambda = mean(x))
    ),
    class = c("PoissonDist", "EuclideanGOFDist", "GOFDist")
  )
  validate_par(dist)
  set_composite_class(dist)
}


##### Skew-Normal?

##### Bernoulli

#' @title Create a Bernoulli distribution object for energy testing
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against a Bernoulli
#'   distribution. Only simple tests are implemented.
#' @inherit normal_dist  return author
#' @param prob Same as [rbinom()], but must be length 1.
#' @examples
#'
#' d <- bernoulli_dist(.5)
#' egofd(rbinom(10, 1, .5), d, 0)
#'
#' @export
bernoulli_dist <- function(prob = 0.5) {
  dist <- structure(
    list(
      name = "Bernoulli",
      composite_p = composite_not_allowed(prob),
      par = list(prob = prob),
      sampler_par = list(prob = prob),
      par_domain = function(par) {
        all(
          par$prob > 0 && par$prob < 1,
          length(par$prob) == 1
        )
      },
      support = function(x, par) all(x %in% c(0L, 1L)),
      sampler = function(n, par) {
        rbinom(n, size = 1, prob = par$prob)
      },
      EYY = function(par) {
        2 * par$prob * (1 - par$prob)
      },
      EXYhat = function(x, par) {
        h <- sum(x)
        n <- length(x)
        (h * (1 - par$prob) + (n - h) * par$prob) / n
      },
      statistic = function(x) list(prob = mean(x))
    ),
    class = c("BernoulliDist", "EuclideanGOFDist", "SimpleGOFDist", "GOFDist")
  )
  validate_par(dist)
  dist
}

##### Binomial

#' @title Create a Binomial distribution object for energy testing
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against a Binomial
#'   distribution. Only a simple GOF test is supported.
#' @inherit normal_dist return author
#' @param prob Same as [stats::rbinom()], but must be length 1.
#' @param size Same as [stats::rbinom()], but must be length 1.
#' @examples
#'
#' d <- binomial_dist(1, 0.5)
#' egofd(rbinom(10, 1, .5), d, 0)
#'
#' @export
binomial_dist <- function(size = 1, prob = 0.5) {
  dist <- structure(
    list(
      name = "Binomial",
      composite_p = composite_not_allowed(prob),
      par = list(size = size, prob = prob),
      sampler_par = list(size = size, prob = prob),
      par_domain = function(par) {
        all(
          length(par$size) == 1,
          length(par$prob) == 1,
          par$prob > 0 && par$prob < 1,
          par$size >= 1,
          par$size == floor(par$size)
        )
      },
      support = function(x, par) all(x %in% 0:par$size),
      sampler = function(n, par) {
        rbinom(n, size = par$size, prob = par$prob)
      },
      EYY = function(par) {
        i <- 0:par$size
        probs <- dbinom(i, size = par$size, prob = par$prob)
        diffmat <- abs(outer(i, i, "-"))
        sum(diffmat * outer(probs, probs))
      },
      EXYhat = function(x, par) {
        i <- 0:par$size
        probs <- dbinom(i, size = par$size, prob = par$prob)
        mean(sapply(x, function(xi) sum(abs(xi - i) * probs)))
      },
      statistic = function(x) list(prob = mean(x))
    ),
    class = c("BinomialDist", "EuclideanGOFDist", "SimpleGOFDist", "GOFDist")
  )
  validate_par(dist)
  dist
}


##### Beta

#' @title Create a beta distribution object for energy testing
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against a beta
#'   distribution. If shape1 and shape2 are NULL, a composite test is
#'   performed, otherwise a simple test is performed.
#' @inherit normal_dist return author
#' @param shape1 Same as [rbeta()], but must be length 1.
#' @param shape2 Same as [rbeta()], but must be length 1.
#'
#' @examples
#' d <- beta_dist(5, 5)
#' egofd(rbeta(10, 5, 5), d, 0)
#'
#' @export
beta_dist <- function(shape1 = NULL, shape2 = NULL) {
  cp <- is_composite(shape1, shape2)
  dist <- structure(
    list(
      name = "Beta",
      composite_p = cp,
      par = list(shape1 = shape1, shape2 = shape2),
      sampler_par = if (!cp) {
        list(shape1 = shape1, shape2 = shape2)
      } else {
        list(shape1 = 1, shape2 = 1)
      },
      par_domain = function(par) {
        all(
          (par$shape1 > 0 && length(par$shape1) == 1) || is.null(par$shape1),
          (par$shape2 > 0 && length(par$shape2) == 1) || is.null(par$shape2)
        )
      },
      sampler = function(n, par) {
        rbeta(n, shape1 = par$shape1, shape2 = par$shape2)
      },
      support = function(x, par) all(x < 1) && all(x > 0),
      EYY = function(par) {
        integrand <- function(x, par) {
          a <- par$shape1
          b <- par$shape2
          ExY <- 2 *
            x *
            pbeta(x, a, b) -
            x +
            a / (a + b) -
            2 * pbeta(x, a + 1, b) * beta(a + 1, b) / beta(a, b)
          ExY * dbeta(x, a, b)
        }
        integrate(integrand, 0, 1, par)$value
      },
      # incorrect in book
      EXYhat = function(x, par) {
        a <- par$shape1
        b <- par$shape2
        mean(
          2 *
            x *
            pbeta(x, a, b) -
            x +
            a / (a + b) -
            2 * beta(a + 1, b) / beta(a, b) * pbeta(x, a + 1, b)
        )
      },
      xform = function(x, par) {
        ## Probability integral transform
        pbeta(x, par$shape1, par$shape2)
      },
      statistic = function(x) {
        as.list(fitdistrplus::mledist(x, "beta")$estimate)
      }
    ),
    class = c("BetaDist", "EuclideanGOFDist", "GOFDist")
  )
  validate_par(dist)
  set_composite_class(dist)
}

##### Dirchlet?

##### Geometric

#' @title Create a geometric distribution object for energy testing
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against a geometric
#'   distribution. Only a simple test is supported.
#' @inherit normal_dist  return author
#' @param prob Same as [rgeom()], but must be length 1.
#' @examples
#'
#' d <- geometric_dist(.5)
#' egofd(rgeom(10, .5), d, 0)
#'
#' @export
geometric_dist <- function(prob = 0.5) {
  dist <- structure(
    list(
      name = "Geometric",
      composite_p = composite_not_allowed(prob),
      par = list(prob = prob),
      sampler_par = list(prob = prob),
      par_domain = function(par) {
        all(
          par$prob > 0 && par$prob < 1,
          length(par$prob) == 1
        )
      },
      support = function(x, par) all(x == floor(x)) && all(x >= 0),
      sampler = function(n, par) rgeom(n, par$prob),
      EYY = function(par) {
        q <- 1 - par$prob
        (2 * q) / (1 - q^2)
      },
      EXYhat = function(x, par) {
        mean(x + 1 + (1 - 2 * pgeom(x, par$prob)) / par$prob)
      }
    ),
    class = c("GeometricDist", "EuclideanGOFDist", "SimpleGOFDist", "GOFDist")
  )
  validate_par(dist)
  dist
}


##### Negative Binomial?

##### Half-Normal
#' @title Create a half-normal distribution object for energy testing
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against a
#'   half-normal distribution. If scale is NULL, a composite test is performed.
#' @inherit normal_dist return author
#' @param scale NULL, or a positive scale parameter, like sd in [rnorm()]. Must
#'   be length 1.
#' @description This is exactly the distribution of \eqn{|X|}, where \eqn{X ~
#'   N(0,\theta = scale)}
#'
#' @examples
#'
#' d <- halfnormal_dist(4)
#' egofd(abs(rnorm(10, 4)), d, 0)
#'
#' @export
halfnormal_dist <- function(scale = NULL) {
  cp <- is.null(scale)
  dist <- structure(
    list(
      name = "Half-Normal",
      composite_p = cp,
      par = list(scale = scale),
      sampler_par = if (!cp) {
        list(scale = scale)
      } else {
        list(scale = 1)
      },
      par_domain = function(par) {
        (par$scale > 0 && length(par$scale) == 1) || is.null(par$scale)
      },
      support = function(x, par) all(x > 0),
      sampler = function(n, par) {
        abs(rnorm(n, 0, sd = par$scale))
      },
      EXYhat = function(x, par) {
        scale <- par$scale
        mean(
          2 *
            x *
            (2 * pnorm(x, 0, scale) - 1) -
            x +
            scale * sqrt(2 / pi) -
            2 * sqrt(2 / pi) * scale * (1 - exp(-x^2 / (2 * scale^2)))
        )
      },
      EYY = function(par) {
        par$scale * 2 * (2 - sqrt(2)) / sqrt(pi)
      },
      xform = function(x, par) x / par$scale,
      statistic = function(x) {
        n <- length(x)
        list(scale = sqrt(sum(x^2) / n))
      }
    ),
    class = c("HalfNormalDist", "EuclideanGOFDist", "GOFDist")
  )
  validate_par(dist)
  set_composite_class(dist)
}

##### Laplace

#' @title Create a Laplace distribution object for energy testing
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against a Laplace
#'   distribution. If location and scale are both NULL, a composite test is
#'   performed.
#' @inherit normal_dist  return author
#' @param location NULL, or the median of the distribution
#' @param scale NULL or a positive scale parameter
#' @description This is exactly the distribution corresponding to the PDF
#'
#' \deqn{f(x|\mu, b) = \frac{1}{2b} \exp \left(-\frac{|x - \mu|}{b} \right), }
#'
#' where `location` = \eqn{\mu} and `scale` = \eqn{b}.
#'
#' @examples
#' d <- laplace_dist(1, 1)
#'
#' x <- d$sampler(10, d$par)
#'
#' egofd(x, d, 0)
#'
#' @export
laplace_dist <- function(location = NULL, scale = NULL) {
  cp <- is_composite(location, scale)
  dist <- structure(
    list(
      name = "Laplace",
      composite_p = cp,
      par = list(location = location, scale = scale),
      sampler_par = if (cp) {
        list(location = 0, scale = 1)
      } else {
        list(location = location, scale = scale)
      },
      par_domain = function(par) {
        all(
          (par$scale > 0 && length(par$scale) == 1) || is.null(par$scale),
          (length(par$location) == 1) || is.null(par$location)
        )
      },
      support = function(x, par) all(is.finite(x)),
      sampler = function(n, par) {
        par$location + sign(runif(n) - 0.5) * rexp(n, 1 / par$scale)
      },
      EXYhat = function(x, par) {
        mean(
          abs(x - par$location) +
            par$scale * exp(-abs(x - par$location) / par$scale)
        )
      },
      EYY = function(par) {
        par$scale * (3 / 2)
      },
      statistic = function(x) {
        list(location = median(x), scale = mean(abs(x - median(x))))
      },
      xform = function(x, par) {
        (x - par$location) / par$scale
      }
    ),
    class = c("LaplaceDist", "EuclideanGOFDist", "GOFDist")
  )
  validate_par(dist)
  set_composite_class(dist)
}

##### Log-Normal

#' @title Create a lognormal distribution object for energy testing
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against a lognormal
#'   distribution. If `meanlog` and `sdlog` are both `NULL`, a composite test is
#'   performed.
#' @inherit normal_dist return author
#' @aliases lnorm_dist
#' @param meanlog NULL or as in [rlnorm()], must be length 1.
#' @param sdlog NULL or as in [rlnorm()], must be length 1.
#'
#' d <- lognormal_dist(0, 1)
#' x <- d$sampler(10, d$par)
#'
#' egofd(x, d, 0)
#'
#' @export
lognormal_dist <- function(meanlog = NULL, sdlog = NULL) {
  cp <- is_composite(meanlog, sdlog)
  dist <- structure(
    list(
      name = "Lognormal",
      composite_p = cp,
      par = list(meanlog = meanlog, sdlog = sdlog),
      sampler_par = if (cp) {
        list(meanlog = 0, sdlog = 1)
      } else {
        list(meanlog = meanlog, sdlog = sdlog)
      },
      support = function(x, par) {
        all(x > 0) && all(is.finite(x))
      },
      par_domain = function(par) {
        all(
          par$sdlog > 0 || is.null(par$sdlog),
          is.finite(par$meanlog) || is.null(par$meanlog)
        )
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
        integrate(integrand, lower = 0, upper = Inf, par = par)$value # TODO
        # should be Inf
      },
      xform = function(x, par) exp((log(x) - par$meanlog) / par$sdlog), # ~ lognormal(0, 1)
      statistic = function(x) {
        as.list(fitdistrplus::mledist(x, "lnorm")$estimate)
      }
    ),
    class = c("LogNormalDist", "EuclideanGOFDist", "GOFDist")
  )
  validate_par(dist)
  set_composite_class(dist)
}


##### Asymmetric Laplace
## TODO fix broken estimator
#' @title Create an asymmetric Laplace distribution object for energy testing
#' @inherit normal_dist return author
#' @aliases alaplace_dist
#' @param location NULL, or a location parameter
#' @param scale NULL, or a positive scale parameter
#' @param skew NULL, or a positive skewness parameter. Skew = 1 corresponds to
#'   a symmetric Laplace distribution (though note the difference between the
#'   PDF in this description and the one in [laplace_dist()]).
#' @aliases alaplace_dist
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against an asymmetric
#'   Laplace distribution. If all three parameters are NULL, perform a
#'   composite test. This is exactly the distribution corresponding to the PDF
#'
#' \deqn{
#'   f(x | \theta, \sigma, \kappa) =
#'   \frac{\sqrt{2}\kappa}{\sigma(1 + \kappa^2)}
#'   \begin{cases}
#'     \exp\Big( -\frac{\sqrt{2} \kappa |x - \theta|}{\sigma} \Big),
#'     & x \ge \theta, \\[6pt]
#'     \exp\Big( -\frac{\sqrt{2} |x - \theta|}{\kappa \sigma} \Big),
#'     & x < \theta.
#'   \end{cases}
#', }
#'
#' where \eqn{\theta} = `location`, \eqn{\sigma} = `scale`, and \eqn{\kappa} = `skew`.
#'
#' @examples
#'
#' d <- asymmetric_laplace_dist(0, 1, .5)
#' x <- d$sampler(10, d$par)
#'
#' egofd(x, d, 0)
#'
#'
#' @export
asymmetric_laplace_dist <- function(
  location = NULL,
  scale = NULL,
  skew = NULL
) {
  dist <- structure(
    list(
      name = "Asymmetric Laplace",
      composite_p = is_composite(location, scale, skew),
      par = list(location = location, scale = scale, skew = skew),
      sampler_par = list(location = location, scale = scale, skew = skew),
      par_domain = function(par) {
        all(
          (par$scale > 0 && length(par$scale) == 1) || is.null(par$scale),
          (par$skew > 0 && length(par$skew) == 1) || is.null(par$skew),
          length(par$location) == 1 || is.null(par$location)
        )
      },
      support = function(x, par) all(is.finite(x)),
      sampler = function(n, par) {
        loc <- par$location
        scale <- par$scale
        k <- par$skew
        loc + scale * log(runif(n)^k / runif(n)^(1 / k)) / sqrt(2)
      },
      EXYhat = function(x, par) {
        loc <- par$location
        scale <- par$scale
        k <- par$skew
        mu <- (1 / k - k) / sqrt(2)
        lam <- sqrt(2) * k / scale
        beta <- sqrt(2) / (k * scale)
        pk <- 1 / (1 + k^2)
        qk <- 1 - pk
        mean(
          ifelse(
            x >= loc,
            x -
              loc -
              mu +
              (2 * pk / lam) *
                exp(-lam * abs(x - loc)),
            -x +
              loc +
              mu +
              (2 * qk / beta) *
                exp(-beta * abs(x - loc))
          )
        )
      },
      EYY = function(par) {
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
      xform = function(x, par) {
        (x - par$location) / par$scale # ~ AL(0, 1, kappa)
      },
      statistic = function(x) {
        AML_EM(x)
      }
    ),
    class = c("AsymmetricLaplaceDist", "EuclideanGOFDist", "GOFDist")
  )
  validate_par(dist)
  set_composite_class(dist)
}

#' @rdname asymmetric_laplace_dist
alaplace_dist <- asymmetric_laplace_dist

##### F Distribution (Fisher Distribution)
#' @title Create an F distribution object for energy testing
#' @inherit normal_dist return author
#' @param df1 Positive.
#' @param df2 Must be greater than 2.
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against a F
#'   distribution. Only simple tests are supported.
#' @examples
#'
#' d <- f_dist(3, 3)
#' egofd(rf(10, 3, 3), d, 0)
#'
#' @export
f_dist <- function(df1 = 3, df2 = 3) {
  cp <- composite_not_allowed(df1, df2)
  dist <- structure(
    list(
      name = "F",
      composite_p = cp,
      par = list(df1 = df1, df2 = df2),
      sampler_par = list(df1 = df1, df2 = df2),
      support = function(x, par) {
        all(x > 0)
      },
      par_domain = function(par) {
        all(
          (par$df1 > 0 && length(par$df1) == 1),
          (par$df2 > 2 && length(par$df2) == 1)
        )
      },
      sampler = function(n, par) {
        rf(n, df1 = par$df1, df2 = par$df2)
      },
      EXYhat = function(x, par) {
        df1 <- par$df1
        df2 <- par$df2
        EX <- df2 / (df2 - 2)
        u <- df1 * x / (df1 * x + df2)
        mean(
          x *
            (2 * pf(x, df1, df2) - 1) +
            EX * (1 - 2 * pbeta(u, df1 / 2 + 1, df2 / 2 - 1))
        )
      },
      ## Wrong
      EYY = function(par) {
        df1 <- par$df1
        df2 <- par$df2
        EY <- df2 / (df2 - 2)
        integrand <- function(t, df1, df2) {
          a <- df1 / 2
          b <- df2 / 2
          (1 - pbeta(t, a, b))^2 / (1 - t)^2
        }
        I <- integrate(integrand, 0, 1, df1 = df1, df2 = df2)$value
        2 * EY - 2 * df2 / df1 * I
      }
    ),
    class = c("FDist", "EuclideanGOFDist", "GOFDist")
  )
  validate_par(dist)
  dist
}


##### Weibull
#' @title Create a Weibull distribution object for energy testing
#' @inherit normal_dist return author
#' @param shape NULL, or if specified, same as the shape parameter in [stats::rweibull()]
#' @param scale NULL, or if specified, same as the scale parameter in
#' [stats::rweibull()]
#' #' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against a Weibull
#'   distribution. If `shape` and `scale` are both NULL, perform a composite test.
#' @examples
#'
#' d <- weibull_dist(3, 3)
#' egofd(rweibull(10, 3, 3), d, 0)
#'
#' @export
weibull_dist <- function(shape = NULL, scale = NULL) {
  cp <- is_composite(shape, scale)
  dist <- structure(
    list(
      name = "Weibull",
      composite_p = cp,
      par = list(shape = shape, scale = scale),
      sampler_par = list(shape = shape, scale = scale),
      support = function(x, par) {
        all(x > 0)
      },
      par_domain = function(par) {
        all(
          (par$shape > 0 && length(par$shape) == 1) || is.null(par$shape),
          (par$scale > 0 && length(par$scale) == 1) || is.null(par$shape)
        )
      },
      sampler = function(n, par) {
        rweibull(n, shape = par$shape, scale = par$scale)
      },
      EXYhat = function(x, par) {
        if (!cp) {
          z = (x / par$scale)^par$shape
          mean(
            2 *
              x *
              pweibull(x, par$shape, par$scale) -
              x +
              par$scale *
                gamma(1 + 1 / par$shape) *
                (1 - 2 * pgamma(z, 1 + 1 / par$shape, 1))
          )
        } else {
          mean(x + 1 * (1 - 2 * pexp(x, 1)))
        }
      },
      EYY = function(par) {
        if (!cp) {
          # scale = lambda
          # shape = k
          (2 * par$scale / par$shape) *
            gamma(1 / par$shape) *
            (1 - 2^(-1 / par$shape))
        } else {
          1
        }
      },
      statistic = function(x) {
        as.list(fitdistrplus::mledist(x, "weibull")$estimate)
      },
      xform = function(x, par) {
        (x / par$scale)^par$shape # ~ exp(1)
      }
    ),
    class = c("WeibullDist", "EuclideanGOFDist", "GOFDist")
  )
  validate_par(dist)
  set_composite_class(dist)
}

##### Gamma

#' @title Create a gamma distribution object for energy testing
#' @inherit normal_dist return author
#' @param shape Same shape parameter in [stats::rgamma()] (must be length 1)
#' @param rate Same rate parameter in [stats::rgamma()] (must be length 1)
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against a Gamma
#'   distribution. Only simple tests are supported.
#' @examples
#' d <- gamma_dist(4, 4)
#' egofd(rgamma(10, 4, 4), d, 0)
#'
#' @export
gamma_dist <- function(shape = 1, rate = 1) {
  cp <- composite_not_allowed(shape, rate)
  dist <- structure(
    list(
      name = "Gamma",
      composite_p = cp, # will figure this out later
      par = list(shape = shape, rate = rate),
      sampler_par = list(shape = shape, rate = rate),
      support = function(x, par) {
        all(x > 0)
      },
      par_domain = function(par) {
        all(
          (par$shape > 0 && length(par$shape) == 1) || is.null(par$shape),
          (par$rate > 0 && length(par$rate) == 1) || is.null(par$rate)
        )
      },
      sampler = function(n, par) {
        rgamma(n, shape = par$shape, rate = par$rate)
      },
      EXYhat = function(x, par) {
        if (!cp) {
          a <- par$shape
          b <- par$rate
          mean(
            (a / b) -
              (2 * a / b) * pgamma(x, a + 1, b) +
              x * (2 * pgamma(x, a, b) - 1)
          )
        } else {
          df <- 2 * par$shape # will be the statistic from xform_dist.GammaDist
          mean(2 * x * pchisq(x, df) - x + df - 2 * df * pchisq(x, df + 2))
        }
      },
      EYY = function(par) {
        if (!cp) {
          a <- par$shape
          b <- par$rate
          2 * gamma(a + 1 / 2) / (b * gamma(a) * sqrt(pi))
        } else {
          df <- 2 * par$shape
          4 * gamma((df + 1) / 2) / gamma(df / 2) / sqrt(pi)
        }
      },
      statistic = function(x) {
        as.list(fitdistrplus::mledist(x, "gamma")$estimate)
      },
      xform = function(x, par) {
        x / 2 * par$rate # ~ gamma (a/2, 1/2) ~ chisq (a)
      }
    ),
    class = c("GammaDist", "EuclideanGOFDist", "GOFDist")
  )
  validate_par(dist)
  set_composite_class(dist)
}

##### Inverse Gamma?

##### Chi-Square

#' @title Create a Chi-squared distribution object for energy testing
#' @inherit normal_dist return author
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against a Chi-squared
#'   distribution. Only simple tests are supported.
#' @param df Same as in [stats::rchisq()].
#' @examples
#' d <- chisq_dist(4)
#' egofd(rchisq(10, 4), d, 0)
#'
#' @export
chisq_dist <- function(df = 2) {
  dist <- structure(
    list(
      name = "Chi-Squared",
      composite_p = composite_not_allowed(df),
      par = list(df = df),
      sampler_par = list(df = df),
      support = function(x, par) {
        all(x > 0) && all(is.finite(x))
      },
      par_domain = function(par) {
        par$df > 0 && length(par$df) == 1
      },
      sampler = function(n, par) {
        rchisq(n, df = par$df, ncp = 0)
      },
      EXYhat = function(x, par) {
        v <- par$df
        mean(2 * x * pchisq(x, v) - x + v - 2 * v * pchisq(x, v + 2))
      },
      EYY = function(par) {
        v <- par$df
        4 * gamma((v + 1) / 2) / gamma(v / 2) / sqrt(pi)
      }
    ),
    class = c("ChiSquaredDist", "EuclideanGOFDist", "SimpleGOFDist", "GOFDist")
  )
  validate_par(dist)
  dist
}


##### Inverse Gaussian

#' @title Create an inverse Gaussian distribution object for energy testing
#' @inherit normal_dist return author
#' @aliases invgauss_dist
#' @param mean NULL or a positive mean parameter
#' @param shape NULL or a positive shape parameter
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against an inverse
#'   Gaussian distribution. If `mean` and `shape` are both NULL, perform a
#'   composite test. This is exactly the distribution corresponding to the PDF
#'
#' \deqn{
#'   f(x | \mu, \lambda) =
#'   \left( \frac{\lambda}{2 \pi x^3} \right)^{1/2}
#'   \exp \left( -\frac{\lambda (x - \mu)^2}{2 \mu^2 x} \right),
#'   \qquad x > 0,
#' }
#'
#' where `mean` is \eqn{\mu} and `shape` is \eqn{\lambda}.
#'
#'@details
#' This distribution requires an intense amount of numerical integration for
#' the simple (known parameters) case, and the implementation seems to be fine
#' for samples up to 1000. For the composite case, the data are transformed to
#' a Chi-squared distribution (conditional on the parameter estimates), and the
#' performance is much better, as there is no numerical integration in this
#' case.
#'
#' @examples
#' d <- inverse_gaussian_dist(4, 4)
#' x <- d$sampler(10, d$par)
#'
#' egofd(x, d, 0)
#'
#' @export
inverse_gaussian_dist <- function(mean = NULL, shape = NULL) {
  cp <- is_composite(mean, shape)
  dist <- structure(
    list(
      name = "Inverse Gaussian",
      composite_p = cp,
      par = list(mean = mean, shape = shape),
      support = function(x, par) {
        all(x > 0)
      },
      par_domain = function(par) {
        all(
          (par$mean > 0 && length(par$mean) == 1) || is.null(par$mean),
          (par$shape > 0 && length(par$shape) == 1) || is.null(par$shape)
        )
      },
      sampler = function(n, par) {
        m <- par$mean
        lam <- par$shape
        statmod::rinvgauss(n, m, lam)
      },
      sampler_par = {
        list(mean = mean, shape = shape)
      },
      EXYhat = function(x, par) {
        if (!cp) {
          m <- par$mean
          lam <- par$shape
          Msave <- numeric(length(x))
          M_closed <- function(x, mu, lambda) {
            a <- lambda / (2 * mu^2)
            b <- lambda / 2
            z <- 2 * sqrt(a * b)
            pref <- sqrt(lambda / (2 * pi)) * exp(lambda / mu)
            integrand <- function(y) y^(-0.5) * exp(-a * y - b / y)
            val <- integrate(
              integrand,
              lower = 0,
              upper = x,
              rel.tol = 1e-5
            )$value
            Mx <- pref * val
            Mx
          }
          Msave <- sapply(x, M_closed, mu = m, lambda = lam)
          mean(2 * x * statmod::pinvgauss(x, m, lam) - x + m - 2 * Msave)
        } else {
          d <- chisq_dist(1)
          d$EXYhat(x, list(df = 1))
        }
      },
      EYY = function(par) {
        if (!cp) {
          m <- par$mean
          lam <- par$shape
          integrand <- function(t, m, lam) {
            phi <- sqrt(lam / m)
            erf <- function(w) 2 * pnorm(w * sqrt(2)) - 1
            8 * exp(-t^2) * t * erf(t) / sqrt(pi) / sqrt(t^2 + 2 * phi^2)
          }
          m * integrate(integrand, 0, Inf, m = m, lam = lam)$value
        } else {
          d <- chisq_dist(1)
          d$EYY(list(df = 1))
        }
      },
      statistic = function(x) {
        list(mean = mean(x), shape = {
          n <- length(x)
          x.tilde <- n / sum(1 / x)
          x.bar <- mean(x)
          ilam.hat <- 1 / x.tilde - 1 / x.bar
          1 / ilam.hat
        })
      },
      xform = function(x, par) {
        ## The half-normal transformation seems to not be sensitive?
        m <- par$mean
        lam <- par$shape
        ## abs(sqrt(lam / x) * (x - m) / m)
        lam * (x - m)^2 / m^2 / x # ~ chisq(1)
      }
    ),
    class = c("InverseGaussianDist", "EuclideanGOFDist", "GOFDist")
  )
  validate_par(dist)
  set_composite_class(dist)
}

#' @rdname inverse_gaussian_dist
invgauss_dist <- inverse_gaussian_dist

##### Inverse Gamma?

##### Gumbel?

##### Generalized Lambda?

#### Generalized Goodness-of-fit Tests

##### Pareto

#' @title Create a Pareto (type I) distribution object for energy testing
#' @inherit normal_dist return author
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the energy goodness-of-fit test against a Pareto
#'   distribution. If `scale` and `shape` are both NULL, perform a composite
#'   test.
#' @param scale NULL or a positive scale parameter
#' @param shape NULL or a positive shape parameter. If shape > 1, shape is used
#'   to transform x
#' @param pow Optional exponent of the energy test. Pow must be less than
#'   shape. If shape > 1 and pow != 1, pow will be scaled down.
#' @details If shape > 1, the energy test is more difficult, so data are
#'   transformed to data^shape ~ Pareto(scale^shape, 1).
#' @examples
#' d <- pareto_dist(1, .5)
#' x <- d$sampler(10, d$par)
#' egofd(x, d, 0)
#'
#' @export
pareto_dist <- function(scale = NULL, shape = NULL, pow = NULL) {
  if (is.null(pow)) {
    if (is.null(shape)) {
      pow <- 1 / 4
    } else if (shape > 1) {
      pow <- 1
    } else if (shape == 1) {
      pow <- 0.5
    } else if (shape < 1) {
      pow <- shape / 2
    }
  }
  dist <- structure(
    list(
      name = "Pareto (Type I)",
      composite_p = is_composite(scale, shape),
      par = list(scale = scale, shape = shape, pow = pow),
      sampler_par = list(scale = 1, shape = 1, pow = 0.25),
      support = function(x, par) {
        all(x > par$scale)
      },
      par_domain = function(par) {
        all(
        (par$scale > 0 && length(par$scale) == 1) || is.null(par$scale),
        (par$shape > 0 && length(par$shape) == 1) || is.null(par$shape),
        par$pow > 0,
        any(
          is.null(par$shape),
          par$shape == 1 && par$pow < 1,
          par$shape < 1 && par$pow < par$shape,
          par$shape > 1 && par$pow == 1
        )
        )
      },
      sampler = function(n, par) {
        u <- runif(n)
        par$scale / u^(1 / par$shape)
      },
      EXYhat = function(x, par) {
        shape <- par$shape
        scale <- par$scale
        pow <- par$pow
        x0 <- (x - scale) / x
        if (shape == 1) {
          A <- (x - scale)^pow
          B <- scale * pow * x^(pow - 1)
          C <- x0^pow /
            pow +
            x0^(pow + 1) / (pow + 1) * gsl::hyperg_2F1(1, pow + 1, pow + 2, x0)
          D <- scale * x^(pow - 1) * beta(pow + 1, 1 - pow)
          mean(A - B * C + D)
        } else if (shape > 1 && pow == 1) {
          mean(
            x +
              (2 * scale^shape * x^(1 - shape) - shape * scale) /
              (shape - 1)
          )
        } else {
          ## Shape < 1 and pow < 1/2
          mean(
          (x - scale)^pow -
            scale^shape *
            (pow *
               beta(pow, 1 - shape) *
               pbeta(x0, pow, 1 - shape) -
               shape * beta(shape - pow, pow + 1)) /
            x^(shape - pow)
          )
        }
      },
      EYY = function(par) {
        shape <- par$shape
        scale <- par$scale
        pow <- par$pow
        if (shape == 1) {
          2 * scale^pow / (2 - pow) * beta(1 - pow, pow + 1)
        } else if (shape > 1 && pow == 1) {
          2 * shape * scale / (shape - 1) / (2 * shape - 1)
        } else {
          ## Shape < 1 and pow < 1/2
          #2 * shape^2 * scale^pow * beta(shape - pow, pow + 1) / (2 * shape -
          #pow)
          ## I thought this was unstable, so i wrote on log scale, but I no
          ## longer believe it to be unstable.
          L <- log(2) +
            2 * log(shape) +
            pow * log(scale) +
            lbeta(shape - pow, pow + 1) -
            log(2 * shape - pow)
          exp(L)
        }
      },
      statistic = function(x) {
        list(scale = min(x), shape = {
          n <- length(x)
          n / (sum(log(x / min(x))))
        })
      },
      xform = function(x, par) {
        x^par$shape
      },
      notes = {
        if (!is.null(shape) && shape > 1 && pow != 1) {
          message(
            "\n Note: Shape > 1 and pow != 1. Transforming data by data^shape to conduct energy GOF test.\n"
          )
        }
      }
    ),
    class = c("ParetoDist", "GeneralizedGOFDist", "GOFDist")
  )
  validate_par(dist)
  ## Specific to Pareto: transform the params if necessary.
  dist <- pareto_set_sampler_par(dist)
  set_composite_class(dist)
}

pareto_set_sampler_par <- function(dist) {
  ## X ~ P(scale, shape) -> X^shape ~ P(newscale = scale^shape, newshape = 1)
  initpar <- dist$par
  initshape <- initpar$shape
  initscale <- initpar$scale
  initpow <- initpar$pow
  ## New par
  xshape <- 1
  xscale <- initscale^initshape
  xpow <- ifelse(initpow < xshape, initpow, 0.5)
  xpar <- list(scale = xscale, shape = xshape, pow = xpow)
  if (!is.null(dist$par$shape)) {
    if (dist$par$shape > 1 && (dist$par$pow != 1)) {
      dist$sampler_par <- xpar
    } else {
      dist$sampler_par <- dist$par
    }
  }
  validate_par(dist)
  dist
}


##### Cauchy

#' @title Create a Cauchy distribution object for energy testing
#' @inherit normal_dist return author
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the generalized energy goodness-of-fit test
#'   against a Cauchy distribution. If `location` and `scale` are both NULL,
#'   perform a composite test.
#' @param location NULL, or same as in [stats::rcauchy()]
#' @param scale NULL, or same as in [stats::rcauchy()]
#' @param pow Optionally set the exponent of the energy test. 0 < pow < 1 is
#'   required for the Cauchy distribution. Default is 0.5.
#' @examples
#' d <- cauchy_dist(4, 4)
#' x <- rcauchy(10, 4, 4)
#' egofd(x, d, 0)
#'
#' @export
cauchy_dist <- function(location = NULL, scale = NULL, pow = 0.5) {
  dist <- structure(
    list(
      name = "Cauchy",
      composite_p = is_composite(location, scale),
      par = list(location = location, scale = scale, pow = pow),
      sampler_par = list(location = 0, scale = 1, pow = 0.5),
      support = function(x, par) {
        all(is.finite(x))
      },
      par_domain = function(par) {
        all(
          (par$scale > 0 && length(par$scale) == 1) || is.null(par$scale),
          (length(par$location) == 1) || is.null(par$location),
          par$pow < 1,
          par$pow > 0
        )
      },
      sampler = function(n, par) {
        rcauchy(n, location = 0, scale = 1)
      },
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
      },
      statistic = function(x) {
        as.list(fitdistrplus::mledist(x, "cauchy")$estimate)
      }
    ),
    class = c("CauchyDist", "GeneralizedGOFDist", "GOFDist")
  )
  validate_par(dist)
  set_composite_class(dist)
}

##### Stable

#' @title Create a stable distribution object for energy testing
#' @inherit normal_dist return author
#' @param location Same as in [stats::rcauchy()]
#' @param scale Same as in [stats::rcauchy()]
#' @param skew -1 < skew < 1 is required
#' @param stability The tail index or stability index. Controls the fatness of
#'   the tails. 0 < stability <= 2 is required.
#' @description Create an S3 object that sets all the required data needed by
#'   energyGOFdist to execute the generalized energy goodness-of-fit test against a
#'   stable distribution. Only simple tests are supported.
#' @param pow Exponent of the energy test. 0 < stability/2 is required.
#' @details This is a very slow test due to an onerous amount of numerical
#'   integration required.
#' @noRd
## stable_dist <- function(location = 0, scale = 1,
##                         skew = 0, stability = 1,
##                         pow = stability / 4) {
##   dist <- structure(
##     list(
##       name = "Stable",
##       composite_p = composite_not_allowed(location, scale, skew, stability),
##       par = list(location = location, scale = scale, skew = skew,
##                  stability = stability, pow = pow),
##       sampler_par = list(location = 0, scale = 1, skew = skew,
##                          stability = stability, pow = pow),
##       support = function(x, par) {
##         if (par$stability < 1 && par$skew == 1) {
##           all(x > par$location, is.finite(x))
##         } else if (par$stability < 1 && par$skew == -1) {
##           all(x < par$location, is.finite(x))
##         } else
##           all(is.finite(x))
##       },
##       par_domain = function(par) {
##         all(
##           length(par$stability) == 1,
##           length(par$skew) == 1,
##           length(par$scale) == 1,
##           length(par$location) == 1,
##           par$stability > 0 && par$stability <= 2,
##           is.finite(par$location),
##           par$scale > 0,
##           par$skew >= -1 && par$skew <= 1,
##           par$pow < par$stability / 2,
##           par$pow > 0
##         )
##       },
##       sampler = function(n, par) {
##         d <- par$location
##         s <- par$scale
##         b <- par$skew
##         a <- par$stability
##         stabledist::rstable(n, alpha = a, beta = b,
##                             gamma = s, delta = d)
##       },
##       EXYhat = function(x, par) {
##         n <- length(x)
##         a <- par$stability
##         b <- par$skew
##         pow <- par$pow
##         if (a == 1 && b != 0) {
##           # skewed and cauchy-like
##           A <- 2 / pi * gamma(pow + 1)
##           B <- sinpi(pow / 2)
##           integrand <- function(t, x, par) {
##             a <- par$stability
##             b <- par$skew
##             pow <- par$pow
##             (1 - exp(-t^a) * cos(b * t^a * log(t) + x * t)) / t^(pow + 1)
##           }
##           I <- 0
##           for (i in 1:n) {
##             if (x[i] > 2000) {
##               I <- I + abs(x[i])^pow
##             } else {
##               I <- I + integrate(integrand, 0, 1e5, x = x[i], par = par,
##                                  subdivisions = 1000)$value
##             }
##             return (A * B * I / n)
##           }
##         } else if (a != 1 && b != 0) {
##           # General Asymmetric case
##           A <- 2 / pi * gamma(pow + 1)
##           B <- sinpi(pow / 2)
##           integrand <- function(t, x, par, pow) {
##             a <- par$stability
##             b <- par$skew
##             (1 - exp(-t^a) * cos(b * t^a * tanpi(a / 2) - x * t)) / t^(pow + 1)
##           }
##           I <- 0
##           for (i in 1:n) {
##             I <- if (x[i] > 2000) {
##               I + abs(x[i])
##             } else {
##               I + integrate(integrand, 0, Inf, x = x[i], par = par,
##                             subdivisions = 1000)$value
##             }
##           }
##           return(A * B * I / n)
##         } else if (a != 1 && b == 0){
##           #General Symmetric Case
##           integrand <- function(t, x, par) {
##             a <- par$stability
##             pow <- par$pow
##             (1 - exp(-t^a) * cos(x * t)) / t^(pow + 1)
##           }
##           I <- 0
##           for (i in 1:n) {
##             I <- if (x[i] > 2000) {
##               I + abs(x[i])
##             } else {
##               I + integrate(integrand, 0, 1e5, x = x[i],
##                             par = par, subdivisions = 1000,
##                             rel.tol = 1e-4, abs.tol = 1e-4)$value
##             }
##           }
##           return(I / n)
##         } else {
##           # Cauchy Case: a == 1 and b == 0
##           mean((1 + x^2)^(pow / 2) * cos(pow * atan(x)) / cospi(pow / 2))
##         }
##       },
##       EYY = function(par) {
##         a <- par$stability
##         b <- par$skew
##         pow <- par$pow
##         2^(pow / a + 1) * gamma(1 - pow / a) *
##           gamma(pow) * sinpi(pow / 2) / pi
##       },
##       xform = function(x, par) {
##         # Must transform in simple case.
##         d <- par$location
##         s <- par$scale
##         (x - d) / s
##       }
##     ), class = c("StableDist", "GeneralizedGOFDist", "SimpleGOFDist", "GOFDist")
##   )
##   validate_par(dist)
##   dist
## }
