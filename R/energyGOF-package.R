### energyGOF: Goodness-of-fit tests via the energy of data

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


#' @title energyGOF: Goodness-of-fit tests via the energy of data
#' @author John T. Haman
#'
#' @description Perform a parametric goodness-of-fit test of univariate data.
#'   Both simple (known parameter) and composite (unknown parameter) cases are
#'   supported, but not all distributions allow for a composite test. In all
#'   cases, *p*-values are determined via parametric bootstrap. The energyGOF
#'   package supporting testing certain distribution families that lack finite
#'   first absolute moment (currently just Cauchy and Pareto).
#'
#' @section Getting Started:
#'
#' The main entry point is [energyGOF.test], but you may alternatively use
#'   [energyGOF], which is a different, non-standard interface.
#'
#' Here is a simple example
#'
#' x <- rnorm(10)
#'
#' ## Composite energy goodness-of-fit test (test for Normality with unknown
#' ## parameters)
#'
#' energyGOF.test(x, "normal", nsim = 10)
#'
#' ## Simple energy goodness-of-fit test (test for Normality with known
#' ## parameters). eg is an alias for energyGOF.test.
#'
#' egof.test(x, "normal", nsim = 10, mean = 0, sd = 1)
#'
#'
#' @section Distributions Supported:
#'  The following distributions are supported
#'
#' \tabular{llll}{
#'   \strong{Distribution} \tab \strong{Function} \tab \strong{Parameter} \tab \strong{Composite Allowed} \cr
#'   Normal                \tab normal_dist            \tab mean, sd                              \tab Yes\cr
#'   Uniform               \tab uniform_dist           \tab min, max                              \tab No \cr
#'   Exponential           \tab exponential_dist       \tab rate                                  \tab Yes\cr
#'   Poisson               \tab poisson_dist           \tab lambda                                \tab Yes\cr
#'   Bernoulli             \tab bernoulli_dist         \tab prob                                  \tab No \cr
#'   Binomial              \tab binomial_dist          \tab prob                                  \tab Yes\cr
#'   Beta                  \tab beta_dist              \tab shape1, shape2                        \tab No \cr
#'   Half-Normal           \tab halfnormal_dist        \tab theta                                 \tab No \cr
#'   Laplace               \tab laplace_dist           \tab location, scale                       \tab No \cr
#'   Lognormal             \tab lognormal_dist         \tab meanlog, sdlog                        \tab No \cr
#'   Asymmetric Laplace    \tab asymmetriclaplace_dist \tab location, scale, skew                 \tab No \cr
#'   F                     \tab d_dist                 \tab df1, df2                              \tab No \cr
#'   Weibull               \tab weibull_dist           \tab shape, scale                          \tab No \cr
#'   Gamma                 \tab gamma_dist             \tab shape, rate                           \tab No \cr
#'   Chi Squared           \tab chisq_dist             \tab df                                    \tab No \cr
#'   Inverse Gaussian      \tab inversegaussian_dist   \tab mean, shape                           \tab No \cr
#'   Pareto                \tab pareto_dist            \tab scale, shape, pow, r                  \tab No \cr
#'   Cauchy                \tab cauchy_dist            \tab location, scale, pow                  \tab No \cr
#'   Stable                \tab stable_dist            \tab location, scale, skew, stability, pow \tab No
#' }
#'
#' Note: Although it may be proper to hyphenate some distribution names in
#' text, there are no hyphens in any `*_dist` function names.
#'
#'
#' @section About Energy:
#'
#' Székely, G. J., & Rizzo, M. L. (2023) provide the motivation:
#'
#' "Data energy is a real number (typically a non-negative number) that depends
#' on the distances between data. This concept is based on the notion of
#' Newton’s gravitational potential energy, which is also a function of the
#' distance between bodies. The idea of data energy or energy statistics is to
#' consider statistical observations (data) as heavenly bodies governed by the
#' potential energy of data, which is zero if and only if an underlying
#' statistical hypothesis is true."
#'
#' The notation X' indicates that X' is an independent and identically
#' distributed copy of X.
#'
#' If X and Y are independent and E(|X|^s + |Y|^s) is finite, then for 0 < s <
#' 2,
#'
#' \deqn{2E|X-Y|^s - E|X-X'|^s - E|Y-Y'|^s \ge 0.}
#'
#' Equality is attained if and only if X and Y are identically distributed. The
#' left side of the equation is the energy between X and Y. Energy can be
#' generalized to multivariate data and even more exotic data types, but in
#' this R package, we only treat univariate data.
#'
#' The concept of data energy between two random variables can be adapted to
#' the one-sample goodness-of-fit problem. The one-sample *s*-energy is
#'
#' \deqn{E^* = \frac{2}{n} \sum_i E|x_i - Y|^s - E|Y-Y'|^s - \frac{1}{n^2}
#' \sum_i \sum_j |x_i - x_j|^s,}
#'
#' when \eqn{0 < s < 2} and \eqn{E|X|^s, E|Y|^s < \infty.}
#'
#' In most tests in the energyGOF package s = 1. In some cases (Pareto, Cauchy,
#' Stable), E|Y| is not finite, so we need to use an s < 1. This is done by
#' passing `pow` into `...` (but in all tests a default `pow` is provided).
#' These tests are called generalized energy goodness-of-fit tests.
#'
#' In the one-sample goodness-of-fit regime, we wonder if \eqn{x_i ~ X} (where
#' the distribution of X is hidden) follows the same distribution as Y, which
#' is specified. If X and Y have the same distribution, then \eqn{Q = nE^*} is
#' a quadratic form of centered Gaussian random variables with expected value
#' \eqn{E|Y-Y'|^s}. If X and Y differ, then Q goes to Inf. So, Q provides a
#' consistent goodness-of-fit test, even in some situations where E|Y| is not
#' finite. And that's what energyGOF.test does. Asymptotic theory of
#' V-statistics can be applied to prove that tests based on Q are statistically
#' consistent goodness-of-fit tests.
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
#' Yang, G. (2012). The Energy Goodness-of-fit Test for Univariate Stable
#' Distributions (Doctoral dissertation, Bowling Green State University).
#'
#'
#' @docType package
#' @name energyGOF
#' @keywords internal
"_PACKAGE"
