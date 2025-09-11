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
#'
#' @section Getting Started:
#' The main entry point is energyGOF.test. The only documentation you need to
#' read is [energyGOF.test] and [energyGOF-package].
#'
#'
#' Here is a simple example to get you going
#'
#' ```
#' x <- rnorm(10)
#'
#' ## Composite energy goodness-of-fit test (test for Normality with unknown
#' ## parameters)
#'
#' energyGOF.test(x, "normal", nsim = 1e5)
#'
#' ## Simple energy goodness-of-fit test (test for Normality with known
#' ## parameters). egof.test is an alias for energyGOF.test.
#'
#' egof.test(x, "normal", nsim = 1e5, mean = 0, sd = 1)
#'
#'
#' ## Two-sample test
#' y <- rt(10, 1)
#' egof.test(x, y, nsim = 1e5)
#'
#' ## Test agaist any distribution function by transforming data to uniform
#' egof.test(y, pt, nsim = 1e5)
#'
#' ```
#'
#' You may alternatively use the [energyGOFdist] function, which is a different
#' interface using S3 objects, but it provides the same result. There is a lot
#' of documentation in this package for the various S3 constructors that are
#' needed by [energyGOFdist], **BUT** if you just want to do some testing and
#' use the standard interface, you can probably ignore all of that and just
#' read the page for [energyGOF.test].
#'
#' @section Distributions Supported:
#'  The following distributions are supported.
#'
#' ## TODO: Update table
#'
#' \tabular{llll}{
#'   \strong{Distribution} \tab \strong{Function} \tab \strong{Parameter} \tab \strong{Composite Allowed} \cr
#'   Normal                \tab `normal_dist`            \tab mean, sd                              \tab Yes\cr
#'   Uniform               \tab `uniform_dist`           \tab min, max                              \tab No \cr
#'   Exponential           \tab `exponential_dist`       \tab rate                                  \tab Yes\cr
#'   Poisson               \tab `poisson_dist`           \tab lambda                                \tab Yes\cr
#'   Bernoulli             \tab `bernoulli_dist`         \tab prob                                  \tab No \cr
#'   Binomial              \tab `binomial_dist`          \tab prob                                  \tab No\cr
#'   Beta                  \tab `beta_dist`              \tab shape1, shape2                        \tab Yes \cr
#'   Half-Normal           \tab `halfnormal_dist`        \tab theta                                 \tab No \cr
#'   Laplace               \tab `laplace_dist`           \tab location, scale                       \tab Yes \cr
#'   Lognormal             \tab `lognormal_dist`         \tab meanlog, sdlog                        \tab Yes \cr
#'   Asymmetric Laplace    \tab `asymmetriclaplace_dist` \tab location, scale, skew                 \tab Yes \cr
#'   F                     \tab `f_dist`                 \tab df1, df2                              \tab No \cr
#'   Weibull               \tab `weibull_dist`           \tab shape, scale                          \tab Yes \cr
#'   Gamma                 \tab `gamma_dist`             \tab shape, rate                           \tab No \cr
#'   Chi Squared           \tab `chisq_dist`             \tab df                                    \tab No \cr
#'   Inverse Gaussian      \tab `inversegaussian_dist`   \tab mean, shape                           \tab Yes \cr
#'   Pareto                \tab `pareto_dist`            \tab scale, shape, pow, r                  \tab Yes \cr
#'   Cauchy                \tab `cauchy_dist`            \tab location, scale, pow                  \tab Yes \cr
#' }
#'
#'
#' @section Simple and Composite Testing:
#'
#' There are two types of goodness-of-fit tests covered by the energyGOF
#' package, *simple* and *composite*. It's important to know the difference
#' because they yield different results. Simple GOF tests test the data `x`
#' against a specific distribution with _known parameters_ that you must pass
#' to energyGOF.test in the ellipsis argument (...). You should use a simple
#' GOF test if you wish to test questions like
#' "my data is Normal with mean 1 and sd 2". `energyGOF` can also conduct
#' _some_ composite GOF tests. A composite test is performed if no parameters
#' are passed in the ellipsis argument (...). You should conduct a composite
#' test if your research question is
#' "my data are Normal, but I don't know what the parameters are." Obviously,
#' this composite question is much more common in practice.
#'
#' All the composite tests in energyGOF assume that *none* of the parameters
#' are known. So while there is a statistical test of Normality with known mean
#' and unknown sd, this is not implemented in the energyGOF package, you will
#' get an error is you try ask for that. So, either pass all the distribution
#' parameters or none of them. (In the special case of the Normal distribution,
#' you can use the energy package to test the GOF hypothesis with any
#' combination of known and unknown parameters.)
#'
#' For each test, `energyGOF.test` calculates the test statistic and a
#' *p*-value. In all cases the *p*-value is calculated via parametric
#' bootstrap. For large `nsim`, the *p*-values should be reasonably honest in
#' small-ish samples. You should set nsim to be a very large number in
#' practice. I recommend at least 10,000.
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
#' The notation \eqn{X'} indicates that \eqn{X'} is an independent and
#' identically distributed copy of \eqn{X}.
#'
#' If \eqn{X} and \eqn{Y} are independent and \eqn{E(|X|^s + |Y|^s)} is finite,
#' then for \eqn{0 < s < 2},
#'
#' \deqn{2E|X-Y|^s - E|X-X'|^s - E|Y-Y'|^s \ge 0.}
#'
#' Equality is attained if and only if \eqn{X} and \eqn{Y} are identically
#' distributed. The left side of the equation is the energy between \eqn{X} and
#' \eqn{Y}. Energy can be generalized to multivariate data and even more exotic
#' data types, but in this R package, we only treat univariate data.
#'
#' The concept of data energy between two random variables can be adapted to
#' the one-sample goodness-of-fit problem. The one-sample \eqn{s}-energy is
#'
#' \deqn{E^* = \frac{2}{n} \sum_i E|x_i - Y|^s - E|Y-Y'|^s - \frac{1}{n^2}
#' \sum_i \sum_j |x_i - x_j|^s,}
#'
#' when \eqn{0 < s < 2} and \eqn{E|X|^s, E|Y|^s < \infty.}
#'
#' In most tests in the `energyGOF` package \eqn{s = 1}. In some cases (Pareto
#' and Cauchy), \eqn{E|Y|} is not finite, so we need to use an \eqn{s < 1}.
#' This is done by passing `pow` into `...` (but in all tests a default `pow`
#' is provided). These tests are called *generalized* energy goodness-of-fit
#' tests in this package as well as in Székely, G. J., & Rizzo, M. L. (2023).
#'
#' To connect energy back to GOF testing, in the one-sample goodness-of-fit
#' regime, we test if a sample \eqn{x_1, \ldots, x_n \sim X} (where the
#' distribution of \eqn{X} is hidden) follows the same distribution as \eqn{Y},
#' which is specified. If \eqn{X} and \eqn{Y} have the same distribution, then
#' the distribution of \eqn{Q = nE^*} is a quadratic form of centered Gaussian
#' random variables with expected value \eqn{E|Y-Y'|^s}. If \eqn{X} and \eqn{Y}
#' differ, then \eqn{Q \to \infty} with \eqn{n}. So, \eqn{Q} provides a
#' consistent goodness-of-fit test, even in some situations where \eqn{E|Y|} is
#' not finite. And that's what `energyGOF.test` does. Asymptotic theory of
#' V-statistics can be applied to prove that tests based on \eqn{Q} are
#' statistically consistent goodness-of-fit tests.
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
#' @aliases energyGOF-package

"_PACKAGE"
