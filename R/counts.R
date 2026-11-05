## counts.R -- moments of, and calibration to, the zero-truncated count distributions used
## for the mosquito uptake / inoculation draws (LAMBDA_OOCYSTS / K_OOCYSTS, LAMBDA_SPOROZOITES).

#' Realised mean and variance of a zero-truncated count distribution
#'
#' `LAMBDA_OOCYSTS` / `LAMBDA_SPOROZOITES` are **pre-truncation** rates, so they are *not*
#' the mean count you would read off a paper: conditioning away the zeros pulls the realised
#' mean up. This returns what the model actually produces for a given `lambda` (and dispersion
#' `k`), so you can check a parameter choice against a reported figure.
#'
#' With base mean \eqn{\lambda}, dispersion \eqn{k} and zero-probability \eqn{p_0}, the
#' truncated moments are \eqn{E[X] = \lambda/(1-p_0)} and
#' \eqn{Var[X] = E[X^2]/(1-p_0) - E[X]^2}. `k <= 0` is Poisson
#' (\eqn{p_0 = e^{-\lambda}}); `k > 0` is negative binomial with
#' \eqn{Var = \lambda + \lambda^2/k} and \eqn{p_0 = (k/(k+\lambda))^k}.
#'
#' @param lambda pre-truncation mean (> 0).
#' @param k dispersion; `0` (default) = Poisson, `> 0` = negative binomial.
#' @return named numeric of `mean` and `var` (the realised, zero-truncated values).
#' @examples
#' uptake_moments(2)          # Poisson lambda = 2 -> realised mean ~2.31
#' uptake_moments(2, k = 1)   # NB, same lambda   -> realised mean ~3.0, far more spread
#' @seealso [uptake_calibrate()] for the inverse (target moments -> lambda, k).
#' @export
uptake_moments <- function(lambda, k = 0) {
  if (lambda <= 0) return(c(mean = 1, var = 0))
  ## q = 1 - p0, by expm1() and log1p(): 1 - exp(...) rounds to 0 once lambda or lambda/k is
  ## below about 1e-16. k log(1 + lambda/k) is written as lambda log1p(r)/r, whose ratio tends to 1
  ## (the Poisson limit) where r = lambda/k underflows.
  if (k <= 0) { q <- -expm1(-lambda); ex2 <- lambda + lambda^2 }
  else {
    r <- lambda / k
    q <- -expm1(-lambda * (if (r < 1e-8) 1 - r / 2 else log1p(r) / r))
    ex2 <- lambda + lambda^2 / k + lambda^2
  }
  m <- lambda / q
  c(mean = unname(m), var = unname(ex2 / q - m^2))
}

#' Convert a target realised mean (and variance) into lambda and k
#'
#' The inverse of [uptake_moments()]: give it the mean -- and optionally the variance -- you want
#' the model to *produce* (e.g. a mean oocyst count per infected mosquito reported in a
#' feeding study) and it returns the `LAMBDA_*` / `K_OOCYSTS` values to set.
#'
#' If `var` is omitted, or is not larger than the Poisson variance at that mean, a Poisson
#' (`k = 0`) solution is returned -- the negative binomial can only *add* dispersion. The mean is
#' always matched. A zero-truncated count with a given mean can be only so dispersed: its variance
#' is largest as `k` approaches 0 (about 11.1 at a mean of 3, 271 at a mean of 10), and
#' `k = 0.001`, the smallest the model draws, comes within 0.5% of that limit. A larger target
#' variance returns `k = 0.001` with a warning.
#'
#' @param mean target realised (zero-truncated) mean; must be > 1 and at most 1000, the largest
#'   mean [vivax_params()] accepts.
#' @param var optional target realised variance.
#' @return list with `lambda`, `k`, and the `achieved` moments from [uptake_moments()].
#' @examples
#' uptake_calibrate(2.5)                 # Poisson: lambda giving a realised mean of 2.5
#' uptake_calibrate(2.5, var = 6)        # overdispersed: lambda + k
#' @export
uptake_calibrate <- function(mean, var = NULL) {
  if (mean <= 1) stop("a zero-truncated count has mean > 1; got ", mean, call. = FALSE)
  if (mean > 1000)
    stop("the model accepts a mean of at most 1000 per uptake or bite; got ", mean, call. = FALSE)
  ## The lambda giving this mean at dispersion k. The truncated mean rises with lambda, from 1 as
  ## lambda -> 0, and is at least lambda, so the root lies in (0, mean].
  lambda_at <- function(k)
    stats::uniroot(function(l) uptake_moments(l, k)[["mean"]] - mean, c(1e-300, mean),
                   tol = 1e-15 * mean)$root
  lam0 <- lambda_at(0)
  if (is.null(var) || var <= uptake_moments(lam0, 0)[["var"]] + 1e-9)
    return(list(lambda = lam0, k = 0, achieved = uptake_moments(lam0, 0)))
  ## With the mean held, the variance falls as k rises, so k is a one-dimensional search between
  ## 0.001, the most dispersed count the model draws, and 1e6, which is all but Poisson.
  var_at <- function(logk) uptake_moments(lambda_at(exp(logk)), exp(logk))[["var"]]
  lo <- log(1e-3); hi <- log(1e6)
  k <- if (var >= var_at(lo)) 1e-3 else if (var <= var_at(hi)) 1e6 else
    exp(stats::uniroot(function(lk) log(var_at(lk) / var), c(lo, hi), tol = 1e-10)$root)
  lam <- lambda_at(k); ach <- uptake_moments(lam, k)
  if (abs(ach[["var"]] / var - 1) > 0.01)
    warning("var = ", var, " is more than a zero-truncated count with mean ", mean, " can have: ",
            "the most dispersed count the model draws (k = 0.001) has variance ",
            signif(ach[["var"]], 4), ". Returning that; see $achieved.", call. = FALSE)
  list(lambda = lam, k = k, achieved = ach)
}
