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
  if (k <= 0) { p0 <- exp(-lambda); ex2 <- lambda + lambda^2 }
  else        { p0 <- (k / (k + lambda))^k; ex2 <- lambda + lambda^2 / k + lambda^2 }
  m <- lambda / (1 - p0)
  c(mean = m, var = ex2 / (1 - p0) - m^2)
}

#' Convert a target realised mean (and variance) into lambda and k
#'
#' The inverse of [uptake_moments()]: give it the mean -- and optionally the variance -- you want
#' the model to *produce* (e.g. a mean oocyst count per infected mosquito reported in a
#' feeding study) and it returns the `LAMBDA_*` / `K_OOCYSTS` values to set.
#'
#' If `var` is omitted, or is not larger than the Poisson variance at that mean, a Poisson
#' (`k = 0`) solution is returned -- the negative binomial can only *add* dispersion.
#'
#' @param mean target realised (zero-truncated) mean; must be > 1.
#' @param var optional target realised variance.
#' @return list with `lambda`, `k`, and the `achieved` moments from [uptake_moments()].
#' @examples
#' uptake_calibrate(2.5)                 # Poisson: lambda giving a realised mean of 2.5
#' uptake_calibrate(2.5, var = 6)        # overdispersed: lambda + k
#' @export
uptake_calibrate <- function(mean, var = NULL) {
  if (mean <= 1) stop("a zero-truncated count has mean > 1; got ", mean, call. = FALSE)
  ## Poisson branch: one unknown, solve the mean directly (monotone in lambda).
  lam0 <- stats::uniroot(function(l) uptake_moments(l, 0)["mean"] - mean,
                         interval = c(1e-10, 20 * mean), tol = 1e-12)$root
  if (is.null(var) || var <= uptake_moments(lam0, 0)["var"] + 1e-9)
    return(list(lambda = lam0, k = 0, achieved = uptake_moments(lam0, 0)))
  ## NB branch: two unknowns, two moments. Solved on the log scale (a nested root-find
  ## degenerates as k -> 0, where the truncated mean stops depending on lambda).
  obj <- function(par) {
    a <- uptake_moments(exp(par[1]), exp(par[2]))
    (log(a[["mean"]] / mean))^2 + (log(a[["var"]] / var))^2
  }
  fit <- stats::optim(c(log(lam0), log(max(mean, 1))), obj, method = "L-BFGS-B",
                      lower = c(log(1e-6), log(1e-3)), upper = c(log(1e4), log(1e6)))
  lam <- exp(fit$par[1]); k <- exp(fit$par[2]); ach <- uptake_moments(lam, k)
  if (abs(ach[["mean"]] / mean - 1) > 0.01 || abs(ach[["var"]] / var - 1) > 0.01)
    warning("could not match the target moments closely; see $achieved", call. = FALSE)
  list(lambda = lam, k = k, achieved = ach)
}
