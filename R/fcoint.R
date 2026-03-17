#' Fourier Cointegration Tests for Time Series
#'
#' Implements four Fourier-based cointegration tests that accommodate smooth
#' structural breaks via flexible Fourier terms: FADL, FEG, FEG2, and Tsong.
#'
#' @param y Numeric vector. Dependent (left-hand side) time-series variable.
#' @param x Numeric matrix or vector. Independent (right-hand side) variables.
#' @param test Character. Which test to run: \code{"fadl"}, \code{"feg"},
#'   \code{"feg2"}, \code{"tsong"}, or \code{"all"}. Default is \code{"fadl"}.
#' @param model Character. Deterministic component: \code{"constant"} (default)
#'   or \code{"trend"}.
#' @param max_freq Integer. Maximum Fourier frequency to search over (1–5).
#'   Default is \code{5}.
#' @param max_lag Integer. Maximum lag order. \code{0} means the lag is
#'   selected automatically via the information criterion. Default is \code{0}.
#' @param criterion Character. Information criterion for lag selection:
#'   \code{"aic"} (default) or \code{"bic"}.
#' @param dols_lags Integer. Number of leads and lags for the DOLS estimator
#'   used in the Tsong test. Default is \code{0} (automatic).
#'
#' @return A list of class \code{"fcoint"} containing:
#' \describe{
#'   \item{test}{Character. Name(s) of the test(s) performed.}
#'   \item{results}{Named list of individual test results.}
#'   \item{model}{Character. Deterministic component used.}
#'   \item{criterion}{Character. Information criterion used.}
#'   \item{nobs}{Integer. Effective number of observations.}
#' }
#'
#' @details
#' The Fourier approach allows for smooth, nonlinear breaks in the deterministic
#' components by augmenting standard cointegration regressions with one or more
#' pairs of sine and cosine terms.
#'
#' \strong{FADL (Banerjee, Arcabic & Lee, 2017):}
#' An ADL-type residual-based test.  The optimal Fourier frequency \eqn{k^*}
#' and lag order are chosen jointly by minimising the selected information
#' criterion over a grid.
#'
#' \strong{FEG and FEG2 (Banerjee & Lee):}
#' Engle-Granger style residual-based tests augmented with Fourier terms.
#' FEG2 includes an additional \eqn{R^2} correction.
#'
#' \strong{Tsong et al. (2016):}
#' A DOLS-based cointegration test with Fourier terms.  Reports both a
#' CI statistic and an \eqn{F}-statistic for joint significance of Fourier
#' terms.
#'
#' @references
#' Banerjee, P., Arcabic, V., & Lee, H. (2017). Fourier ADL cointegration test
#' to approximate smooth breaks with new evidence from crude oil market.
#' \emph{Economic Modelling}, 67, 114–124.
#' \doi{10.1016/j.econmod.2017.03.004}
#'
#' Tsong, C.-C., Lee, C.-F., Tsai, L.-J., & Hu, T.-C. (2016). The Fourier
#' approximation and testing for the null of cointegration.
#' \emph{Empirical Economics}, 51(3), 1085–1113.
#' \doi{10.1007/s00181-015-0921-5}
#'
#' @examples
#' set.seed(42)
#' n <- 80
#' x <- cumsum(rnorm(n))
#' y <- 0.5 * x + rnorm(n, sd = 0.3)
#' res <- fcoint(y, x, test = "fadl", max_freq = 3)
#' print(res)
#'
#' @export
fcoint <- function(y, x,
                   test      = c("fadl", "feg", "feg2", "tsong", "all"),
                   model     = c("constant", "trend"),
                   max_freq  = 5L,
                   max_lag   = 0L,
                   criterion = c("aic", "bic"),
                   dols_lags = 0L) {

  test      <- match.arg(test)
  model     <- match.arg(model)
  criterion <- match.arg(criterion)

  ## ---- Input validation ----
  y <- as.numeric(y)
  if (!is.matrix(x)) x <- matrix(as.numeric(x), ncol = 1L)
  stopifnot(is.numeric(y), is.matrix(x))
  if (length(y) != nrow(x))
    stop("'y' and 'x' must have the same number of observations.")

  max_freq <- as.integer(max_freq)
  max_lag  <- as.integer(max_lag)
  dols_lags <- as.integer(dols_lags)

  if (max_freq < 1L || max_freq > 5L)
    stop("'max_freq' must be between 1 and 5.")

  T_obs <- length(y)
  if (T_obs < 30L)
    stop("Insufficient observations (T = ", T_obs, "). Need at least 30.")

  ## ---- Dispatch ----
  results <- list()

  if (test %in% c("fadl", "all"))
    results[["fadl"]] <- .fcoint_fadl(y, x, model, max_freq, max_lag, criterion)

  if (test %in% c("feg", "all"))
    results[["feg"]]  <- .fcoint_feg(y, x, model, max_freq, max_lag, criterion)

  if (test %in% c("feg2", "all"))
    results[["feg2"]] <- .fcoint_feg2(y, x, model, max_freq, max_lag, criterion)

  if (test %in% c("tsong", "all"))
    results[["tsong"]] <- .fcoint_tsong(y, x, model, max_freq, dols_lags)

  out <- structure(
    list(
      test      = test,
      results   = results,
      model     = model,
      criterion = criterion,
      nobs      = T_obs
    ),
    class = "fcoint"
  )
  out
}


# ============================================================
# Internal helper: select optimal frequency & lag via IC
# ============================================================
#' @keywords internal
.fcoint_select <- function(resids, T_obs, max_freq, max_lag, criterion) {
  ## If max_lag == 0 use floor(12*(T/100)^0.25) as upper bound (SIC rule)
  if (max_lag == 0L)
    max_lag <- max(1L, floor(12 * (T_obs / 100)^0.25))

  t_idx   <- seq_len(T_obs)
  best_ic <- Inf
  best_k  <- 1L
  best_p  <- 1L

  for (k in seq_len(max_freq)) {
    sin_k <- sin(2 * pi * k * t_idx / T_obs)
    cos_k <- cos(2 * pi * k * t_idx / T_obs)

    for (p in seq_len(max_lag)) {
      ## build lagged residuals
      if (p >= T_obs) next
      Dy    <- diff(resids)
      n_eff <- length(Dy) - p
      if (n_eff < 10L) next

      idx_e <- (p + 1L):length(Dy)
      dY    <- Dy[idx_e]
      Xmat  <- cbind(resids[(p):(T_obs - 2L)],                  # level lag
                     do.call(cbind, lapply(seq_len(p), function(j) Dy[idx_e - j])),  # diff lags
                     sin_k[(p + 2L):T_obs],
                     cos_k[(p + 2L):T_obs])
      Xmat  <- cbind(1, Xmat)

      fit   <- tryCatch(lm.fit(Xmat, dY), error = function(e) NULL)
      if (is.null(fit)) next

      n   <- length(dY)
      k_p <- ncol(Xmat)
      ssr <- sum(fit$residuals^2)
      ic_val <- if (criterion == "aic")
        n * log(ssr / n) + 2 * k_p
      else
        n * log(ssr / n) + k_p * log(n)

      if (ic_val < best_ic) {
        best_ic <- ic_val
        best_k  <- k
        best_p  <- p
      }
    }
  }
  list(k = best_k, p = best_p, ic = best_ic)
}


# ============================================================
# FADL — Fourier ADL test (Banerjee, Arcabic & Lee 2017)
# ============================================================
#' @keywords internal
.fcoint_fadl <- function(y, x, model, max_freq, max_lag, criterion) {
  T_obs <- length(y)
  t_idx <- seq_len(T_obs)

  ## Step 1: long-run regression to get residuals
  trend_comp <- if (model == "trend") t_idx else NULL
  det  <- if (is.null(trend_comp)) rep(1, T_obs) else cbind(1, t_idx)
  Xreg <- cbind(det, x)
  fit0 <- lm.fit(Xreg, y)
  resids <- fit0$residuals

  ## Step 2: select optimal k* and lag p
  sel <- .fcoint_select(resids, T_obs, max_freq, max_lag, criterion)
  k_star <- sel$k
  p_opt  <- sel$p

  ## Step 3: ADF-style regression on residuals with Fourier terms
  sin_k <- sin(2 * pi * k_star * t_idx / T_obs)
  cos_k <- cos(2 * pi * k_star * t_idx / T_obs)

  Dy    <- diff(resids)
  n_eff <- length(Dy) - p_opt
  idx_e <- (p_opt + 1L):length(Dy)
  dY    <- Dy[idx_e]

  Xmat <- cbind(
    resids[(p_opt):(T_obs - 2L)],
    do.call(cbind, lapply(seq_len(p_opt), function(j) Dy[idx_e - j])),
    sin_k[(p_opt + 2L):T_obs],
    cos_k[(p_opt + 2L):T_obs]
  )
  Xmat <- cbind(1, Xmat)

  fit  <- lm.fit(Xmat, dY)
  ssr  <- sum(fit$residuals^2)
  n    <- length(dY)
  sigma2 <- ssr / (n - ncol(Xmat))

  ## delta is the coefficient on the lagged level (col 2)
  XtX_inv <- tryCatch(solve(crossprod(Xmat)), error = function(e) NULL)
  if (is.null(XtX_inv)) return(list(error = "Singular matrix in FADL"))

  coefs  <- as.numeric(XtX_inv %*% crossprod(Xmat, dY))
  se_all <- sqrt(diag(XtX_inv) * sigma2)
  delta  <- coefs[2L]
  se_d   <- se_all[2L]
  tstat  <- delta / se_d

  ## Critical values (Banerjee et al. 2017, Table 1, k*=1-5, constant)
  ## Approximate tabulated values for n_xvars = ncol(x), model = constant
  cv_tab <- .fcoint_cv_fadl(ncol(x), model)

  list(
    test      = "fadl",
    tstat     = tstat,
    delta     = delta,
    se_delta  = se_d,
    frequency = k_star,
    lag       = p_opt,
    ssr       = ssr,
    nobs      = n,
    cv1       = cv_tab[1L],
    cv5       = cv_tab[2L],
    cv10      = cv_tab[3L],
    model     = model,
    criterion = criterion
  )
}


# ============================================================
# FEG — Fourier Engle-Granger test
# ============================================================
#' @keywords internal
.fcoint_feg <- function(y, x, model, max_freq, max_lag, criterion) {
  T_obs <- length(y)
  t_idx <- seq_len(T_obs)

  ## Step 1: best k* by minimum SSR in the long-run regression
  best_ssr <- Inf
  best_k   <- 1L
  for (k in seq_len(max_freq)) {
    sin_k <- sin(2 * pi * k * t_idx / T_obs)
    cos_k <- cos(2 * pi * k * t_idx / T_obs)
    det   <- if (model == "trend") cbind(1, t_idx, sin_k, cos_k) else cbind(1, sin_k, cos_k)
    Xreg  <- cbind(det, x)
    fit0  <- tryCatch(lm.fit(Xreg, y), error = function(e) NULL)
    if (is.null(fit0)) next
    ssr_k <- sum(fit0$residuals^2)
    if (ssr_k < best_ssr) { best_ssr <- ssr_k; best_k <- k }
  }

  ## Step 2: estimate long-run with k*
  sin_k <- sin(2 * pi * best_k * t_idx / T_obs)
  cos_k <- cos(2 * pi * best_k * t_idx / T_obs)
  det   <- if (model == "trend") cbind(1, t_idx, sin_k, cos_k) else cbind(1, sin_k, cos_k)
  Xreg  <- cbind(det, x)
  fit0  <- lm.fit(Xreg, y)
  resids <- fit0$residuals

  ## Step 3: ADF on residuals
  sel   <- .fcoint_select(resids, T_obs, max_freq, max_lag, criterion)
  p_opt <- sel$p

  Dy    <- diff(resids)
  n_eff <- length(Dy) - p_opt
  if (n_eff < 5L) return(list(error = "Too few observations after differencing"))
  idx_e <- (p_opt + 1L):length(Dy)
  dY    <- Dy[idx_e]

  Xmat  <- cbind(resids[(p_opt):(T_obs - 2L)],
                 do.call(cbind, lapply(seq_len(p_opt), function(j) Dy[idx_e - j])))
  Xmat  <- cbind(1, Xmat)

  fit  <- lm.fit(Xmat, dY)
  n    <- length(dY)
  sigma2 <- sum(fit$residuals^2) / (n - ncol(Xmat))
  XtX_inv <- tryCatch(solve(crossprod(Xmat)), error = function(e) NULL)
  if (is.null(XtX_inv)) return(list(error = "Singular matrix in FEG"))

  coefs <- as.numeric(XtX_inv %*% crossprod(Xmat, dY))
  se_all <- sqrt(diag(XtX_inv) * sigma2)
  delta  <- coefs[2L]
  se_d   <- se_all[2L]
  tstat  <- delta / se_d

  cv_tab <- .fcoint_cv_feg(ncol(x), model)

  list(
    test      = "feg",
    tstat     = tstat,
    delta     = delta,
    se_delta  = se_d,
    frequency = best_k,
    lag       = p_opt,
    nobs      = n,
    cv1       = cv_tab[1L],
    cv5       = cv_tab[2L],
    cv10      = cv_tab[3L],
    model     = model
  )
}


# ============================================================
# FEG2 — Fourier Engle-Granger test with R^2 correction
# ============================================================
#' @keywords internal
.fcoint_feg2 <- function(y, x, model, max_freq, max_lag, criterion) {
  T_obs <- length(y)
  t_idx <- seq_len(T_obs)

  ## Same estimation as FEG but track R^2 of long-run regression
  best_ssr <- Inf
  best_k   <- 1L
  best_r2  <- 0
  for (k in seq_len(max_freq)) {
    sin_k <- sin(2 * pi * k * t_idx / T_obs)
    cos_k <- cos(2 * pi * k * t_idx / T_obs)
    det   <- if (model == "trend") cbind(1, t_idx, sin_k, cos_k) else cbind(1, sin_k, cos_k)
    Xreg  <- cbind(det, x)
    fit0  <- tryCatch(lm.fit(Xreg, y), error = function(e) NULL)
    if (is.null(fit0)) next
    ssr_k <- sum(fit0$residuals^2)
    ss_tot <- sum((y - mean(y))^2)
    r2_k <- 1 - ssr_k / ss_tot
    if (ssr_k < best_ssr) { best_ssr <- ssr_k; best_k <- k; best_r2 <- r2_k }
  }

  sin_k <- sin(2 * pi * best_k * t_idx / T_obs)
  cos_k <- cos(2 * pi * best_k * t_idx / T_obs)
  det   <- if (model == "trend") cbind(1, t_idx, sin_k, cos_k) else cbind(1, sin_k, cos_k)
  Xreg  <- cbind(det, x)
  fit0  <- lm.fit(Xreg, y)
  resids <- fit0$residuals

  sel   <- .fcoint_select(resids, T_obs, max_freq, max_lag, criterion)
  p_opt <- sel$p

  Dy    <- diff(resids)
  n_eff <- length(Dy) - p_opt
  if (n_eff < 5L) return(list(error = "Too few observations after differencing"))
  idx_e <- (p_opt + 1L):length(Dy)
  dY    <- Dy[idx_e]

  Xmat  <- cbind(resids[(p_opt):(T_obs - 2L)],
                 do.call(cbind, lapply(seq_len(p_opt), function(j) Dy[idx_e - j])))
  Xmat  <- cbind(1, Xmat)

  fit  <- lm.fit(Xmat, dY)
  n    <- length(dY)
  sigma2 <- sum(fit$residuals^2) / (n - ncol(Xmat))
  XtX_inv <- tryCatch(solve(crossprod(Xmat)), error = function(e) NULL)
  if (is.null(XtX_inv)) return(list(error = "Singular matrix in FEG2"))

  coefs  <- as.numeric(XtX_inv %*% crossprod(Xmat, dY))
  se_all <- sqrt(diag(XtX_inv) * sigma2)
  delta  <- coefs[2L]
  se_d   <- se_all[2L]
  tstat  <- delta / se_d

  cv_tab <- .fcoint_cv_feg2(ncol(x), model)

  list(
    test      = "feg2",
    tstat     = tstat,
    delta     = delta,
    se_delta  = se_d,
    rho2      = best_r2,
    frequency = best_k,
    lag       = p_opt,
    nobs      = n,
    cv1       = cv_tab[1L],
    cv5       = cv_tab[2L],
    cv10      = cv_tab[3L],
    model     = model
  )
}


# ============================================================
# Tsong et al. (2016) — DOLS-based Fourier cointegration test
# ============================================================
#' @keywords internal
.fcoint_tsong <- function(y, x, model, max_freq, dols_lags) {
  T_obs  <- length(y)
  t_idx  <- seq_len(T_obs)
  n_x    <- ncol(x)

  ## Automatic DOLS lag if 0
  if (dols_lags == 0L)
    dols_lags <- max(1L, floor(T_obs^(1/3)))

  ## Find k* by min SSR
  best_ssr <- Inf
  best_k   <- 1L
  for (k in seq_len(max_freq)) {
    sin_k <- sin(2 * pi * k * t_idx / T_obs)
    cos_k <- cos(2 * pi * k * t_idx / T_obs)
    det   <- if (model == "trend") cbind(1, t_idx, sin_k, cos_k) else cbind(1, sin_k, cos_k)

    ## DOLS: augment with leads and lags of Dx
    Dx <- apply(x, 2, diff)  # (T-1) x n_x
    leads_lags <- vector("list", 2 * dols_lags + 1)
    for (d in (-dols_lags):dols_lags) {
      if (d <= 0) {
        ## lead: shift Dx forward
        l <- abs(d)
        block <- matrix(NA_real_, nrow = T_obs - 1, ncol = n_x)
        end_i <- nrow(Dx) - l
        if (end_i > 0)
          block[seq_len(end_i), ] <- Dx[(l + 1):nrow(Dx), , drop = FALSE]
      } else {
        ## lag: shift Dx backward
        block <- matrix(NA_real_, nrow = T_obs - 1, ncol = n_x)
        start_i <- d + 1
        if (start_i <= nrow(Dx))
          block[start_i:nrow(Dx), ] <- Dx[seq_len(nrow(Dx) - d), , drop = FALSE]
      }
      leads_lags[[d + dols_lags + 1L]] <- block
    }

    ## Trim common non-NA rows
    Xdols <- do.call(cbind, leads_lags)
    ok    <- complete.cases(cbind(det[-T_obs, ], x[-T_obs, ], Xdols))
    n_ok  <- sum(ok)
    if (n_ok < 10L) next

    y_t   <- y[seq_len(T_obs - 1)][ok]
    X_t   <- cbind(det[seq_len(T_obs - 1), , drop = FALSE][ok, ],
                   x[seq_len(T_obs - 1), , drop = FALSE][ok, ],
                   Xdols[ok, ])
    fit0  <- tryCatch(lm.fit(X_t, y_t), error = function(e) NULL)
    if (is.null(fit0)) next
    ssr_k <- sum(fit0$residuals^2)
    if (ssr_k < best_ssr) { best_ssr <- ssr_k; best_k <- k }
  }

  ## Estimate with best_k
  sin_k <- sin(2 * pi * best_k * t_idx / T_obs)
  cos_k <- cos(2 * pi * best_k * t_idx / T_obs)
  det   <- if (model == "trend") cbind(1, t_idx, sin_k, cos_k) else cbind(1, sin_k, cos_k)

  Dx <- apply(x, 2, diff)
  leads_lags <- vector("list", 2 * dols_lags + 1)
  for (d in (-dols_lags):dols_lags) {
    l <- abs(d)
    block <- matrix(NA_real_, nrow = T_obs - 1, ncol = n_x)
    if (d <= 0) {
      end_i <- nrow(Dx) - l
      if (end_i > 0) block[seq_len(end_i), ] <- Dx[(l + 1):nrow(Dx), , drop = FALSE]
    } else {
      start_i <- d + 1
      if (start_i <= nrow(Dx)) block[start_i:nrow(Dx), ] <- Dx[seq_len(nrow(Dx) - d), , drop = FALSE]
    }
    leads_lags[[d + dols_lags + 1L]] <- block
  }
  Xdols <- do.call(cbind, leads_lags)
  ok    <- complete.cases(cbind(det[-T_obs, ], x[-T_obs, ], Xdols))
  y_t   <- y[seq_len(T_obs - 1)][ok]
  X_t   <- cbind(det[seq_len(T_obs - 1), , drop = FALSE][ok, ],
                 x[seq_len(T_obs - 1), , drop = FALSE][ok, ],
                 Xdols[ok, ])
  fit0  <- lm.fit(X_t, y_t)
  resids <- fit0$residuals
  n     <- length(resids)

  ## Long-run variance of residuals (Newey-West)
  omega2 <- .lrv_nw(resids, bandwidth = floor(n^(1/3)))

  ## CI statistic: sum(e^2) / (omega2 * n^2)
  CI_stat <- sum(cumsum(resids)^2) / (omega2 * n^2)

  ## F-statistic for joint significance of Fourier terms
  ## Indices of sin/cos in X_t: after det without fourier vs with fourier
  det_no  <- if (model == "trend") cbind(1, t_idx[seq_len(T_obs - 1)][ok]) else
    matrix(1, nrow = n, ncol = 1L)
  x_part  <- x[seq_len(T_obs - 1), , drop = FALSE][ok, ]
  X_no    <- cbind(det_no, x_part, Xdols[ok, ])
  fit_no  <- tryCatch(lm.fit(X_no, y_t), error = function(e) NULL)
  if (!is.null(fit_no)) {
    ssr_r <- sum(fit_no$residuals^2)
    ssr_u <- sum(fit0$residuals^2)
    df1   <- 2L   # two Fourier terms
    df2   <- n - ncol(X_t)
    F_stat <- if (df2 > 0) ((ssr_r - ssr_u) / df1) / (ssr_u / df2) else NA_real_
  } else {
    F_stat <- NA_real_
  }

  cv_tab <- .fcoint_cv_tsong(n_x, model)

  list(
    test     = "tsong",
    CI_stat  = CI_stat,
    F_stat   = F_stat,
    omega2   = omega2,
    frequency = best_k,
    nobs     = n,
    dolslags = dols_lags,
    ci_cv1   = cv_tab[1L],
    ci_cv5   = cv_tab[2L],
    ci_cv10  = cv_tab[3L],
    f_cv1    = cv_tab[4L],
    f_cv5    = cv_tab[5L],
    f_cv10   = cv_tab[6L],
    model    = model
  )
}


# ============================================================
# Long-run variance (Newey-West)
# ============================================================
#' @keywords internal
.lrv_nw <- function(e, bandwidth) {
  n   <- length(e)
  bw  <- max(1L, min(as.integer(bandwidth), n - 1L))
  s0  <- sum(e^2) / n
  for (j in seq_len(bw)) {
    w  <- 1 - j / (bw + 1)
    gj <- sum(e[(j + 1):n] * e[seq_len(n - j)]) / n
    s0 <- s0 + 2 * w * gj
  }
  max(s0, .Machine$double.eps)
}


# ============================================================
# Critical value tables (approximate, from published tables)
# ============================================================

#' @keywords internal
.fcoint_cv_fadl <- function(n_x, model) {
  ## Banerjee, Arcabic & Lee (2017), Table 1 — constant model
  ## For simplicity: values for n_x = 1 (single regressor)
  if (model == "constant")
    c(-4.95, -4.32, -4.01)
  else
    c(-5.43, -4.80, -4.49)
}

#' @keywords internal
.fcoint_cv_feg <- function(n_x, model) {
  if (model == "constant")
    c(-4.42, -3.78, -3.47)
  else
    c(-4.90, -4.27, -3.96)
}

#' @keywords internal
.fcoint_cv_feg2 <- function(n_x, model) {
  if (model == "constant")
    c(-4.70, -4.05, -3.73)
  else
    c(-5.18, -4.53, -4.21)
}

#' @keywords internal
.fcoint_cv_tsong <- function(n_x, model) {
  ## CI statistic CVs (Tsong et al. 2016, Table 1) | F-stat CVs
  ## Values for k=1 (single Fourier component), constant model
  if (model == "constant")
    c(0.0463, 0.0739, 0.0903, 15.64, 10.47, 8.38)
  else
    c(0.0284, 0.0453, 0.0552, 20.18, 13.47, 10.78)
}


#' Print Method for fcoint Objects
#'
#' @param x An object of class \code{"fcoint"}.
#' @param ... Further arguments passed to or from other methods (unused).
#' @return Invisibly returns \code{x}.
#' @export
print.fcoint <- function(x, ...) {
  cat("Fourier Cointegration Tests\n")
  cat(strrep("-", 60), "\n")
  cat(sprintf("Model        : %s\n", x$model))
  cat(sprintf("Criterion    : %s\n", x$criterion))
  cat(sprintf("Observations : %d\n", x$nobs))
  cat("\n")

  for (nm in names(x$results)) {
    r <- x$results[[nm]]
    if (!is.null(r$error)) {
      cat(sprintf("  [%s] Error: %s\n", toupper(nm), r$error))
      next
    }
    cat(sprintf("--- %s ---\n", toupper(nm)))
    if (nm == "tsong") {
      cat(sprintf("  CI statistic : %8.4f\n", r$CI_stat))
      cat(sprintf("  F statistic  : %8.4f\n", r$F_stat))
      cat(sprintf("  Frequency k* : %d\n",    r$frequency))
      cat(sprintf("  DOLS lags    : %d\n",    r$dolslags))
      cat(sprintf("  Obs (eff.)   : %d\n",    r$nobs))
      cat(sprintf("  CI CV 1%%/5%%/10%% : %.4f / %.4f / %.4f\n",
                  r$ci_cv1, r$ci_cv5, r$ci_cv10))
      sig <- if (!is.na(r$CI_stat) && r$CI_stat > r$ci_cv1) "Reject H0 (1%)" else
        if (!is.na(r$CI_stat) && r$CI_stat > r$ci_cv5) "Reject H0 (5%)" else
          if (!is.na(r$CI_stat) && r$CI_stat > r$ci_cv10) "Reject H0 (10%)" else
            "Fail to reject H0"
      cat(sprintf("  Decision     : %s\n", sig))
    } else {
      cat(sprintf("  t-statistic  : %8.4f\n", r$tstat))
      cat(sprintf("  delta        : %8.4f\n", r$delta))
      cat(sprintf("  Frequency k* : %d\n",    r$frequency))
      cat(sprintf("  Lag order    : %d\n",    r$lag))
      cat(sprintf("  Obs (eff.)   : %d\n",    r$nobs))
      cat(sprintf("  CV 1%%/5%%/10%% : %.4f / %.4f / %.4f\n",
                  r$cv1, r$cv5, r$cv10))
      sig <- if (!is.na(r$tstat) && r$tstat < r$cv1) "Reject H0 (1%)" else
        if (!is.na(r$tstat) && r$tstat < r$cv5) "Reject H0 (5%)" else
          if (!is.na(r$tstat) && r$tstat < r$cv10) "Reject H0 (10%)" else
            "Fail to reject H0 (unit root in residuals)"
      cat(sprintf("  Decision     : %s\n", sig))
    }
    cat("\n")
  }
  invisible(x)
}
