\name{fcoint}
\alias{fcoint}
\title{Fourier Cointegration Tests for Time Series}
\usage{
fcoint(y, x, test = c("fadl", "feg", "feg2", "tsong", "all"),
  model = c("constant", "trend"), max_freq = 5L, max_lag = 0L,
  criterion = c("aic", "bic"), dols_lags = 0L)
}
\arguments{
\item{y}{Numeric vector. Dependent (left-hand side) time-series variable.}

\item{x}{Numeric matrix or vector. Independent (right-hand side) variables.}

\item{test}{Character. Which test to run: \code{"fadl"}, \code{"feg"},
  \code{"feg2"}, \code{"tsong"}, or \code{"all"}. Default is \code{"fadl"}.}

\item{model}{Character. Deterministic component: \code{"constant"} (default)
  or \code{"trend"}.}

\item{max_freq}{Integer. Maximum Fourier frequency to search over (1--5).
  Default is \code{5}.}

\item{max_lag}{Integer. Maximum lag order. \code{0} selects automatically
  via the information criterion. Default is \code{0}.}

\item{criterion}{Character. Information criterion for lag selection:
  \code{"aic"} (default) or \code{"bic"}.}

\item{dols_lags}{Integer. Leads and lags for the DOLS estimator used in
  the Tsong test. Default is \code{0} (automatic).}
}
\value{
A list of class \code{"fcoint"} containing:
\describe{
  \item{test}{Character. Name(s) of the test(s) performed.}
  \item{results}{Named list of individual test results.}
  \item{model}{Character. Deterministic component used.}
  \item{criterion}{Character. Information criterion used.}
  \item{nobs}{Integer. Effective number of observations.}
}
}
\description{
Implements four Fourier-based cointegration tests that accommodate smooth
structural breaks via flexible Fourier terms: FADL, FEG, FEG2, and Tsong.
}
\details{
The Fourier approach allows for smooth, nonlinear breaks in the deterministic
components by augmenting standard cointegration regressions with one or more
pairs of sine and cosine terms.

\strong{FADL (Banerjee, Arcabic & Lee, 2017):}
An ADL-type residual-based test. The optimal Fourier frequency \eqn{k^*}
and lag order are chosen jointly by minimising the selected information
criterion over a grid.

\strong{FEG and FEG2 (Banerjee & Lee):}
Engle-Granger style residual-based tests augmented with Fourier terms.
FEG2 includes an additional \eqn{R^2} correction.

\strong{Tsong et al. (2016):}
A DOLS-based cointegration test with Fourier terms. Reports both a
CI statistic and an \eqn{F}-statistic for joint significance of Fourier
terms.
}
\examples{
set.seed(42)
n <- 80
x <- cumsum(rnorm(n))
y <- 0.5 * x + rnorm(n, sd = 0.3)
res <- fcoint(y, x, test = "fadl", max_freq = 3)
print(res)
}
\references{
Banerjee, P., Arcabic, V., & Lee, H. (2017). Fourier ADL cointegration test
to approximate smooth breaks with new evidence from crude oil market.
\emph{Economic Modelling}, 67, 114--124.
\doi{10.1016/j.econmod.2017.03.004}

Tsong, C.-C., Lee, C.-F., Tsai, L.-J., & Hu, T.-C. (2016). The Fourier
approximation and testing for the null of cointegration.
\emph{Empirical Economics}, 51(3), 1085--1113.
\doi{10.1007/s00181-015-0921-5}
}
