\name{print.fcoint}
\alias{print.fcoint}
\title{Print Method for fcoint Objects}
\usage{
\method{print}{fcoint}(x, ...)
}
\arguments{
\item{x}{An object of class \code{"fcoint"}.}

\item{...}{Further arguments passed to or from other methods (unused).}
}
\value{
Invisibly returns \code{x}.
}
\description{
Prints a formatted summary of Fourier cointegration test results.
}
\examples{
set.seed(1)
n <- 60
x <- cumsum(rnorm(n))
y <- 0.5 * x + rnorm(n)
res <- fcoint(y, x, test = "fadl", max_freq = 2)
print(res)
}
