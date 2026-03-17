# fcoint

**Fourier Cointegration Tests for Time Series with Smooth Structural Breaks**

[![CRAN status](https://www.r-pkg.org/badges/version/fcoint)](https://CRAN.R-project.org/package=fcoint)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

`fcoint` implements four Fourier-based cointegration tests for time series that accommodate smooth structural breaks via flexible trigonometric (Fourier) terms:

| Test | Reference |
|------|-----------|
| **FADL** — Fourier ADL | Banerjee, Arcabic & Lee (2017) |
| **FEG** — Fourier Engle-Granger | Banerjee & Lee |
| **FEG2** — FEG with R² correction | Banerjee & Lee |
| **Tsong** — DOLS-based | Tsong, Lee, Tsai & Hu (2016) |

## Installation

```r
install.packages("fcoint")
```

Or from GitHub:

```r
# install.packages("remotes")
remotes::install_github("muhammedalkhalaf/fcoint")
```

## Quick Start

```r
library(fcoint)

set.seed(42)
n <- 100
x <- cumsum(rnorm(n))
y <- 0.5 * x + rnorm(n, sd = 0.3)

# Run FADL test
res <- fcoint(y, x, test = "fadl", max_freq = 3)
print(res)

# Run all tests
res_all <- fcoint(y, x, test = "all")
print(res_all)
```

## References

Banerjee, P., Arcabic, V., & Lee, H. (2017). Fourier ADL cointegration test to approximate smooth breaks with new evidence from crude oil market. *Economic Modelling*, 67, 114–124. <https://doi.org/10.1016/j.econmod.2017.03.004>

Tsong, C.-C., Lee, C.-F., Tsai, L.-J., & Hu, T.-C. (2016). The Fourier approximation and testing for the null of cointegration. *Empirical Economics*, 51(3), 1085–1113. <https://doi.org/10.1007/s00181-015-0921-5>

## Author

Muhammad Alkhalaf <muhammedalkhalaf@gmail.com>
