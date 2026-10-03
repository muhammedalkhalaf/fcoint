# fcoint

> **This repository is superseded and no longer maintained.**
> At the request of the CRAN team, this package was merged into the CRAN package
> [cointests](https://cran.r-project.org/package=cointests). The function `fcoint()` is maintained there,
> with corrections that are not in this repository. The code here is an older version
> and should not be used for new work.
>
> ```r
> install.packages("cointests")
> ```

**Fourier Cointegration Tests for Time Series with Smooth Structural Breaks**

[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

`fcoint` implements four Fourier-based cointegration tests for time series that accommodate smooth structural breaks via flexible trigonometric (Fourier) terms:

| Test | Reference |
|------|-----------|
| **FADL**: Fourier ADL | Banerjee, Arcabic and Lee (2017) |
| **FEG**: Fourier Engle-Granger | Banerjee and Lee |
| **FEG2**: FEG with R² correction | Banerjee and Lee |
| **Tsong**: DOLS-based | Tsong, Lee, Tsai and Hu (2016) |

## Installation

```r
# Old version, GitHub only (see the notice at the top of this page)
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

Banerjee, P., Arcabic, V. and Lee, H. (2017). Fourier ADL cointegration test to approximate smooth breaks with new evidence from crude oil market. *Economic Modelling*, 67, 114–124. <https://doi.org/10.1016/j.econmod.2016.11.004>

Tsong, C.-C., Lee, C.-F., Tsai, L.-J. and Hu, T.-C. (2016). The Fourier approximation and testing for the null of cointegration. *Empirical Economics*, 51(3), 1085–1113. <https://doi.org/10.1007/s00181-015-1028-6>

## Author

Muhammad Alkhalaf <muhammedalkhalaf@gmail.com>
