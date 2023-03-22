# coint-for-r

A suite of R functions for nonstationary time series and model selection.

This is a project to implement some of the procedures provided by Peter C.B. Phillips and Sam Ouliaris in R.

There are already some where good R-packages for the analysis of nonstationary and cointegrated time series such as [urca](https://cran.r-project.org/web/packages/urca/index.html) and [cointReg](https://cran.r-project.org/web/packages/cointReg/index.html) available on [CRAN](https://cran.r-project.org/).

However, there were still some specifications and procedures missing that were useful for my own work and that I would like to share here.

## FM_OLS

Function which computes the Phillips (1995) "Fully Modified" OLS estimator for single and multi equation cointegrated regression models.
