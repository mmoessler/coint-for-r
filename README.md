# coint-for-r

A suite of R functions for nonstationary time series and model selection.

This is a project to implement some of the procedures provided by Peter C.B. Phillips and Sam Ouliaris in R.

There are already some where good R-packages for the analysis of nonstationary and cointegrated time series such as [urca](https://cran.r-project.org/web/packages/urca/index.html) and [cointReg](https://cran.r-project.org/web/packages/cointReg/index.html) available on [CRAN](https://cran.r-project.org/).

However, there were still some specifications and procedures missing that were useful for my own work and that I would like to share here.

## FM_OLS

Function which computes the Phillips (1995) "Fully Modified" OLS estimator for single and multi equation cointegrated regression models.

```
fm_ols(y, x, d, l)
```

Inputs:

* `y`: Dependent variables $\left(nobs \times n\right)$
* `x`: Explanatory variables $\left(nobs \times m\right)$
* `d`: Deterministic components $\left(nobs \times m_d\right)$
* `l`: Number of autocovariance terms to compute the spectrum at frequency zero, i.e., the long-run variance

Outputs:

* `beta`: $\left(m + m_d\right)n$ vector containing the parameter estimates
* `vc`: Variance matrix for the parameter estimates

## References

Phillips, P. C. (1995). Fully modified least squares and vector autoregression. *Econometrica: Journal of the Econometric Society*, 1023-1078.
