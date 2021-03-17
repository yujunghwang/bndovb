bndovb
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

## Introduction

The R package **bndovb** implements a Hwang(2021) estimator to bound
omitted variable bias using auxiliary data. The basic assumption is that
the main data includes a dependent variable and every regressor but one
omitted variable. So there is an omitted variable bias in the OLS result
from the main data. However, if there is auxiliary data that includes
every right-hand side regressor (or its noisy proxies), it can bound the
omitted variable bias. Hwang(2021) provides a more general estimator for
when there is more than one omitted variable and when the auxiliary data
does not contain every right-hand side regressor (or its noisy proxies).

This package implements a simple estimator when the number of omitted
variables is just one, and the auxiliary data contains every right-hand
side regressor but the only omitted variable. This package provides two
different functions.

The function ‘bndovb’ can be used when the auxiliary data contain every
right-hand side variable without measurement errors. The function
‘bndovbme’ can be used when noisy proxies for the omitted variable
exist in the auxiliary data. Other regressors in the auxiliary data are
assumed measurement error-free. The function requires another R package,
‘factormodel,’ written by the same author.

When using ‘bndovb,’ a user should specify a method for density
estimation as part of an estimation procedure, either 1 (parametric
normal density assumption) or 2 (nonparametric kernel density
estimation). In general, it is strongly recommended to use method 1
(parametric normal density assumption), particularly when data is large
or the regression model is large. Method 2 calls an R package “np” (Li
and Racine, 2008; Li, Lin and Racine, 2013) but this method is very
slow, as emphasized in their vignette file. I recommend using their
method only when (i) data is small and (ii) there is only one common
variable, which makes the conditional CDF and quantile function
univariate. The default method is 1, using parametric normal density
assumption.

## Installation

You can install a package **bndovb** using either CRAN or github.

``` r
install.packages("bndovb")
```

or

``` r
# install.packages("devtools")
devtools::install_github("yujunghwang/bndovb")
```

## Example 1 : bndovb

The below example shows how to use a function ‘bndovb’ when the
auxiliary data contain the omitted variable from main data without any
measurement error. The code first simulates fake data using the same DGP
in both main data and auxiliary data. Next, the main data omits one
variable. The function ‘bndovb’ provides a bound on regression
coefficients by using both main data and auxiliary data.

``` r
library(bndovb)
#> Warning: replacing previous import 'MASS::select' by 'dplyr::select' when
#> loading 'bndovb'
#> Warning: replacing previous import 'dplyr::filter' by 'stats::filter' when
#> loading 'bndovb'
#> Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading
#> 'bndovb'
library(MASS)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following object is masked from 'package:MASS':
#> 
#>     select
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

# sample size
Nm <- 100000 # main data
Na <- 50000 # auxiliary data

# use same DGP in maindat and auxdat
maindat <- mvrnorm(Nm,mu=c(2,3,1),Sigma=rbind(c(2,1,1),c(1,2,1),c(1,1,2)))
auxdat  <- mvrnorm(Na,mu=c(2,3,1),Sigma=rbind(c(2,1,1),c(1,2,1),c(1,1,2)))

maindat <- as.data.frame(maindat)
auxdat <- as.data.frame(auxdat)

colnames(maindat) <- c("x1","x2","x3")
colnames(auxdat) <- c("x1","x2","x3")

# this is a true parameter which we try to get bounds on
truebeta <- matrix(c(2,1,3,2),ncol=1)

# generate a dependent variable
maindat$y <- as.matrix(cbind(rep(1,Nm),maindat[,c("x1","x2","x3")]))%*%truebeta

# main data misses one omitted variable "x1"
maindat <- maindat[,c("x2","x3","y")]

# use "bndovb" function assuming parametric "normal" distribution for the CDF and quantile function (set method=1)
# see Hwang(2021) for further details
oout <- bndovb(maindat=maindat,auxdat=auxdat,depvar="y",ovar="x1",comvar=c("x2","x3"),method=1)
print(oout)
#> $hat_beta_l
#>        con         x1         x2         x3 
#>  1.9261145 -0.9231628  2.9851454  2.0125252 
#> 
#> $hat_beta_u
#>      con       x1       x2       x3 
#> 3.280312 1.059706 3.640417 2.664897

# use "bndovb" function using nonparametric estimation of the CDF and quantile function (set method=2)
# for nonparametric density estimator, the R package "np" was used. See Hayfield and Racine (2008), Li and Racine (2008), Li, Lin and Racine (2013)
#### The next line takes very long because of large sample size. You can try using a smaller sample and run the next line.
#oout <- bndovb(maindat=maindat,auxdat=auxdat,depvar="y",ovar="x1",comvar=c("x2","x3"),method=2)
#print(oout)
```

## Conclusion

This vignette showed how to use functions in \`bndovb’ R package.

## References

  - [Hayfield, Tristen, and Jeffrey S. Racine. “Nonparametric
    econometrics: The np package.” Journal of statistical software 27,
    no. 5 (2008): 1-32.](https://doi.org/10.18637/jss.v027.i05)

  - [Hwang, Yujung (2021). Bounding Omitted Variable Bias Using
    Auxiliary Data. Working
    Paper.](https://sites.google.com/view/yujunghwang/research?authuser=0)

  - [Li, Qi, and Jeffrey S. Racine. “Nonparametric estimation of
    conditional CDF and quantile functions with mixed categorical and
    continuous data.” Journal of Business & Economic Statistics 26,
    no. 4
    (2008): 423-434.](https://www.tandfonline.com/doi/abs/10.1198/073500107000000250?casa_token=tUdlEvmt_fMAAAAA:IALIN23XRN8bpnL0ZBEQp8MIUO9Ie_SkHzdBTBVa-DD1oYaxdqaqrr_EyBD7IMKjxmXIanGfHIBW)

  - [Li, Qi, Juan Lin, and Jeffrey S. Racine. “Optimal bandwidth
    selection for nonparametric conditional distribution and quantile
    functions.” Journal of Business & Economic Statistics 31, no. 1
    (2013): 57-65.](https://www.tandfonline.com/doi/full/10.1080/07350015.2012.738955?casa_token=8S1iR7ki1qkAAAAA%3AQQI6lwjAKgYfv5dmpbvmbCZckpscxYXFUSvZlXQ64Gz8D45E1yYEh0BPQF5DSg0chfTkcvG6HMim)
