
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ERPM

<!-- badges: start -->
<!-- badges: end -->

ERPM is the software implementation of the statistical model outlined
in:

Hoffman, M., Block, P., & Snijders, T. A. (2023). Modeling partitions of
individuals. Sociological Methodology, 53(1), 1-41.

[Link to the study in Sociological
Methodology](https://doi.org/10.1177/00811750221145166).

The model allows for the analysis of emergent group compositions in
partitions, i.e., sets of non-overlapping groups. For a given partition
of individuals (or nodes, to follow the language of network analysis),
we can use this model to understand group formation processes that led
to the observation of this partition, based on individual attributes,
relations between individuals, and size-related factors. It can be seen
as an extension of Exponential Random Graph Models (ERGMs) to the case
of partition objects. With this package, one can either simulate this
model or estimate its parameters for a given dataset.

This package also provides a longitudinal extension of the model, for a
list of partitions, where each partition depends on the previous
partitions. It follows the definition proposed in:

Hoffman, M., & Chabot, T. (2023). The role of selection in socioeconomic
homophily: Evidence from an adolescent summer camp. Social Networks, 74,
259-274.

The **ERPM Manual** is available here on github in the documentation
folder.

# Note from the developers

The package and the documentation might still have bugs or errors, or
you might not be able to do what you want. In that case, or if you are
unsure please create an issue here or directly send an email to the
package maintainer (marion.hoffman\[at\]iast.fr).

We are currently (Mar 2024) on version 0.1.0 on github.

# Installation

You can install ERPM either from GitHub or from CRAN.

``` r
# from GitHub:
# install.packages("remotes")
remotes::install_github("stocnet/ERPM")
```

``` r
# from CRAN:
install.packages("ERPM")
```

# Cross-sectional example

In this section, we outline a simple example with synthetic data.

``` r
library(ERPM)
```

## The Data

Let us define a set of n = 6 nodes with three attributes (label, gender,
and age), and an arbitrary covariate matrix (friendship). We create the
following dataframe.

``` r
n <- 6 
nodes <- data.frame(label = c("A","B","C","D","E","F"),
                    gender = c(1,1,2,1,2,2),
                    age = c(20,22,25,30,30,31)) 
friendship <- matrix(c(0, 1, 1, 1, 0, 0,
                       1, 0, 0, 0, 1, 0,
                       1, 0, 0, 0, 1, 0,
                       1, 0, 0, 0, 0, 0,
                       0, 1, 1, 0, 0, 1,
                       0, 0, 0, 0, 1, 0), 6, 6, TRUE) 
```

We consider a partition for these 6 individuals. We define a vector with
six elements, indicating the id of each individual’s group.

``` r
partition <- c(1,1,2,2,2,3)
```

## Model specification

First, we need to choose the effects (i.e., explaining variables) we
want to include. For example we set four (which is of course not
reasonable for 6 nodes): 1. “num_groups”: tendency to form more groups
2. “same” (for the gender covariate): tendency to form groups with
individuals with the same gender 3. “diff” (for the age covariate):
tendency to form groups with individuals with high age differences 4.
“tie” (for the friendship covariate): tendency to form groups with
friends

``` r
effects <- list(names = c("num_groups","same","diff","tie"),
                objects = c("partition","gender","age","friendship"))
objects <- list()
objects[[1]] <- list(name = "friendship", object = friendship)
```

The effect objects should contain names of pre-written functions in the
package (see manual for all effect names) as well as objects they are
referring to (either the partition, or covariates). When objects are not
individual covariates, we need to create an additional list to store
these extra objects.

## Estimation

The parameters of the model can be estimated using Maximum-likelihood
estimation (equivalent to the Method of Moment estimation). For more
details on the parametrization of the estimation algorithm, see the
manual.

``` r
estimation <- estimate_ERPM(partition, 
                          nodes, 
                          objects, 
                          effects, 
                          startingestimates = c(-1.5,0.2,-0.2,0.2), 
                          burnin = 100, 
                          thining = 20,
                          length.p1 = 500, # number of samples in phase 1
                          multiplication.iter.p2 = 20,  # multiplication factor for the number of iteration in phase 2 subphases 
                          num.steps.p2 = 4, # number of phase 2 subphases
                          length.p3 = 1000) # number of samples in phase 3
#> [1] "Observed statistics"
#> [1]  3  2 12  2
#> [1] "Burn-in"
#> [1] 100
#> [1] "Thining"
#> [1] 20
#> [1] "Autocorrelations in phase 1:"
#> [1] 0.11580435 0.07567835 0.13095529 0.10673909
#> [1] "Covariance matrix"
#>            [,1]       [,2]      [,3]       [,4]
#> [1,]  0.4893146 -0.3568016 -2.100393 -0.4514068
#> [2,] -0.3568016  0.9933908  2.641619  0.6344008
#> [3,] -2.1003928  2.6416192 41.244313  2.1706293
#> [4,] -0.4514068  0.6344008  2.170629  0.8254068
#> [1] "Invert scaling matrix"
#>            [,1]        [,2]         [,3]         [,4]
#> [1,] 3.26174503  0.18897837  0.074621129  1.153866020
#> [2,] 0.18897837  1.54410755 -0.038320247 -0.786131520
#> [3,] 0.07462113 -0.03832025  0.029495348 -0.005843019
#> [4,] 1.15386602 -0.78613152 -0.005843019  2.212018064
#> [1] "Estimated statistics after phase 1"
#> [1]  3.108  1.686 10.224  1.838
#> [1] "Estimates after phase 1"
#>            [,1]
#> [1,] -1.4734758
#> [2,]  0.4690300
#> [3,] -0.1686545
#> [4,]  0.1765069
#> Difference to estimated statistics after phase 2, step 1NULL
#> [1]  0.04166667 -0.04166667  0.58333333  0.04166667
#> Estimates after phase 2, step 1NULL
#>                
#> [1,] -1.8338234
#> [2,]  0.7115367
#> [3,] -0.2248016
#> [4,] -0.2998742
#> Difference to estimated statistics after phase 2, step 2NULL
#> [1]  0.00000000  0.01960784  0.43137255 -0.05882353
#> Estimates after phase 2, step 2NULL
#>                
#> [1,] -1.9462692
#> [2,]  0.6089256
#> [3,] -0.2275217
#> [4,] -0.1412861
#> Difference to estimated statistics after phase 2, step 3NULL
#> [1]  0.01550388 -0.02325581 -0.41860465 -0.02325581
#> Estimates after phase 2, step 3NULL
#>                
#> [1,] -1.9945672
#> [2,]  0.6111215
#> [3,] -0.2188468
#> [4,] -0.1050062
#> Difference to estimated statistics after phase 2, step 4NULL
#> [1] -0.01869159  0.04672897  0.05919003  0.03115265
#> Estimates after phase 2, step 4NULL
#>                 
#> [1,] -1.84065662
#> [2,]  0.53564516
#> [3,] -0.20753277
#> [4,] -0.09767082
#> [1] "Autocorrelations in phase 3:"
#> [1] 0.03606953 0.03474974 0.01362077 0.04370177
#> [1] "Estimated statistics after phase 3"
#> [1]  3.031  1.955 11.167  1.851
#> [1] "Estimates after phase 3"
#>                 
#> [1,] -1.84065662
#> [2,]  0.53564516
#> [3,] -0.20753277
#> [4,] -0.09767082
#>       effect     object         est   std.err sig          t        conv
#> 1 num_groups  partition -1.84065662 2.1976543     -0.8375551  0.04509781
#> 2       same     gender  0.53564516 1.4070600      0.3806840 -0.03910384
#> 3       diff        age -0.20753277 0.1715036     -1.2100785 -0.11187049
#> 4        tie friendship -0.09767082 1.9009242     -0.0513807 -0.16260731
estimation$results
#>       effect     object         est   std.err        conv
#> 1 num_groups  partition -1.84065662 2.1976543  0.04509781
#> 2       same     gender  0.53564516 1.4070600 -0.03910384
#> 3       diff        age -0.20753277 0.1715036 -0.11187049
#> 4        tie friendship -0.09767082 1.9009242 -0.16260731
```

### Simulation

We can check how the model reproduces statistics of the observed data by
simulating the estimated model. We can also use this function to
simulate theoretical models.

``` r
nsimulations <- 1000
simulations <- draw_Metropolis_single(theta = estimation$results$est, 
                          first.partition = partition, 
                          nodes = nodes, 
                          effects = effects, 
                          objects = objects, 
                          burnin = 100, 
                          thining = 20, 
                          num.steps = nsimulations, 
                          neighborhood = c(1,1,1), 
                          sizes.allowed = 1:n,
                          sizes.simulated = 1:n,
                          return.all.partitions = T)
```

### Log-likelihood and AIC

Finally, we can estimate the log-likelihood and AIC of the model (useful
to compare two models for example). First we need to estimate the ML
estimates of a simple model with only one parameter for number of groups
(this parameter should be in the model!).

``` r
likelihood_function <- function(x){ exp(x*max(partition)) / compute_numgroups_denominator(n,x)}
curve(likelihood_function, from=-2, to=0)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

``` r
parameter_base <- optimize(likelihood_function, interval=c(-2, 0), maximum=TRUE)
parameters_basemodel <- c(parameter_base$maximum,0,0,0)
```

Then we can get our estimated logL and AIC.

``` r
logL_AIC <- estimate_logL(partition,
                         nodes,
                         effects, 
                         objects,
                         theta = estimation$results$est,
                         theta_0 = parameters_basemodel,
                         M = 3,
                         num.steps = 200,
                         burnin = 100,
                         thining = 20)
#> [1] "step 1"
#> [1] "step 2"
#> [1] "step 3"
logL_AIC$logL
#>           [,1]
#> [1,] -4.380505
logL_AIC$AIC
#>           [,1]
#> [1,] 0.7610098
```

# More …

For more details on the longitudinal version of the model or other
functions, have a look at the manual in the documentation folder or the
example script in the scripts folder.
