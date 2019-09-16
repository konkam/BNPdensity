
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BNPdensity

Bayesian nonparametric density estimation modeling mixtures by a
Ferguson-Klass type algorithm for posterior normalized random measures.

## Installation

You can install BNPdensity from github with:

``` r
# install.packages("devtools")
devtools::install_github("konkam/BNPdensity")
```

You will need to have the CRAN package `devtools`
installed.

## How to select the parameters of the Normalized Generalized Gamma process

We suggest the Normalized stable process, which corresponds to setting
`Alpha = 1, Kappa = 0` in the `MixNRMIx` functions. The stable process
is a convenient model because its parameter γ has a convenient
interpretation: it can be used to tune how informative the prior on the
number of components is. Small values of `Gama` bring the process closer
to a Dirichlet process, where the prior on the number of components is a
relatively peaked distribution around \(\alpha \log n\). Larger values
of `Gama` make this distribution flatter. More guidelines on how to
choose the parameters may be found in Lijoi et al. (2007b), notably by
considering the expected prior number of components.

We provide a function to compute the expected number of components for a
normalized stable process:

``` r
library(BNPdensity)
expected_number_of_components_stable(100, 0.8)
#> 1 'mpfr' number of precision  5300   bits 
#> [1] 42.7094763347433064145139642374761953935253068800586645704061217932416521673216907205393634497387972550423672434248859853361487805996952132428560926301735979533962638142919404245491320418681146784117675560180180997897282151292943142921297120817286344277056532757201595800714590399801677580702348825476449555104774274897690717000458400780607335224957505642484513472533071356437832922468218811058624371598679991922247578979527622066152240312741621582983449454587742014738374261209928436865458013868459635836514415533408847942283699058466218465661304372077549346754754275121921814735346212725066178828226717075277232621519641895705391136625948612074142339580470342092448293688243480992850978992403497102847376421266264478072691640608266093779281787084116016688298257625264151386576179078300622956769752677608717950109926087638199943458107790372665551885265327164919674565823322119685775409502437206065096591074028037310013746802805763272997327361578989028373338101309163355379821181170147389233251989072209330180918520528393655319318639895074310114218241697995900064982050710406525820477683889052900277360670969502511993685089281685384529730226668468252315463797567163871977846036005950883423952365091248314552153703995699839012485091395887925071273380982174320743129811746851032817149541351629118801652203235928207673097109874143882905678985903497701028455475882246459109961351322459122178215039563311710681870538381442385456732732345958405084458427383775726002258077957996564219300933845643215966461641784547945616409630066291778349690565899480469847731894606043486473685568118112287972667409028307
```

This number may be compared to the prior number of components induced by
a Dirichlet process with `Alpha = 1`:

``` r
expected_number_of_components_Dirichlet(100, 1.)
#> [1] 5.187378
```

We also provide a way to visualise the prior distribution on the number
of components:

``` r
plot_prior_number_of_components(50, 0.4) 
#> Computing the prior probability on the number of clusters for the Dirichlet process
#> Computing the prior probability on the number of clusters for the Stable process
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- --> \#\# How to
fit a dataset

We illustrate the package by estimating the distribution of the
`acidity` dataset.

``` r
library(BNPdensity)
data(acidity)
str(acidity)
#>  num [1:155] 2.93 3.91 3.73 3.69 3.82 ...
```

``` r
hist(acidity)
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
library(BNPdensity)
data(acidity)
fit = MixNRMI1(acidity, Nit = 3000)
#> MCMC iteration 500 of 3000 
#> MCMC iteration 1000 of 3000 
#> MCMC iteration 1500 of 3000 
#> MCMC iteration 2000 of 3000 
#> MCMC iteration 2500 of 3000 
#> MCMC iteration 3000 of 3000 
#>  >>> Total processing time (sec.):
#>    user  system elapsed 
#>  98.980   0.082  99.081
```

`MixNRMI1()` creates an object of class `MixNRMI1`, for which we provide
common S3 methods.

``` r
print(fit)
#> Fit of a semiparametric normal mixture model on 155 data points.
#> The MCMC algorithm was run for 1500 iterations with 10 % discarded for burn-in.
```

``` r
plot(fit)
```

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
summary(fit)
#> Density estimation using a Normalized stable process with stability parameter Gamma = 0.4
#> A semiparametric normal mixture model was used.
#> There were 155 data points.
#> The MCMC algorithm was run for 1500 iterations with 10% discarded for burn-in.
#> To obtain information on the estimated number of clusters, please use summary(object, number_of_clusters = TRUE).
```

## How to use the convergence diagnostics

We also provide an interface to run several chains in parallel, using
the functions `multMixNRMI1()`. We interface our package with the `coda`
package by providing a conversion method for the output this function.
This allows for instance to compute the convergence diagnostics included
in `coda`.

One detail is that due to the Nonparametric nature of the model, the
number of parameters which could potentially be monitored for
convergence of the chains varies. The location parameter of the
clusters, for instance, vary at each iteration, and even the labels of
the clusters vary, which makes them tricky to follow. However, it is
possible to monitor the log-likelihood of the data along the iterations,
the value of the latent variable `u`, the number of components and for
the semi-parametric model, the value of the common scale parameter.

``` r
library(BNPdensity)
library(coda)
data(acidity)
fitlist = multMixNRMI1(acidity, Nit = 5000)
mcmc_list = as.mcmc(fitlist)
coda::traceplot(mcmc_list)
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-11-4.png)<!-- -->

``` r
coda::gelman.diag(mcmc_list)
#> Potential scale reduction factors:
#> 
#>                 Point est. Upper C.I.
#> ncomp                 1.06       1.17
#> Sigma                 1.14       1.37
#> Latent_variable       1.06       1.14
#> log_likelihood        1.09       1.26
#> 
#> Multivariate psrf
#> 
#> 1.12
```

## How to use the Goodness of fit plots

### Non censored data

``` r
library(BNPdensity)
data(acidity)
fit = MixNRMI1(acidity, extras = TRUE)
#> MCMC iteration 500 of 1500 
#> MCMC iteration 1000 of 1500 
#> MCMC iteration 1500 of 1500 
#>  >>> Total processing time (sec.):
#>    user  system elapsed 
#>  40.746   0.035  40.783
GOFplots(fit)
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

### Censored data

``` r
library(BNPdensity)
data(salinity)
fit = MixNRMI1cens(salinity$left,salinity$right, extras = TRUE)
#> MCMC iteration 500 of 1500 
#> MCMC iteration 1000 of 1500 
#> MCMC iteration 1500 of 1500 
#>  >>> Total processing time (sec.):
#>    user  system elapsed 
#>  61.339   0.013  61.358
GOFplots(fit)
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

## Posterior analysis of the clustering structure

The MCMC algorithm provides a sample of the posterior distribution on
the space of all clusterings. This is a very large discrete space, which
is not ordered. This means that for any reasonably sized problem, each
configuration in the posterior will have been explored no more than once
or twice, and that many potentially good configurations will not be
present in the MCMC sample. Moreover, the lack of ordering makes it not
trivial to summarize the posterior by an optimal clustering and to
provide credible sets.

We suggest using the approach developed in S. Wade and Z. Ghahramani,
“Bayesian cluster analysis: Point estimation and credible balls (with
discussion),” Bayesian Anal., vol. 13, no. 2, pp. 559–626, 2018.

The main proposal from this paper is to summarize the posterior on all
possible clusterings by an optimal clustering where optimality is
defined as minimizing the posterior expectation of a specific loss
function, the Variation of Information. Credible sets are also
available.

We use the implementation described in R. Rastelli and N. Friel,
“Optimal Bayesian estimators for latent variable cluster models,”
Stat. Comput., vol. 28, no. 6, pp. 1169–1186, Nov. 2018, which is faster
and implemented in the CRAN package `GreedyEPL`.

Using this approach requires installing the R package `GreedyEPL`, which
can be achieved with the following command:

``` r
install.packages("GreedyEPL")
```

Note that investigating the clustering makes more sense for the fully
Nonparametric NRMI model than for the Semiparametric. This is because to
use a single scale parameters for all the clusters, the Semiparametric
model may favor numerous small clusters, for flexibility. The larger
number of clusters may render interpretation of the clusters more
challenging.

The clustering structure may be visualized as follows:

``` r
data(acidity)
out <- MixNRMI2(acidity,  extras = TRUE)
#> MCMC iteration 500 of 1500 
#> MCMC iteration 1000 of 1500 
#> MCMC iteration 1500 of 1500 
#>  >>> Total processing time (sec.):
#>    user  system elapsed 
#>  22.760   0.044  22.804
clustering = compute_optimal_clustering(out)
plot_clustering_and_CDF(out, clustering)
```

![](README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->
