# BNPdensity

A fork of the BNPdensity package from CRAN (https://cran.r-project.org/web/packages/BNPdensity/index.html)


## How to use the Goodness of fit plots

### Non censored data
```
library(BNPdensity)
data(acidity)
fit = MixNRMI1(acidity, extras = TRUE)
plotGOF(fit)
```

![](GOFplot_noncensored.png)

### Censored data
```
library(BNPdensity)
data(salinity)
fit = MixNRMI1cens(salinity$left,salinity$right, extras = TRUE)
plotGOF(fit)
```

![](GOFplot_censored.png)

