normalregMix
============================================

Overview
--------
normalregMix is an R package that estimates the number of components in normal regression mixture models by the method suggested in [Testing the Number of Components in Normal Mixture Regression Models](http://dx.doi.org/10.1080/01621459.2014.986272) by Kasahara and Shimotsu (2015). 

This is the main feature of this package that distinguishes normalregMix from other existing R packages available from [CRAN](http://CRAN.R-project.org/), such as [mixtools](https://cran.r-project.org/web/packages/mixtools/index.html), [mclust](https://cran.r-project.org/web/packages/mclust/index.html), and [flexmix](https://cran.r-project.org/web/packages/flexmix/index.html), as these existing packages do not provide tools to determine the number of components in normal regression mixture models except for model selection procedures that rely on information criteria (e.g., AIC and BIC).

Installing normalregMix
------------

Open R console and run the following script.
``` r
# install.packages("devtools") # if devtools has not been installed, run this line.
devtools::install_github("hkasahara/normalregmix")
```

Quick Tour
-----
The following code tests the null hypothesis that `eruptions` is a homogeneous normal random variable:
``` r
data(faithful)
attach(faithful)
library(normalregMix)
normalmixMEMtest(y = eruptions, m = 1, crit.method = 'boot')
```
Similarly, the following code can be used to test the null hypothesis that `eruptions` is a two component mixture normal random variable.
``` r
normalmixMEMtest(y = eruptions, m = 2, crit.method = 'boot')
```
Running this code repeatedly for different values of `m` can be cumbersome. The following code can be run instead to test the null hypothesis of `m= 1`, `m=2`, `m=3`, and `m=4`:
``` r
normalmixMEMtestSeq(y = eruptions, maxm = 4)
```
which also gives model selection information based on AIC and BIC as well.

Learning normalregMix
-----
Detailed tutorials can be found in the vignette in our [github repository](https://github.com/hkasahar/normalregMix/tree/master/vignettes).