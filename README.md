# Peer effects in the linear-in-means model may be inestimable even when identified

Linear-in-means models are often used to investigate peer effects. Estimating peer effects in linear-in-means models requires care, as peer effects may be subject to the ``reflection problem'', an identification failure in the form of perfect collinearity. In many settings, well-known identification conditions guarantee that perfect collinearity is not an issue. However, these identifying conditions are not sufficient to guarantee that peer effects are estimable. Even when identifying conditions guarantee that peer effect terms are not collinear, peer effects can become increasingly collinear as sample size grows larger. We show that asymptotic collinearity occurs whenever nodal covariates are independent of the network and the minimum degree of the network is growing. Asymptotic collinearity can cause estimates of peer effects to be inconsistent or to converge at slower than expected rates. We also show that dependence between nodal covariates and network structure can alleviate collinearity issues in random dot product graphs. These results suggest that linear-in-means models are less reliable for studying peer influence than previously believed.

## To replicate our computational results

We use [`renv`](https://rstudio.github.io/renv/) to record package dependencies and [`targets`](https://books.ropensci.org/targets/) to coordinate our simulation study and data analysis.

To replicate our results, clone this Github repository. Once you have the repository cloned locally, re-create the project library by calling

``` r
# install.packages("renv")
renv::restore()
```

At this point, you should be ready to replicate our simulation results.

``` r
tar_make()
```

This is a somewhat computationally intensive project. It's likely that, at some point, the targets pipeline will crash for some sundry computational reason. If this is the case, simply re-run `tar_make()`. `tar_make()` will only attempt to re-run incomplete portions of the build pipeline.

Results will appear as image files in the `figures/` folder.
