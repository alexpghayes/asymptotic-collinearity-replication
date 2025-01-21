# Peer effects in the linear-in-means model may be inestimable even when identified

Linear-in-means models are often used to investigate peer effects. Estimating peer effects in linear-in-means models requires care, as peer effects may be subject to the ``reflection problem'', an identification failure in the form of perfect collinearity. In many settings, well-known identification conditions guarantee that perfect collinearity is not an issue. However, these identifying conditions are not sufficient to guarantee that peer effects are estimable. Even when identifying conditions guarantee that peer effect terms are not collinear, peer effects can become increasingly collinear as sample size grows larger. We show that asymptotic collinearity occurs whenever nodal covariates are independent of the network and the minimum degree of the network is growing. Asymptotic collinearity can cause estimates of peer effects to be inconsistent or to converge at slower than expected rates. We also show that dependence between nodal covariates and network structure can alleviate collinearity issues in random dot product graphs. These results suggest that linear-in-means models are less reliable for studying peer influence than previously believed.

## TODO:

renv::snapshot() once things are working

## To replicate our computational results

We use [`renv`](https://rstudio.github.io/renv/) to record package dependencies and [`targets`](https://books.ropensci.org/targets/) to coordinate our simulation study and data analysis.

To replicate our results, clone this Github repository to your local computer. Once you have the repository cloned locally, re-create the project library by calling

``` r
# install.packages("renv")
renv::restore()
```

At this point, you should be ready to replicate our simulation results.

``` r
tar_make()
```

This is a computationally intensive project. It's likely that, at some point, the targets pipeline will crash for some sundry computational reason. If this is the case, simply re-run `tar_make()`. `tar_make()` will only attempt to re-run incomplete portions of the build pipeline.

To fully replicate the paper, increase the number of replications in the misspecification and simulation studies by setting:

- `tar_target(num_chunks, 1)` -> `tar_target(num_chunks, 10)` on line 73 of `_misspecification.R`

otherwise the misspecification and simulation studies will run using fewer replicates than in the paper (for computational purposes). We recommend you make sure the code runs with the default values before increasing them to the full replication size.

Results will appear as imagine files in the `figures/` folder.
