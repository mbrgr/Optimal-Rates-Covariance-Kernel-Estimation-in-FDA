
# Optimal Rates Covariance Kernel Estimation in FDA

<!-- badges: start -->
<!-- badges: end -->

This repository contains the code for the Figures of the paper “Optimal
rates for estimating the covariance kernel from synchronously sampled
functional data” by Hajo Holzmann and Max Berger. Note the this code
needs the package “biLocPol” of which you can get the development
version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mbrgr/biLocPol")
# library(biLocPol)
```

## Structure

The R-File “functions.R” contains additional functions for the bandwidth
comparison, error decomposition and simulation of the cross validation
procedure for the simulations that are not contained in the “biLocPol”
package. No data is produced here. The file “illustrations.R” contains
the code for Figures 1, 2 and 7 as well as additional grafics with the
processes used. All the code that is produced by this file is saved in
the file “illustrations.Rdata” in the “data” folder. The file
“bw_cv_error_decomposition.r” contains the code for most of the repeated
simulations as the bandwidth comparison, the cross validation and the
error decomposition. The results are stored in different .Rdata files in
the “data” folder. The evaluation of these simulations can be found in
the “simulation_evaluation.R” file (Figure 3, 4, 5, 6 and 8). The
comparison of our proposed estimator (only estimating on one side of the
diagonal and mirroring) compared to a vanilla bivariate local polynomial
estimator can be found in “estimator_comparison.R”. Again the evaluation
of the results can be found in “simulation_evaluation.R” and the “data”
folder. Finally the Figures 9, 10, 11 of the real data example are
generated with the code in “temperature_analysis.R” along with some
additional grafics.

## Comments

Some of the “.Rdata” files are large since they contain the actual
weights of some of the estimations. Note that the calculation of the
weights and the simulations are paralallized with the “future.apply” and
“future” package.
