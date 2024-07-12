#### Simulation for Covariance Estimation with other Process #####
# Z(x) = sqrt(2)/6 * (N^2 - 1) * sin(pi * x)  + 2/3 * (M - 1) * (x - 0.5)
library(tidyverse)
library(plotly)
library(biLocPol)

z_two_rv = function(n = 1, p = 100, t = NULL){
  if( is.null(t) ) {
    t = (1:p - 1/2)/p
  }
  N = matrix(rnorm(2*n, 0, 1), n, 2)
  apply(N, 1, function(m){
  2/3 * m[1] * sin(pi * t) + 
    sqrt(2)*2/3 * m[2] * cos(1.25*pi*t)
  }) |> t()
}


cov_z_2rv = function(t){
  4/9*sin(pi*t[1])*sin(pi*t[2]) + 
    8/9*cos(1.25*pi*t[1])*cos(1.25*pi*t[2])
}


