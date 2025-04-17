library(biLocPol) # please install this package from Github first. See "README.md" file for instructions
library(future.apply)
library(tidyverse)
library(tictoc)



##### Illustration #####

del10_cov_OU = function(t, sigma = 3, theta = 2){
  if(t[1] > t[2]){
    return(sigma^2/2 * (exp(-theta*(t[1] - t[2])) + exp(-theta*(t[1]+t[2])) ))
  }
  if(t[1] == t[2]){return(sigma^2 * exp(-theta*2*t[1]))}
  if(t[1] < t[2]){
    return(sigma^2/2 * (-exp(theta*(t[1] - t[2])) + exp(-theta*(t[1]+t[2])) ))
  }
}

del01_cov_OU = function(t, sigma = 3, theta = 2){
  return(del10_cov_OU(c(t[2], t[1]), sigma, theta))
}


x = observation_grid(p.eval, comp = "full")[1:p.eval, 1]
cov10.ou.eval = matrix(apply(observation_grid(p.eval, comp = "full"), 1, del10_cov_OU), p.eval, p.eval)
plot_ly() |> 
  add_surface(x = x, y = x, z = cov10.ou.eval, colorscale = c("blue", "yellow"), alpha= .5) |> 
  add_surface(x = x, y = x, z = eval_weights(w, Z)[,,3], alpha = .5,
              colors = c('#BF382A', '#0C4B8E')) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = ""))) 

##### Bandwidth Evaluation Function #####

bandwidth_evaluation_derivative_OU = function(h, p, p.eval, n.seq, N){
  
  x.design = (1:p - 0.5)/p
  s = length(h)
  x.eval.grid = observation_grid(p.eval, comp = "full")
  G10 = apply(x.eval.grid, 1, del10_cov_OU)
  
  sup.err10 = numeric(length(n.seq))
  w = local_polynomial_weights(p, h, p.eval, F, m = 2, grid.type = "less", del = 1)
  for(k in 1:length(n.seq)){
    n = n.seq[k]
    sup.err[k] =  mean(replicate(N, {
      Y = biLocPol::OU(n, x.design, sigma = sigma, alpha = theta) + matrix(rnorm(n * p, 0, 0.75), n, p) # n x p
      Z = observation_transformation(Y)
      trim = which(x.design > 0.1 & x.design < (1-0.1))
      estimate = eval_weights(w, Z)[,,3]
      max(abs(estimate - G10)[trim, trim])
    }))
  }
  rm(w)
  matrix(c(n.seq, rep(c(p, h), each = length(n.seq)), sup.err), ncol =  4)
}


##### Bandwidth Comparison for Covariance derivative estimation #####
n = 400
N = 1000
n.seq = c(50, 100, 200, 400)
p.seq = c(15, 25, 50, 75, 100)

p.eval = 100
H = lapply(1:length(p.seq), function(l){seq(.7, 3/p.seq[l], -0.05)})

# Parameter OU Process
theta = 2; sigma = 3
# Standard deviation for additional errors
sd = 0.75

bw_comparison = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::multisession)

for(l in 1:length(p.seq)){
  tic()
  bw_comparison[[l]] = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         cov.ou, list(theta = theta, sigma = sigma), 
                         OU, list(alpha = theta, sigma = sigma, x0 = 0), 
                         eps.arg = list(sd = sd))
  },
  future.seed = T)),ncol = 4)
  cat("p =", p.seq[l], "done.")
  toc()
}

parallel::stopCluster(cl)