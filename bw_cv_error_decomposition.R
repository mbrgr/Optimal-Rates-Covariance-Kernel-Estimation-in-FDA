#### Packages ####

library(biLocPol)
library(future.apply)
library(tidyverse)
library(tictoc)

#### Bandwidth Comparison ####

N = 1000
n.seq = c(50, 100, 200, 400)
p.seq = c(15, 25, 50, 75, 100)
# p.seq = 25

p.eval = 100
H = lapply(1:length(p.seq), function(l){seq(1, 3/p.seq[l], -0.05)})
# H = lapply(1:length(p.seq), function(l){seq(.6, 3/p.seq[l], -0.2)})

# Parameter OU Process
theta = 2; sigma = 3
# Standard deviation for additional errors
sd = 0.75

bw_comparison = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

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



#### Cross Validation ####
##### One-Fold Cross Validation #####

one_fold_cv = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  one_fold_cv[[l]] =  lopocv_sim(N, n = 100, p.seq[l], H[[l]], 
               m = 1, w.parallel = T, 
               theta = theta, sigma = sigma, sd = sd)
  cat("p =", p.seq[l], "done.")
  toc()
}

parallel::stopCluster(cl)

##### 5-Fold Cross Validation #####

five_fold_cv = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  five_fold_cv[[l]] = k_fold_cv_simulation(N, n = 100, p.seq[l], H[[l]], 
               K = 5, m = 1, w.parallel = T, 
               theta = theta, sigma = sigma, sd = sd)
  cat("p =", p.seq[l], "done.")
  toc()
}

parallel::stopCluster(cl)



##### One-Fold Cross Validation #####
# different n

one_fold_cv_n400 = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:3){
  tic()
  one_fold_cv_n400[[l]] =  lopocv_sim(N, n = 400, p.seq[l], H[[l]], 
                                 m = 1, w.parallel = T, 
                                 theta = theta, sigma = sigma, sd = sd)
  cat("p =", p.seq[l], "done.")
  toc()
}

parallel::stopCluster(cl)

##### 5-Fold Cross Validation #####

five_fold_cv = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  five_fold_cv[[l]] = k_fold_cv_simulation(N, n = 100, p.seq[l], H[[l]], 
                                           K = 5, m = 1, w.parallel = T, 
                                           theta = theta, sigma = sigma, sd = sd)
  cat("p =", p.seq[l], "done.")
  toc()
}

parallel::stopCluster(cl)

##### 5-Fold Cross Validation #####

five_fold_cv_n400 = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  five_fold_cv_n400[[l]] = k_fold_cv_simulation(N, n = 400, p.seq[l], H[[l]], 
                                           K = 5, m = 1, w.parallel = T, 
                                           theta = theta, sigma = sigma, sd = sd)
  cat("p =", p.seq[l], "done. ")
  toc()
}

parallel::stopCluster(cl)

##### 2-Fold Cross Validation #####

two_fold_cv = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  two_fold_cv[[l]] = k_fold_cv_simulation(N, n, p.seq[l], H[[l]], 
                         K = 5, m = 1, w.parallel = T, 
                         theta = theta, sigma = sigma, sd = sd)
  cat("p =", p.seq[l], "done.")
  toc()
}



#### Error Decomposition ####

# use best bandwidths from first simulation 

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

error_decomp_arr = array(0, c(length(p.seq), 6, length(n.seq)))


for(i in 1:length(n.seq)){
  h_min = min_h_tibble |> filter(n == n.seq[i]) |> pull(h)
  for (j in 1:length(p.seq)) {
    tic()
    weight                 = local_polynomial_weights(p.seq[j], h_min[j], 
                                                      p.eval, parallel = T, 
                                                      m = 1, del = 0, 
                                                      parallel.environment = F)
    error_decomp_arr[j,,i] = error_decomposition(weight, n.seq[i], N, parallel = T, 
                                                 parallel.environment = F)
    cat("n =", n.seq[i], "and p =", p.seq[j], "done. ")
    toc()
  }
  
}

parallel::stopCluster(cl)
rm(weight)

#### Bandwidth Comparison ####

N = 1000
n.seq = c(50, 100, 200, 400)
p.seq = c(15, 25, 50, 75, 100)
sd = 1

p.eval = 100
H = lapply(1:length(p.seq), function(l){seq( 1, 3/p.seq[l], -0.05)})
#H = lapply(1:length(p.seq), function(l){seq(.6, 3/p.seq[l], -0.2)})


bw_comparison = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  bw_comparison[[l]] = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         Gamma   = cov_z_2rv, 
                         r.process = z_two_rv, process.arg = list(p = p.seq[l]),
                         eps.arg = list(sd = sd), 
                         m = 1)
  },
  future.seed = T)),ncol = 4)
  cat("p =", p.seq[l], "done.")
  toc()
}

parallel::stopCluster(cl)

##### 2-Fold Cross Validation #####
n = 200 

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)


tic()
two_fold_cv_2rv = k_fold_cv_simulation_2rv(N, n, p.seq[4], H[[4]], 
                                            K = 2, m = 1, w.parallel = T)
cat("p =", p.seq[4], "done.")
toc()


parallel::stopCluster(cl)


