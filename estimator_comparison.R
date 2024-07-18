# Examine a final sample error for mirrored process

N = 1000
n.seq = c(100, 300)
p.seq = c( 50,  75)

p.eval = 75
H = lapply(1:length(p.seq), function(l){seq(1, 3/p.seq[l], -0.05)})


#### OU ####

# Parameter OU Process
theta = 2; sigma = 3
# Standard deviation for additional errors
sd = 0.75

##### m = 0 #####
m = 0

bw_comparison_OU_m0 = list()
bw_comparison_OU_m0_full = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  bw_comparison_OU_m0[[l]]      = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         cov.ou, list(theta = theta, sigma = sigma), 
                         OU, list(alpha = theta, sigma = sigma, x0 = 0), 
                         eps.arg = list(sd = sd), m = m)
  },
  future.seed = T)),ncol = 4)
  cat("p =", p.seq[l], "done.")
  
  bw_comparison_OU_m0_full[[l]] = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         cov.ou, list(theta = theta, sigma = sigma), 
                         OU, list(alpha = theta, sigma = sigma, x0 = 0), 
                         eps.arg = list(sd = sd), m = m, grid.type = "without diagonal")
  },
  future.seed = T)),ncol = 4)
  
  toc()
}

parallel::stopCluster(cl)


##### m = 1 #####
m = 1

bw_comparison_OU_m1 = list()
bw_comparison_OU_m1_full = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  bw_comparison_OU_m1[[l]]      = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         cov.ou, list(theta = theta, sigma = sigma), 
                         OU, list(alpha = theta, sigma = sigma, x0 = 0), 
                         eps.arg = list(sd = sd), m = m)
  },
  future.seed = T)),ncol = 4)
  cat("p =", p.seq[l], "done.")
  
  bw_comparison_OU_m1_full[[l]] = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         cov.ou, list(theta = theta, sigma = sigma), 
                         OU, list(alpha = theta, sigma = sigma, x0 = 0), 
                         eps.arg = list(sd = sd), m = m, grid.type = "without diagonal")
  },
  future.seed = T)),ncol = 4)
  
  toc()
}

parallel::stopCluster(cl)

##### m = 2 #####
m = 2

H2 = lapply(1:length(p.seq), function(l){seq(1, 4/p.seq[l], -0.05)})
bw_comparison_OU_m2 = list()
bw_comparison_OU_m2_full = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  bw_comparison_OU_m2[[l]]      = matrix(t(future_sapply(1:length(H2[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H2[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         cov.ou, list(theta = theta, sigma = sigma), 
                         OU, list(alpha = theta, sigma = sigma, x0 = 0), 
                         eps.arg = list(sd = sd), m = m)
  },
  future.seed = T)),ncol = 4)
  cat("p =", p.seq[l], "done.")
  
  bw_comparison_OU_m2_full[[l]] = matrix(t(future_sapply(1:length(H2[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H2[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         cov.ou, list(theta = theta, sigma = sigma), 
                         OU, list(alpha = theta, sigma = sigma, x0 = 0), 
                         eps.arg = list(sd = sd), m = m, grid.type = "without diagonal")
  },
  future.seed = T)),ncol = 4)
  
  toc()
}

parallel::stopCluster(cl)




#### 2RV process ####

n.seq = 100
p.seq =  50
##### m = 0 #####
m = 0

bw_comparison_OU_m0_2rv = list()
bw_comparison_OU_m0_full_2rv = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  bw_comparison_OU_m0_2rv[[l]]      = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         Gamma = cov_z_2rv, r.process = z_two_rv, 
                         eps.arg = list(sd = sd), m = m)
  },
  future.seed = T)),ncol = 4)
  cat("p =", p.seq[l], "done.")
  
  bw_comparison_OU_m0_full_2rv[[l]] = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N,  
                         Gamma = cov_z_2rv, r.process = z_two_rv, 
                         eps.arg = list(sd = sd), m = m, grid.type = "without diagonal")
  },
  future.seed = T)),ncol = 4)
  
  toc()
}

parallel::stopCluster(cl)
##### m = 1 #####
m = 1

bw_comparison_OU_m1_2rv = list()
bw_comparison_OU_m1_full_2rv = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  bw_comparison_OU_m1_2rv[[l]]      = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         Gamma = cov_z_2rv, r.process = z_two_rv, 
                         eps.arg = list(sd = sd), m = m)
  },
  future.seed = T)),ncol = 4)
  cat("p =", p.seq[l], "done.")
  
  bw_comparison_OU_m1_full_2rv[[l]] = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N,  
                         Gamma = cov_z_2rv, r.process = z_two_rv, 
                         eps.arg = list(sd = sd), m = m, grid.type = "without diagonal")
  },
  future.seed = T)),ncol = 4)
  
  toc()
}

parallel::stopCluster(cl)


##### m = 2 #####
m = 2

bw_comparison_OU_m2_2rv = list()
bw_comparison_OU_m2_full_2rv = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  bw_comparison_OU_m2_2rv[[l]]      = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                         Gamma = cov_z_2rv, r.process = z_two_rv, 
                         eps.arg = list(sd = sd), m = m)
  },
  future.seed = T)),ncol = 4)
  cat("p =", p.seq[l], "done.")
  
  bw_comparison_OU_m2_full_2rv[[l]] = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k)
  {
    bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N,  
                         Gamma = cov_z_2rv, r.process = z_two_rv, 
                         eps.arg = list(sd = sd), m = m, grid.type = "without diagonal")
  },
  future.seed = T)),ncol = 4)
  
  toc()
}

parallel::stopCluster(cl)

