#### Script to Run on Server ####

library(biLocPol)
library(future.apply)
library(tidyverse)
library(tictoc)

##### Bandwidth Evaluation #####
N = 1000
n.seq = c(50, 100, 200, 400)
p.seq = c(10, 15, 25, 35, 50, 75, 100)

p.eval = 50
H = lapply(1:length(p.seq), function(l){seq(1, 3/p.seq[l], -0.05)})

theta = 2; sigma = 3

result_list = list()

cl = parallel::makeCluster(parallel::detectCores( ) - 1)
future::plan(future::cluster)

for(l in 1:length(p.seq)){
  tic()
  result_list[[l]] = matrix(t(future_sapply(1:length(H[[l]]), FUN = function(k){bandwidth_evaluation(H[[l]][k], p.seq[l], p.eval, n.seq, N, 
                                                    cov.ou, list(theta = theta, sigma = sigma), 
                                                    OU, list(alpha = theta, sigma = sigma, x0 = 0), eps.arg = list(sd = 0.75))}, 
             future.seed = T)), ncol = 4)
  cat("p =", p.seq[l], "done.")
  toc()
}


parallel::stopCluster(cl)

result_tibble = Reduce(rbind, result_list) |>
  as_tibble() |>
  rename(n = V1, p = V2, h = V3, sup.err = V4) |>
  mutate(p = as.factor(p), n = as.factor(n))


result_tibble |>
  filter(p != 35, p != 75) |> 
  ggplot() + 
  geom_line(aes(x = h, y = sup.err, col = n, linetype = p)) + 
  ylim(c(0, 0.9))


min_h_tibble = result_tibble |>
  group_by(n, p) |> 
  slice_min(sup.err) 

sup_err_n400_plot = 
  result_tibble |> 
  filter(n == 400) |> 
  ggplot( ) + 
  geom_line(aes(x = h, y = sup.err, col = p, linetype = p)) +
  ylim(c(0, 0.53)) + 
  labs(title = "n = 400") +
  geom_point(aes(x = h + c(0, 0, 0, -0.0025, 0.0025, -0.0025, 0.0025), y = 0, col = p), min_h_tibble |> 
               filter(n == 400))

##### Cross Validation #####
p.seq = c(10, 25, 50, 100)
n = 200
p.eval = 50
H = lapply(1:length(p.seq), function(l){seq(0.8, 4/p.seq[l], -0.05)})
H
N = 1000
theta = 2; sigma = 3


Y = FDA_observation(n, (1:p - 0.5)/p, r.process = OU, process.arg = list(alpha = theta, sigma = sigma, x0 = 0), 
                    r.error = rnorm, eps.arg = list(sd = 0.25))

observation_transformation(Y)


##### Error Decomposition ###
# Use setting from before
# evaluate the optimal bw per setting
min_bw_df = bw_df |> group_by(n, p) |> 
  slice_min(sup.err)
min_bw_df
N = 1000

erg = min_bw_df |> 
  ungroup() |> 
  select(-n)|> 
  mutate(n = rep(c(50, 100, 200, 400), each = 4)) |> 
  relocate(n) |> 
  apply(1, function(par){
  p.eval = 20
  w_tmp = weights.lin(par[2], par[3], p.eval, parallel = T)
  err_tmp = error.decomposition(w_tmp, par[1], N, r.process = OU, process.arg = list(alpha = theta, sigma = sigma, x0 = 0),
                                Gamma = cov.ou, Gamma.arg = list(sigma, theta), parallel = T)
  print(par)
  return(err_tmp)
})
erg
