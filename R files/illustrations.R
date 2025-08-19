library(plotly)
library(tidyverse)
library(biLocPol) # please install this package from Github first. See "README.md" file for instructions

##### Description #####
# Code for images of the two covariance kernels
# Figure 1, 2, 7 and additional images of the paths

##### Load Data #####
# instead of evaluating the functions all again the results can be loaded with
load("data/illustrations.RData")

#### plotly Layout ####
back_layout = function(p, x = 2, y = 1.2, z = .2) {
  p |> layout(
    scene = list(
      camera = list(eye = list(x = x, y = y, z = z)),# controls the angle
      xaxis = list(title = list(text = "", font = list(size = 24)), tickfont = list(size = 14)),
      yaxis = list(title = list(text = "", font = list(size = 24)), tickfont = list(size = 14)),
      zaxis = list(title = list(text = "", font = list(size = 24)), tickfont = list(size = 14))),
    showlegend = F
  )
}

front_layout = function(p, x = -2, y = -1.2, z = .2) {
  p |> layout(
    scene = list(
      camera = list(eye = list(x = x, y = y, z = z)),# controls the angle
      xaxis = list(title = list(text = "", font = list(size = 24)), tickfont = list(size = 14)),
      yaxis = list(title = list(text = "", font = list(size = 24)), tickfont = list(size = 14)),
      zaxis = list(title = list(text = "", font = list(size = 24)), tickfont = list(size = 14))),
    showlegend = F
  )
}
##### Processes #####
# Figure of paths of the processes of the second process 
# not contained in the paper
set.seed(134)
p = 40
n =  5
obs = OU(n, t = c(0,(1:p - 0.5)/p), x0 = 0) 
obs_tibble = tibble(Y = as.vector(t(obs)), x = rep(c(0, (1:p - 0.5)/p), n), n = gl(n, p + 1)) 
obs_tibble |> 
  ggplot(aes(x, Y, col = n)) + 
  geom_line(aes(lty = n), linewidth = .6) + 
  geom_point(size = .9) + 
  geom_abline(slope = 0, intercept = 0, lty = 4) + 
  labs(y = NULL, x = NULL) + 
  theme(legend.position = "none", 
        text = element_text(size = 16)) 


##### Figure of (empirical) covariances #####
set.seed(513)
n = 100
p = 40
x.design = (1:p - 1/2)/p
p.eval = 150
h = 0.5

# Target function: Orntein-Uhlenbeck process with following parameters: 
sigma = 2; theta = 3
x = observation_grid(p.eval, comp = "full")[1:p.eval, 1]
Y = FDA_observation(n, x.design, 
                    r.process = OU, process.arg = list(alpha = theta, sigma = sigma, x0 = 0), 
                    eps.arg = list(sd = 0.5))
Z = observation_transformation(Y)

temp = apply(Y, 1, tcrossprod) - as.vector(tcrossprod(colMeans(Y)))     
Z.all = rowSums(temp)/(n-1)
df.all = data.frame(observation_grid(p, comp = "full"), Z.all)

# empirical covariances only 
# Figure not containes in the paper
plot_ly(df.all, x = ~Var1, y = ~Var2, z = ~Z.all, size = .4) |> add_markers()

# evaluate target covariance kernel
cov.ou.eval = matrix(apply(observation_grid(p.eval, comp = "full"), 1, function(x){
  cov_ou(x, sigma, theta)
}), p.eval, p.eval)

# plot empirical covariances with covariance kernel
###### Figure 1 ######
figure1 = plot_ly(df.all, x = ~Var1, y = ~Var2, z = ~Z.all, size = .4) |> 
  add_markers(showlegend = F) |> 
  add_surface(x = x, y = x, z = cov.ou.eval, 
              colors = c('#BF382A', '#0C4B8E'), alpha = .3, showscale = FALSE) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))
figure1

save_image(figure1 |> front_layout(), 
           file = "grafics/OU_observation_n100p40theta3sigma2sd05_front.pdf", 
           width = 600, height = 750)
save_image(figure1 |> back_layout(), 
           file = "grafics/OU_observation_n100p40theta3sigma2sd05_back.pdf", 
           width = 600, height = 750)


##### Estimation #####
# calculate local polynomial weights
p.eval = 150
H = seq(1, 3/p, -0.05)
#h_cv = k_fold_cv(Y, H, h.parallel = T, h.parallel.environment = T)
W = local_polynomial_weights(p, 0.3, p.eval, T) # Calculation of weights (takes a while, data can be loaded instead)
est = eval_weights(W, Z)

# plot with estimation
cs2 = list(c(0, 1), c("lightblue", "darkred"))

##### Figure 2 #####
figure2a = plot_ly(df.all, x = ~Var1, y = ~Var2, z = ~Z.all, size = .4) |> 
  add_surface(x = x, y = x, z = cov.ou.eval, alpha = .3, showscale = F) |> 
  add_surface(x = x, y = x, z = est, colorscale = cs2, alpha = .3, showscale = F)|> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))
figure2a |> back_layout(x = 2, y = .8, z = .6)
save_image(figure2a |> back_layout(x = 2, y = .8, z = .6), 
           file = "grafics/ou_estimate_m1_h03_sd075.pdf", 
           width = 600, height = 750)


# calculate estimator that does not mirror the results on the diagonal (without diagonal)
#h0_cv =  k_fold_cv(Y, H, m = 0, h.parallel = T, 
#                   h.parallel.environment = T)
W0 = local_polynomial_weights(p, 0.2, p.eval, T, m = 1, grid.type = "without diagonal")
Z0 = observation_transformation(Y, grid.type = "without diagonal")
est_standard = eval_weights(W0, Z0)

# compare estimate without the diagonal and without mirroring with actual kernel
figure2b = plot_ly(df.all, x = ~Var1, y = ~Var2, z = ~Z.all, size = .4) |> 
  add_surface(x = x, y = x, z = cov.ou.eval, alpha = .3, showscale = F) |> 
  add_surface(x = x, y = x, z = est_standard, colorscale = cs2, alpha = .3, showscale = F)|> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))
figure2b |> back_layout(x = 2, y = .8, z = .6)
save_image(figure2b |> back_layout(x = 2, y = .8, z = .6), 
           file = "grafics/ou_estimate_m1_h02_full_sd075.pdf", 
           width = 600, height = 750)


###### Process consisting of two random variables ######
set.seed(24)
n = 10
p = 40
z_two_rv(n, p) |> dim()
path_simulation = tibble(paths = z_two_rv(n, p) |> t() |> as.vector(),
                         x     = rep(((1:p) - 1/2)/p, n),
                         count = gl(n, p))
path_simulation |> 
  ggplot(aes(x = x, y = paths, col = count, lty = count)) + 
  geom_line(size = .6)  + 
  labs(y = NULL, x = NULL) + 
  theme(legend.position = "none", 
        text = element_text(size = 16)) 

# figure of paths of second process
# not included in paper
path_simulation |> 
  summarise(.by = x, mean = mean(paths)) |> 
  pull(mean)


# Visualization of covarinace kernel
p.eval = 150
point_grid = observation_grid(p.eval, comp = "full")
x = (1:p.eval - 0.5)/p.eval
cov_val = apply(point_grid, 1, cov_z_2rv)

# estimation with mirroring
n = 100
p = 40
set.seed(24)
Y2 = z_two_rv(n, p)
# h_2rv = k_fold_cv(Y, H, m = 1, h.parallel = T, h.parallel.environment = T)
Z2 = observation_transformation(Y2, grid.type = "less")
W01 = local_polynomial_weights(p, h = 0.1, p.eval = p.eval, m = 1, parallel = T)
est2 = eval_weights(W01, Z2)

cov_val2 = apply(point_grid, 1, cov_z_2rv)
c_val2 = matrix(cov_val2, sqrt(length(cov_val)))

figure7a = plot_ly() |>
  add_surface(x = ~x,
              y = ~x, 
              z = ~c_val2, 
              size = .4, alpha = .7) |>
  add_surface(x = ~x,
              y = ~x, 
              z = ~est2, 
              size = .4, alpha = .2, colorscale = cs2) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))
figure7a

save_image(figure7a |> front_layout(), 
           file = "grafics/2rv_estimate_m1_h01.pdf", 
           width = 600, height = 750)


# estimation without mirroring
Z_wd = observation_transformation(Y2, grid.type = "without diagonal")
# h_2rv = k_fold_cv(Y, H, m = 1, h.parallel = T, 
#                   h.parallel.environment = T)
W_wd = local_polynomial_weights(p, h = 0.1, p.eval = p.eval, m = 1, parallel = T, 
                             grid.type = "without diagonal")
est_wd = eval_weights(W_wd, Z_wd)

figure7b = plot_ly() |>
  add_surface(x = ~x,
              y = ~x, 
              z = ~c_val2, 
              size = .4, alpha = .7) |>
  add_surface(x = ~x,
              y = ~x, 
              z = ~est_wd, 
              size = .4, alpha = .7, colorscale = cs2) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))
figure7b

save_image(figure7b |> front_layout(), 
           file = "grafics/2rv_estimate_m1_h01_full.pdf", 
           width = 600, height = 750)

