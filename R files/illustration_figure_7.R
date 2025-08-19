
library(plotly)
library(tidyverse)
library(biLocPol) # please install this package from Github first. See "README.md" file for instructions


##### Load Data #####
# instead of evaluating the functions all again the results can be loaded with
load("data/figure_7.RData")
source("R files/functions.R")

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



cs2 = list(c(0, 1), c("lightblue", "darkred"))
figure7a = plot_ly() |>
  add_surface(x = ~x,
              y = ~x, 
              z = ~c_val2, 
              size = .4, alpha = .7, showscale = F) |>
  add_surface(x = ~x,
              y = ~x, 
              z = ~est2, 
              size = .4, alpha = .2, colorscale = cs2, showscale = F) |> 
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
              size = .4, alpha = .7, showscale = F) |>
  add_surface(x = ~x,
              y = ~x, 
              z = ~est_wd, 
              size = .4, alpha = .7, colorscale = cs2, showscale = F) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))
figure7b

save_image(figure7b |> front_layout(), 
           file = "grafics/2rv_estimate_m1_h01_full.pdf", 
           width = 600, height = 750)

save.image("data/figure_7.rdata")
