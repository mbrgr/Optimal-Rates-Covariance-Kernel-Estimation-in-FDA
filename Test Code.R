

#### Check: weights.lin ####
n = 150
p = 45 
x.design = (1:p -1/2)/p
p.eval = 50

W = weights.lin(p, .2, p.eval, T)
W2 = local.polynomial.weights(p, .2, p.eval, T, m = 2)
Z = FDA.observation(n, x.design, transformed = T)
temp.est = eval.weights(W, Z)
quad_tmp_est = eval.weights(W2, Z)

max(abs(temp.est - observation.grid(p.eval, "full")[,1]))
max(abs(quad_tmp_est - observation.grid(p.eval, "full")[,1]))
# plot(temp.est, ylim = c(0,1))
# points(observation.grid(p.eval, "lesseq")[,1], col = "red")

##### Consistency #####
N = 1000
n.seq = c(24, 60, 120)
p.seq = n.seq/3
gamma = 2
h.seq = (log(n.seq*p.seq)/(n.seq*p.seq))^(1/(2*gamma+1))
p.eval = 50

res = matrix(0, N, length(n.seq))
for(l in 1:length(n.seq)){
  W = weights.lin(p.seq[l], h.seq[l], p.eval, T)
  temp = replicate(N, expr = {
    Z = FDA.observation(n.seq[l], (1:p.seq[l] - 0.5)/p.seq[l], transformed = T)
    crossprod(W$weights, Z)
  }, simplify = T)
  res[,l] = supNorm(temp - W$eval[,1])
  rm(temp); rm(W)
  print(l)
}
colMeans(res)*sqrt(n.seq)
colMeans(res)
boxplot(res)


##### Comparison of weights and direct estimate #####
n = 1000
p = 20
x.design = (1:p - 1/2)/p
p.eval = 40
h = 0.3

Y = FDA.observation(n, x.design, transformed = F)
Z = observation.transformation(Y)
temp = apply(Y, 1, tcrossprod) - as.vector(tcrossprod(colMeans(Y)))     
Z.all = rowSums(temp)/(n-1)
df.all = data.frame(observation.grid(p, "full"), Z.all)
plot_ly(df.all, x = ~Var1, y = ~Var2, z = ~Z.all, size = .4) |> add_markers()

W = weights.lin(p, h, p.eval, T)
est1 = eval.weights(W, Z)

est2 = full.bivariate.loc.lin(x.design, p.eval, Y, h)$est
est2
est1 - est2
max(round(est1 - est2, 5))

##### Visual #####
n = 100
p = 40
x.design = (1:p - 1/2)/p
p.eval = 100
h = 0.5

# Target function
sigma = 2; theta = 3
x = observation_grid(p.eval, comp = "full")[1:p.eval, 1]
Y = FDA_observation(n, x.design, 
                    r.process = OU, process.arg = list(alpha = theta, sigma = sigma, x0 = 0), 
                    eps.arg = list(sd = 0.5))
Z = observation_transformation(Y)

temp = apply(Y, 1, tcrossprod) - as.vector(tcrossprod(colMeans(Y)))     
Z.all = rowSums(temp)/(n-1)
df.all = data.frame(observation_grid(p, comp = "full"), Z.all)
plot_ly(df.all, x = ~Var1, y = ~Var2, z = ~Z.all, size = .4) |> add_markers()

cov.ou.eval = matrix(apply(observation_grid(p.eval, comp = "full"), 1, function(x){
  cov.ou(x, sigma, theta)
}), p.eval, p.eval)


plot_ly(df.all, x = ~Var1, y = ~Var2, z = ~Z.all, size = .4) |> 
  add_markers() |> 
  add_surface(x = x, y = x, z = cov.ou.eval, 
              colors = c('#BF382A', '#0C4B8E'), alpha = .3) |> 
  layout(scene = list(xaxis = list(title = "x"), 
                      yaxis = list(title = "y"), 
                      zaxis = list(title = "")))

W = local_polynomial_weights(p, h, p.eval, T)
est = eval_weights(W, Z)


cs2 = list(c(0, 1), c("lightblue", "darkred"))
plot_ly(df.all, x = ~Var1, y = ~Var2, z = ~Z.all, size = .4) |> 
  add_surface(x = x, y = x, z = cov.ou.eval, alpha = .3) |> 
  add_surface(x = x, y = x, z = est, colorscale = cs2, alpha = .3)|> 
  layout(scene = list(xaxis = list(title = "x"), 
                      yaxis = list(title = "y"), 
                      zaxis = list(title = "")))


full_est = full.bivariate.loc.lin(x.design, p.eval, Y, h, parallel = T, exclude.diagonal = T)$est


r_full_est = locpoly_call(x.design, p.eval, Y, h, "epanechnikov")



# bivReg_est = bivReg_call(x.design, p.eval, Y, h)

round(r_full_est - full_est, 4)
max(abs(round(r_full_est - full_est, 4)))

cs1 = list(c(0, 1), c("tan", "blue"))
cs2 = list(c(0, 1), c("lightblue", "darkred"))
cs3 = list(c(0,1 ), c("grey", "black"))
plot_ly(x = x, y = x, z = cov.ou.eval, colorscale = cs3) |> add_surface()



plot_ly(x = x, y = x, z = cov.ou.eval) |> add_surface(colorscale = cs3) |> 
  add_surface(z = est, colorscale = cs2, alpha = .5) |>
  add_markers(x = df.all$Var1, y = df.all$Var2, z = df.all$Z.all, size = .4)


plot_ly(x = ~x.design, y = ~x.design, z = ~Z.all, colors = c('#BF382A', '#0C4B8E'))  |> add_surface()
plot_ly(df.all, x = ~Var1, y = ~Var2, z = ~Z.all) |> add_markers()

plot_ly(x = x, y = x, z = full_est) |> add_surface(alpha = .5, colorscale = cs1) |>
  add_surface(r_full_est, alpha = .5)


plot_ly(x = x, y = x, z = est, colorscale = cs2) |> add_surface() |> labs(title = "Mirror")


plot_ly(x = x, y = x, z = r_full_est) |> add_surface(alpha = .5) |>
  add_surface(z = full_est, colorscale = cs1, alpha = .5) |>
  add_surface(z = est, colorscale = cs2, alpha = .5) |>
  add_surface(z = cov.ou.eval, colorscale = cs3)

plot_ly(x = x, y = x, z = est) |> add_surface() |> 
  add_surface(x = x, y = x,z = cov.ou.eval, colorscale = colorscale)

plot_ly(x = x, y = x, z = cov.ou.eval - est ) |> add_surface()
max(abs(cov.ou.eval - est))


n = 400
Z = FDA.observation(n, x.design, 
                    r.process = OU, process.arg = list(alpha = theta, sigma = sigma, x0 = 0), 
                    transformed = T)
est = eval.weights(W, Z)
plot_ly(x = x, y = x, z = est) |> add_surface() |> 
  add_surface(x = x, y = x,z = cov.ou.eval, colorscale = colorscale)
max(abs(cov.ou.eval - est))


n = 1000
Z = FDA.observation(n, x.design, 
                    r.process = OU, process.arg = list(alpha = theta, sigma = sigma, x0 = 0), 
                    transformed = T)
est = eval.weights(W, Z)
plot_ly(x = x, y = x, z = est) |> add_surface() |> 
  add_surface(x = x, y = x,z = cov.ou.eval, colorscale = colorscale)
max(abs(cov.ou.eval - est))


#### Error Decomp ####
p = 20
p.eval = 25
h = 0.3

weight_object = local_polynomial_weights(p, h, p.eval, parallel = F, m = 1, del = 0)
weight_object$del
error_decomposition(W = weight_object, n = 50, N = 100, parallel = T)

error_decomposition(W, n, N, r.process = OU, process.arg = list(alpha = theta, sigma = sigma, x0 = 0),
                    Gamma = cov.ou, Gamma.arg = list(sigma, theta), 
                    r.error = rt, eps.arg = list(df = 4), parallel = T)

system.time(error_decomposition(W, n, N, r.process = OU, process.arg = list(alpha = theta, sigma = sigma, x0 = 0),
                    Gamma = cov.ou, Gamma.arg = list(sigma, theta), 
                    r.error = rt, eps.arg = list(df = 4), parallel = T))
system.time(error_decomposition(W, n, N, r.process = OU, process.arg = list(alpha = theta, sigma = sigma, x0 = 0),
                                Gamma = cov.ou, Gamma.arg = list(sigma, theta), 
                                r.error = rt, eps.arg = list(df = 4), parallel = F))




#### Check Quadratic Estimator ####
p = 60
n = 150
p.eval = 25
parallel = T
m = 2
diff = 2
h = 0.3
K = epaK2d

x = c(.2, .3)
z = c(.4, .25)
x.design.grid = observation.grid(p)

w2 = weights.point(x, x.design.grid, h, K, m = 2, del = 2)
w1 = weights.point(x, x.design.grid, h, K, m = 1, del = 1)
head(w1)
head(w2)
weights.point(x, x.design.grid, h, K, m = 1, del = 0)


weight_object = local.polynomial.weights(p, h, 2, F, m = 2, del = 2)
weights = weight_object$weights
weights[1:10, , ]
dim(weights)
1140/6
weights_arr = array(t(weights), c(190, 2, 6), dimnames = list(NULL, NULL, c("0", "x", "y", "xx", "xy", "yy")))


W1 = weights.lin(p, h, p.eval)
W1$weights[1:6, 1:6]
W2 = local.polynomial.weights(p, h, p.eval, m = 2, del =2)
dim(W2$weights)

x.design = (1:p -1/2)/p
Z = FDA.observation(n, x.design, transformed = T, eps.arg = list(mean = 0, sd = 0))


##### Estimate only on diagonal ######
diag_weights = local.polynomial.weights(p, h, p.eval, m = 1, del = 0, eval.type = "diagonal")
diag_weights$weights[1:10, 1:6]
diag_weights$grid.type
diag_weights$weights |>  dim()
tmp = eval.weights(diag_weights, Z) # seems to work
dim(tmp)
tmp

# check eval for derivatives
diag_weights = local.polynomial.weights(p, 0.4, p.eval, m = 2, del = 2, eval.type = "diagonal")
diag_weights$weights[1:10, 1:6,]
diag_weights$weights |>  dim()
tmp = eval.deriv.weights(diag_weights, Z)
dim(tmp)
tmp

##### Test different design.grids #####
full_weights = local_polynomial_weights(p, h, p.eval, m = 2, del = 0, 
                                        grid.type = "without diagonal", eval.type = "full")
dim(full_weights$weights) 
Z_full = FDA_observation(n, x.design, transformed = T, 
                         eps.arg = list(mean = 0, sd = 0), grid.type = "without diagonal")
length(Z_full)
eval_weights(full_weights, Z_full) |> dim()

# test fpr derivatives
full_weights = local.polynomial.weights(p, h, p.eval, m = 2, del = 2, 
                                        grid.type = "without diagonal", eval.type = "diagonal")
dim(full_weights$weights) 
Z_full = FDA.observation(n, x.design, transformed = T, eps.arg = list(mean = 0, sd = 0), grid.type = "without diagonal")
length(Z_full)
eval.deriv.weights(full_weights, Z_full)


#### K-Fold cross validation ####

p = 25
p.seq = 50
n = 100
p.eval = 50
H = lapply(1:length(p.seq), function(l){seq(0.8, 4/p.seq[l], -0.05)})
H
N = 1000
theta = 2; sigma = 3


Y = FDA_observation(n, (1:p - 0.5)/p, r.process = OU, process.arg = list(alpha = theta, sigma = sigma, x0 = 0), 
                                r.error = rnorm, eps.arg = list(sd = 0.25))

observation_transformation(Y)

tic()
cl = parallel::makeCluster(parallel::detectCores() - 1)
future::plan(future::cluster)
result = lopocv_sim(N, n, p.seq[[1]], H[[1]], w.parallel = T)
parallel::stopCluster(cl)
toc()

tic()
cl = parallel::makeCluster(parallel::detectCores() - 1)
future::plan(future::cluster)
result = k_fold_cv_simulation(N, n, p.seq[[1]], H[[1]], w.parallel = T)
parallel::stopCluster(cl)
toc()

table(result) %>%  plot()
summary(result)
result %>%  boxplot()


