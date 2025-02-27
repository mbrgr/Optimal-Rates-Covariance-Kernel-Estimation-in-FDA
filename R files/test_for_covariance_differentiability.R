library(biLocPol)
library(future.apply)
library(tidyverse)

N = 100

n = 500
p = 100
p_eval = 50
h = 0.2
x = (1:p-0.5)/p

w  = biLocPol::local_polynomial_weights(p, h, p_eval, parallel = F, m = 2, del = 2)
  

cov_z_2rv = function(x, y){
  4/9*sin(pi*x)*sin(pi*y) + 
    8/9*cos(1.25*pi*x)*cos(1.25*pi*y)
}

partial01_cov_z_2rv = function(x,y){
  4/9*pi*sin(pi*x)*cos(pi*y) -
    10/9*pi*cos(1.25*pi*x)*sin(1.25*pi*y)
}

partial11_cov_z_2rv = function(x,y){
  4/9*pi^2*cos(pi*x)*cos(pi*y) +
    25/18*pi^2*sin(1.25*pi*x)*sin(1.25*pi*y)
}



true_var = function(t,s){
  cov_z_2rv(t,s) * partial11_cov_z_2rv(t,s) + partial01_cov_z_2rv(t,s)*partial01_cov_z_2rv(s,t)
}

t_var = ((1:p-0.5)/p) |> sapply(function(x){true_var(x,x)})
plot((1:p-0.5)/p, t_var, lty = 2, col = 2, type = "l")
  

count = numeric(p_eval)
c2 = 0

for(i in 1:50){
  #Y1 = biLocPol::z_2rv(n, p)
  Y1 = BM(n, (1:p-0.5)/p)
  Y1 = Y1 |> biLocPol::observation_transformation()
  #Y2 = BM(n, (1:p-0.5)/p)
  #Y2 = biLocPol::z_2rv(n, p)
  #Y2 = Y2 |> biLocPol::observation_transformation()
  #Y3 = biLocPol::z_2rv(n, p)
  #Y3 = Y3 |> biLocPol::observation_transformation()
  
  
  estimation1 = w |> biLocPol::eval_weights(Y1)
  #estimation2 = w |> biLocPol::eval_weights(Y2)
  #estimation3 = w |> biLocPol::eval_weights(Y3)
  
  # variance
  var = diag(estimation1[,,1]) * diag(estimation1[,,5]) + diag(estimation1[,,3]/2 + estimation1[,,2]/2)^2
  #var = diag(estimation3[,,1]) * diag(estimation3[,,5]) + diag(estimation3[,,3]/2 + estimation3[,,2]/2)^2
  #plot((1:p-0.5)/p, var, type = "l")
  
  #par(mfrow = c(1,1))
 # 
  #plot(var)
  #plot(diag(estimation2[,,3]))
  plot(diag(estimation1[,,2] - estimation1[,,3]))
  count = count + as.numeric(sqrt(n)*abs(diag(estimation1[,,2] - estimation1[,,3])/sqrt(var)) > qnorm(0.95))
 # c2 = c2 + any( sqrt(n) * abs(diag(estimation1[,,2] - estimation1[,,3]) / sqrt(2*var) ) > m_quant, na.rm = T)
  print(i)
}

count/50
c2/50

cov_test = function (x,y) {
  cov_z_2rv(x, y) * partial11_cov_z_2rv(x, y) + partial01_cov_z_2rv(x,x) * partial01_cov_z_2rv(y,y)
}

cov_mat = outer(x, x, cov_test)
  
  
cor_mat = t(cov_mat / sqrt( diag(cov_mat) ))/ sqrt(diag(cov_mat))
cor_mat[1:5, 1:5]

plotly::plot_ly(x = ~w$x.eval, y = ~w$x.eval, z = ~estimation1[,,1]) |> 
  plotly::add_surface()
plotly::plot_ly(x = ~w$x.eval, y = ~w$x.eval, z = ~estimation1[,,2]) |> 
  plotly::add_surface()
plotly::plot_ly(x = ~w$x.eval, y = ~w$x.eval, z = ~estimation1[,,3]) |> 
  plotly::add_surface()


##### bootstrap #####
N = 5000
bootstrap_sample = MASS::mvrnorm(N, numeric(p), cor_mat)
max_distribution = bootstrap_sample |> apply(1, function(x){max(abs(x))})
m_quant = quantile(max_distribution, p = 0.90)
m_quant




########## 


test_funct