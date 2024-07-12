##### Mean functions for simulations #####
mu = function(x){
  x = 2*x - 1
  return(sin(3*pi*x)*exp(-2*abs(x)))
}


# Generate full data from a process with Brownian motions and normal distributed errors

FDA.observation = function(n, x.design, f = mu, 
                           r.process = BM, process.arg = list(), 
                           r.error = rnorm, eps.arg = list(), transformed = F, 
                           decomposed = F, ...){
  
  p = length(x.design)
  f.eval = f(x.design) # R^p --> not necessary, since it cancels out anyway
  process.arg = c(list(n = n, t = x.design), process.arg)
  process = do.call(r.process, process.arg) # R^{n x p}
  eps.arg = c(list(n = n*p), eps.arg)
  eps = matrix(do.call(r.error, eps.arg), n, p) # R^{n x p}
  Y = f.eval + process + eps
  rm(process.arg); rm(eps.arg)
  
  if(transformed)
    return(observation.transformation(Y, ...) )
  if(decomposed)
    return(list(f.eval, process, eps, Y))
  return(Y)
}

# Function to generate observations from brownian motion
BM = function(n, t, sigma = 1){
  v = diff(c(0, t))
  p = length(v)
  t(replicate(n, cumsum(rnorm(p, 0, sigma*v^0.5)) ))
}

# @t: 2 dimensional
# @sigma0: Variance of starting distribution, 0 if X0 is deterministic
# @sigma: influence of brownian motion to process
# @theta: mean reversion
cov.ou = function(t, sigma, theta, sigma0 = 0){
  sigma^2/(2*theta) * (exp(-theta*abs(t[1]-t[2])) - exp(-theta*(t[1]+t[2]))) + sigma0^2*exp(-theta*(t[1]+t[2]))
}

OU = function(n, t = seq(0, 1, len = 201), mu = 0, alpha = 1, sigma = 1, x0 = 0){
  r_ou(n, t, mu, alpha, sigma, x0)$data
}

################################################################################
###################### Hypothesis Function for optimality ######################
################################################################################


# Generate Data from a process with Covariance Kernel according to the Hypothesis in Berger/Holzmann
g = function(x, p, gamma = 2.5, L = 10, f.max = 1){
  y = numeric(length(x))
  b = abs(x) < 1
  y[b] = L*(p*f.max)^{-gamma/2} * exp( -1/(1-x[b]^2) )
  y
}

##### not needed #####
Gamma = function(z, x.design, f.max = 1, ...){
  if(is.vector(z)){
    z = matrix(z, 1, 2)
  }
  p = length(x.design)
  l = ceiling(p/2)
  x.tilde = (2*p*f.max) * ( z[,1] - (x.design[l] + x.design[l+1])/2)
  y.tilde = (2*p*f.max) * ( z[,2] - (x.design[l] + x.design[l+1])/2)
  g(x.tilde, p, ...) * g(y.tilde, p, ...)
}

process.hyp = function(n, p, sd.e = 1, sd.p = 1, f = mu, f.max = 1, ...){
  x = (1:p - 1/2)/p 
  W.n = rnorm(n, 0, sd.p)
  l = ceiling(p/2)
  x.tilde = (2*p*f.max) * ( x - (x[l] - x[l+1])/2)
  Z.x = tcrossprod(g(x.tilde, p, f.max = f.max,  ...), W.n )        # dim: p x n
  Mn = f(x)                                   # dim: 1 x p
  err = matrix(rnorm(n*p, 0, sd.e), p, n)     # dim p x n
  t(Mn + Z.x + err)
}

