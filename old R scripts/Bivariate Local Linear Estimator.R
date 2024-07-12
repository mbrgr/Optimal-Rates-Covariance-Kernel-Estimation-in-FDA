######### Libraries ###########
library(tidyverse)
library(locpol) # only used for epaK2d

# Parallel R-Code with Windows
library(parallel)
library(future.apply)

# Ornstein-Uhlenbeck Process
library(goffda)

# 3D Plots
library(plotly)

# Standard bivariate LocPol in R
library(interp)



###### Bivariate Loc Linear Estimator for Cov FDA #######

bivariate.loc.lin.point = function(s, x.design, Z, h, derivatives = F){
  x.design.grid = expand.grid(x.design, x.design)
  p = length(x.design)
  index = as.vector(upper.tri(matrix(0, p, p)))                # Dimension: 
  T.s = t(rbind(1, as.matrix(t(x.design.grid[index,]) - s)/h)) # p*(p-1)/2 x 3
  W.s = diag(apply(T.s[,-1], 1, epaK2d))                       # p*(p-1)/2 x p*(p-1)/2
  
  TW.s = t(crossprod(T.s, W.s))                                # p*(p-1)/2 x 3    
  inv = invert3x3(crossprod(TW.s, T.s))                        # 3 x 3
  tail = crossprod(TW.s, Z)                                    # 3 x 1
  erg = crossprod(inv, tail) * c(1, 1/h, 1/h)                  # 3
  if(derivatives){
    rownames(erg) = c("z", "dx", "dy")
    return(erg)
  }else{
    return(erg[1])
  }
}

# x.design as vector (not bivariate data)
bivariate.loc.lin = function(x.design, p.eval, Z, h, parallel = F, derivatives = F){
  x.eval = (1:p.eval - 1/2)/p.eval
  x.eval.grid = observation.grid(p.eval, "lesseq")
  if(parallel){
    cl = makeCluster(detectCores( ) - 1)
    plan(future::cluster)
    res = future_apply(x.eval.grid, 1, bivariate.loc.lin.point, 
                       x.design = x.design, Z = Z, h = h, derivatives = derivatives, 
                       future.seed = T)
    stopCluster(cl)
  }else{
    res = apply(x.eval.grid, 1, bivariate.loc.lin.point, 
          x.design = x.design, Z = Z, h = h, derivatives = derivatives)
  }
  M = matrix(0, p.eval, p.eval)
  M.up = upper.tri(M, T)
  if(derivatives){
    M[M.up] = res[1,]
    M[!M.up] = t(M)[!M.up]
    Mx = matrix(0, p.eval, p.eval)
    My = matrix(0, p.eval, p.eval)
    Mx[M.up] = res[2,]
    Mx[!M.up] = t(My)[!M.up]
    My[M.up] = res[3,]
    My[!M.up] = t(My)[!M.up]
    return(list(x.eval = x.eval, est = M, dx = Mx, dy = My, h = h))
  }
  M[M.up] = res
  M[!M.up] = t(M)[!M.up]
  return(list(x.eval = x.eval, est = M, h = h))
}



#### Calculation of the weights ####


# x is single, two dim point ob evaluation
# x.design.grid are p*(p-1)/2 x 2 dim grid of Design points with x.design.grid[,1] < x.design.grid[,2]
# Resulting vector is sort like observation.grid's first column with "less"
weights.lin.point = function(x, x.design.grid, h, K = epaK2d){
  L = apply(x.design.grid, 1, function(z){tcrossprod(U(z-x, h), U(z-x, h)) * K((z-x)/h)})
  B = matrix(rowSums(L), 3, 3)
  B.inv = invert3x3(B)
  weights = apply(x.design.grid, 1, function(z){ crossprod(B.inv[,1], U(z-x, h)) * K((z-x)/h) })
  weights 
}

# p amount of Design points -> generates equidistant grid
# p.eval are the amount of evaluation points -> best to be chosen high, till now
# only p.eval[,1] <= p.eval[,2] gets generated
# Returns Matrix of Dim p*(p-1)/2 x q (where q = amount of evaluation points)
# Result: rows are sort like observation.grid(p, "less")'s first column
# columns are sort like observation.grid(p.eval, "lesseq")' first column
# result can be evaluated by observations transformed with "observation.transform" by calling tcrossprod(t(W), Z)
weights.lin = function(p, h, p.eval, parallel = F, ...){
  
  x.design.grid = observation.grid(p)
  x.eval.grid = observation.grid(p.eval, "lesseq")
  
  if(parallel){
    cl = makeCluster(detectCores( ) - 1)
    plan(future::cluster)
    w = future_apply(x.eval.grid, 1, weights.lin.point, 
                       x.design.grid = x.design.grid, h = h, ...,
                       future.seed = T)
    stopCluster(cl)
  }else{
    w = apply(x.eval.grid, 1, weights.lin.point, 
                x.design.grid = x.design.grid, h = h, ...)
  }
  
  # res
  L = list(weights = w, design = x.design.grid, eval = x.eval.grid, 
           x.eval = (1:p.eval-1/2)/p.eval, x.design = (1:p-1/2)/p, bandwidth = h, 
           p = p, p.eval = p.eval)
}