
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



# x is single, two dim point ob evaluation
# x.design.grid are p*(p-1)/2 x 2 dim grid of Design points with x.design.grid[,1] < x.design.grid[,2]
# Resulting vector is sort like observation.grid's first column with "less"
weights.point = function(x, x.design.grid, h, K = epaK2d, m = 1, del = 0){
  L = apply(x.design.grid, 1, function(z){tcrossprod(U(z-x, h, m = m), U(z-x, h, m = m)) * K((z-x)/h)})
  if(m == 2){
    B = matrix(rowSums(L), 6, 6) # This differs from local linear estimator
    B.inv = solve(B) # TODO: can this be sped up????
    U_del = switch(del + 1, 
                    c(T, rep(F, 5)), 
                    rep(c(T, F), each = 3), 
                    rep(T, 6))
  }else{ # m = 1
    B = matrix(rowSums(L), 3, 3)
    B.inv = invert3x3(B)
    U_del = switch(del + 1, c(T,F,F), rep(T, 3))
  }
    
  by_h = 1/h
  u = switch(del+1,
             1,
             c(1, by_h, by_h),
             c(1, by_h, by_h, by_h^2, by_h^2, by_h^2))

  t(u * apply(x.design.grid, 1, function(z){
    crossprod(B.inv[,U_del], U(z-x, h, m = m) * K((z-x)/h))}))
}

# p amount of Design points -> generates equidistant grid
# p.eval are the amount of evaluation points -> best to be chosen high, till now
# only p.eval[,1] <= p.eval[,2] gets generated
# Returns array Dim p*(p-1)/2 x q x d (where q = amount of evaluation points, and d the amount of derivatives demanded)
# Result: rows are sort like observation.grid(p, "less")'s first column
# columns are sort like observation.grid(p.eval, "lesseq")' first column
# result can be evaluated by observations transformed with "observation.transform" by calling tcrossprod(t(W), Z)
local.polynomial.weights = function(p, h, p.eval, parallel = F, m = 1, 
                                    del = 0, x.design.grid = NULL,
                                    grid.type = "less", eval.type = "full", ...){
  
  if( !(grid.type %in% c("less", "without diagonal")) ){
    stop("grid type is not feasible - choose -less- oder -without diagonal-.")
  }
  if(!is.null(x.design.grid) & is.vector(x.design.grid)){ # my design grid is a seq(0,1,)
    x.design.grid = design.point.to.grid(x.design.grid, grid.type)
  }else{
    x.design.grid = observation.grid(p, grid.type)
  }
  
  if( !(eval.type %in% c("full", "diagonal")) ){
    stop("evaluation is only possible for lower diagonal -full- oder only on the -diagonal-.")
  }
  x.eval.grid = switch(eval.type, 
                       full = observation.grid(p.eval, "lesseq"), 
                       diagonal = matrix( (1:p.eval - 0.5)/p.eval, p.eval, 2)   )
  
  if(del > m){m = 2; del = 2}
  if(parallel){
    cl = makeCluster(detectCores( ) - 1)
    plan(future::cluster)
    w = future_apply(x.eval.grid, 1, FUN = weights.point, 
                     x.design.grid = x.design.grid, h = h, m = m, del = del, ...,
                     future.seed = T)
    stopCluster(cl)
  }else{
    w = apply(x.eval.grid, 1, weights.point, 
              x.design.grid = x.design.grid, 
              h = h, m = m, del = del, ...)
  } 
  
  if(del == 0){
    weights = w
  }else{
    if(grid.type == "less"){k = 2}else{k = 1}
    weights = slice.matrix(w, p*(p-1)/k, switch(del, 3, 6))
  }
 rm(w)
 L = list(design = x.design.grid, x.design = (1:p-1/2)/p, p = p, grid.type = grid.type, 
           eval = x.eval.grid, x.eval = (1:p.eval-1/2)/p.eval, p.eval = p.eval, eval.type = eval.type, 
           bandwidth = h, m = m, del = del, weights = weights)
}




###### Evaluation for weights ######
eval.weights = function(W, Z){
  if(W$grid.type == "less"){
    if(W$eval.type == "full"){
      M = matrix(0, W$p.eval, W$p.eval)
      M.up = upper.tri(M, T)
      M[M.up] = crossprod(W$weights, Z)
      M[!M.up] = t(M)[!M.up]
    }else{
      if(W$eval.type != "diagonal"){stop("eval.type needs to be -diagonal- or -full-")}
      M = crossprod(W$weights, Z)
    }
  }else{
    if(W$grid.type != "without diagonal"){stop("grid type not implemented")}
    if(W$p*(W$p-1) != length(Z)){stop("Z is not in correct representation")}
    M = crossprod(W$weights, Z)
  }
  return(M)
}

# W_Diff object from local.polynomial.weights 
eval.deriv.weights = function(W, Z){
  if(W$del == 0){
    return(eval.weights(W, Z)) 
  }else{
    if(W$eval.type == "full"){
      if(W$grid.type == "less"){
        weights_eval = apply(W$weights, 3, function(w){M = matrix(0, W$p.eval, W$p.eval)
                                                       M.up = upper.tri(M, T)
                                                       M[M.up] = crossprod(w, Z)
                                                       M[!M.up] = t(M)[!M.up]
                                                       M})
        weights_eval = array(weights_eval, dim = c(W$p.eval, W$p.eval, switch(W$del, 3, 6)))
      }else if(W$grid.type == "without diagonal"){
        weights_eval = apply(W$weights, 3, function(w){crossprod(w, Z)}) |> 
          array(dim = c(W$p.eval, W$p.eval, switch(W$del, 3, 6)))
      }else{
        stop("grid type not implemented")
      }
    }else{
      if(W$eval.type != "diagonal"){stop("eval.type needs to be -diagonal- or -full-")}
      weights_eval = apply(W$weights, 3, function(w){crossprod(w, Z)}) 
    }
  }
  return(weights_eval)
}


##### Transform Observations Y_11, ..., Y_np #####
# according to estimator in Berger/Holzmann (2024), Formula (7) and (8)
# Order corresponds to calculation in weights.lin
# Y Matrix with Observations in R^{n x p}
observation.transformation = function(Y, grid.type = "less"){ 
  p = length(Y[1,])
  n = length(Y[,1])
  Y.means = colMeans(Y)                         # pointwise mean
  Y.2 = apply(Y, 1, tcrossprod)                 # Dim (p^2 x n) (order of the col's are as.vector of the matrices)
  Z = Y.2 - as.vector(tcrossprod(Y.means))      # vector gets subtracted of all columns. Dim (p^2 x n)
  Z.mean = rowSums(Z)/(n-1)                     # Dim p^2. Reduction that is only possible in the case of synchronous observations. 
  if(grid.type == "less"){
    M2 = as.vector(upper.tri(matrix(0, p, p)))    # corresponds to entries of: observation.grid(p, "less)
 }else if(grid.type == "without diagonal"){
    M2 = as.vector(!diag(T, p, p))
 }else if(grid.type == "full"){
    M2 = T 
 }else{stop("grid type not implemented")}
  return(Z.mean[M2])                                    # in R^(p*(p-1)/2) or R^(p*(p-1)) or R^{p x p}
}




#### U ####
# can be evaluated for one single, two dim point
U = function(x, h, m = 1){
  x = x/h
  if(m == 1){
    return(c(1, x[1], x[2])) 
  }else if(m == 2){
    return(c(1, x[1], x[2], x[1]^2/2, x[1]*x[2], x[2]^2/2))
  }else{
    return("m not implemented")
  }
}



########## Computational Functions #########

observation.grid = function(p, comp = "less"){
  x = (1:p - 1/2)/p
  x.grid = expand.grid(x, x)
  if(comp == "less")
    b = x.grid[,1] < x.grid[,2]
  if(comp == "lesseq")
    b = x.grid[,1] <= x.grid[,2]
  if(comp == "gtr")
    b = x.grid[,1] > x.grid[,2]
  if(comp == "gtreq")
    b = x.grid[,1] >= x.grid[,2]
  if(comp == "without diagonal")
    b = !(x.grid[,1] == x.grid[,2])
  if(comp == "full")
    b = T
  x.grid[b, ]
}

design.point.to.grid = function(grid , comp = "less"){
  x.grid = expand.grid(grid, grid)
  if(comp == "less")
    b = x.grid[,1] < x.grid[,2]
  if(comp == "lesseq")
    b = x.grid[,1] <= x.grid[,2]
  if(comp == "gtr")
    b = x.grid[,1] > x.grid[,2]
  if(comp == "gtreq")
    b = x.grid[,1] >= x.grid[,2]
  if(comp == "without diagonal")
    b = !(x.grid[,1] == x.grid[,2])
  if(comp == "full")
    b = T
  x.grid[b, ]
}


invert3x3 = function(a) {
  b = matrix(rep(0,9), ncol=3)
  
  numb11 = (a[2,3]*a[3,2] - a[2,2]*a[3,3])
  numb21 = (a[2,1]*a[3,3] - a[2,3]*a[3,1])
  numb31 = (a[2,1]*a[3,2] - a[2,2]*a[3,1])
  numb12 = (a[1,3]*a[3,2] - a[1,2]*a[3,3])
  numb22 = (a[1,3]*a[3,1] - a[1,1]*a[3,3])
  numb32 = (a[1,1]*a[3,2] - a[1,2]*a[3,1])
  numb13 = (a[1,3]*a[2,2] - a[1,2]*a[2,3])
  numb23 = (a[1,1]*a[2,3] - a[1,3]*a[2,1])
  numb33 = (a[1,1]*a[2,2] - a[1,2]*a[2,1])
  #det1a = (a[1,3]*a[2,2]*a[3,1] - a[1,2]*a[2,3]*a[3,1] - a[1,3]*a[2,1]*a[3,2] + a[1,1]*a[2,3]*a[3,2] + a[1,2]*a[2,1]*a[3,3] - a[1,1]*a[2,2]*a[3,3])
  det1a = a[1,1]*numb11 - a[2,1]*numb12 + a[3,1]*numb13
  if(abs(det1a) < 10^{-6}){
    print("unstable matrix inversion")
  }
  det2a = -det1a
  det1a.1 = 1/det1a
  det2a.1 = 1/det2a
  
  b[1,1] = numb11 * det1a.1
  b[2,1] = numb21 * det1a.1
  b[3,1] = numb31 * det2a.1
  b[1,2] = numb12 * det2a.1
  b[2,2] = numb22 * det1a.1
  b[3,2] = numb32 * det1a.1
  b[1,3] = numb13 * det1a.1
  b[2,3] = numb23 * det1a.1
  b[3,3] = numb33 * det2a.1
  return(b)
}


supNorm = function(M){
  apply(M, 2, function(m){max(abs(m))})
}

rowMin = function(M){
  apply(M, 1, min)
}



slice.matrix = function(M, d1, d2){
  d3 = dim(M)[2]
  if(d1 * d2 != dim(M)[1])
    stop("Slicing is not of correct dimension")
  tmp_arr = array(t(M), dim = c(d3, d1, d2))
  aperm(tmp_arr, perm = c(2, 1, 3), resize = T)
}




