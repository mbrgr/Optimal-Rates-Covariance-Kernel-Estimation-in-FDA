##### Error Decomposition #####
library(biLocPol)

#' @param W weights
#' @param n amount of curves
#' @param N amount of simulations
#' @param Gamma covariance kernel function form [0,1]^2 to R
#' @param Gamma.arg list with further arguments that shall be passed to Gamma
#' @param lower.order.errors shall the lower order errors be computed?
#' @param f mean function -> not needed for FDA covariance estimation with synchronous design
#' @param r.process distribution of the processes Z
#' @param process.arg further arguments for the function of the process
#' @param r.error distribution of the additional errors
#' @param eps.arg further arguments for the additional errors
#' @param parallel shall the code be executed parallely?
#'
#' @return TODO
#' @importFrom stats rnorm
#' @export
#'
#' @examples
#' 0 # TODO
error_decomposition = function(W, n, N, Gamma = cov.ou, Gamma.arg = list(theta = 2, sigma = 3), 
                               lower.order.errors = F, f = mu, 
                               r.process = OU, process.arg = list(alpha = 2, sigma = 3), 
                               r.error = rnorm, eps.arg = list(sd = 0.75), 
                               parallel = T, parallel.environment = T){
  if (W$del != 0) {stop("error decomposition for derivatives not implemented (yet)")}
  
  w = W$weights
  p = W$p
  
  up = as.vector(upper.tri(diag(numeric(p))))
  dw = as.vector(lower.tri(diag(numeric(p))))
  
  code = function(useless){
    # simulate data
    # f.eval = f(x.design) # R^p --> not neccessary, since it cancels out anyway
    process.arg = c(list(n = n, t = W$x.design), process.arg)
    Process = do.call(r.process, process.arg) # R^{n x p}
    eps.arg = c(list(n = n*p), eps.arg)
    eps = matrix(do.call(r.error, eps.arg), n, p) # R^{n x p}
    Z = observation_transformation(Process + eps)  # + f.eval if mean function shall be considered
    
    # calculate parts of error decomposition 
    # (10): 1/n * sum_{i = 1}^n e_{i,j} e_{i,k}
    E = rowMeans(apply(eps, 1, tcrossprod))[up] # length p*(p-1)/2
    # T10 = as.vector(crossprod(w, E))
    T10 = eval_weights(W, E) |> abs() |> max()
      
    # (12): 1/n*sum_{i=1}^n (Z_ij Z_ik - Gamma_jk)
    Gamma.w.args = function(x){do.call(Gamma, append(list(t = x), Gamma.arg))}
    G = apply(W$design, 1, Gamma.w.args)
    ZZ = rowMeans(apply(Process, 1, tcrossprod))[up]
    T12 = W |> eval_weights(ZZ - G) |> abs() |> max()
    # T12 = as.vector(crossprod(w, ZZ - G))
    
    # (13): Mixture term: 1/n sum_{i=1}^n Z_ij e_ik + Z_ik e_ij
    temp = rowMeans(sapply(1:n, function(i){ tcrossprod(Process[i,], eps[i,]) }))
    ZE = temp[up] + temp[dw]
    T13 = W |> eval_weights(ZE) |>abs() |> max()
    # T13 = as.vector(crossprod(w, ZE))
    
    # (11): Discretization Error: Gamma_{j,k} - Gamma(x,y)  
    temp = apply(W$eval, 1, Gamma.w.args)
    T11 = max(abs((crossprod(w, G) - temp)))
    # T11 as.vector.data.frame()# T11 = as.vector(crossprod(w, G) - temp)
    
    # Overall Error 
    SUM = max(abs((crossprod(w, Z) - temp)))
    rm(temp)
    
    T14 = 0
    if(lower.order.errors){
      E2 = rowSums(sapply(1:n, function(i){
        crossprod( matrix(eps[i,], n-1, p, T), eps[-i, ] )
      }))[up]
      Z2 = rowSums(sapply(1:n, function(i){
        crossprod( matrix(Process[i,], n-1, p, T), Process[-i, ] )
      }))[up]
      ZE2 = rowSums(sapply(1:n, function(i){
        crossprod( matrix(Process[i,], n-1, p, T), eps[-i, ] )
      }))[up]
      EZ2 = rowSums(sapply(1:n, function(i){
        crossprod( matrix(eps[i,], n-1, p, T), Process[-i, ] )
      }))[up]
      T14 = as.vector(crossprod(w, (E2 + Z2 + EZ2 + ZE2)/(n*(n-1))))
    }
    
    cbind(err.2 = T10, Discr = T11, Process = T12, Mix = T13, Rest = T14, sup.err = SUM)
  }
  
  if(parallel){
    if(parallel.environment)
    {
      cl = makeCluster(detectCores( ) - 1)
      plan(future::cluster)
    }
    erg = future_sapply(1:N, code, future.seed = T)
    if (parallel.environment) { stopCluster(cl) }
  }else{
    erg = sapply(1:N, code)
  }
  rowMeans(erg)
}
  