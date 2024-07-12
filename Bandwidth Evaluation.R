#######################################################
################ Bandwidth Evaluation #################
#######################################################

bandwidth_evaluation = function(h, p, p.eval, n.seq, N, 
                                Gamma = min, Gamma.arg = list(), 
                                r.process = BM, process.arg = list(), 
                                r.error = rnorm, eps.arg = list(), m = 1, 
                                grid.type = "less"){
  
  x.design = (1:p - 0.5)/p
  s = length(h)
  x.eval.grid = observation_grid(p.eval, comp = "lesseq")
  up.eval = as.vector(upper.tri(matrix(0, p.eval, p.eval), diag = T))
  temp = as.vector(matrix( (1:p.eval-1/2)/p.eval, p.eval, p.eval, byrow = F))[up.eval]
  
  
  Gamma.w.args = function(x){do.call(Gamma, append(list(t = x), Gamma.arg))}
  G = apply(x.eval.grid, 1, Gamma.w.args)
  

    sup.err = numeric(length(n.seq))
    w = local_polynomial_weights(p, h, p.eval, F, m = m, grid.type = grid.type)
    for(k in 1:length(n.seq)){
      n = n.seq[k]
      sup.err[k] =  mean(replicate(N, {
      # f.eval = f(x.design) # R^p --> not neccessary, since it cancels out anyway
        process.arg = c(list(n = n, t = x.design), process.arg)
        Process = do.call(r.process, process.arg) # R^{n x p}
        eps.arg = c(list(n = n*p), eps.arg)
        eps = matrix(do.call(r.error, eps.arg), n, p) # R^{n x p}
        if(grid.type == "less"){
          Z = observation_transformation(Process + eps)
        } else {
          Z = observation_transformation(Process + eps, grid.type = "without diagonal")
        }
        rm(Process); rm(eps)
        apply(crossprod(w$weights, Z) - G, 2, function(m){max(abs(m))})
      }))
    }
  rm(w)
  matrix(c(n.seq, rep(c(p, h), each = length(n.seq)), sup.err), ncol =  4)
}


