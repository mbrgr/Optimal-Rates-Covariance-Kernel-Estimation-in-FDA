#### K-Fold Cross Validation for bivariate local polynomial estimator ####

lopocv = function(Y, h.seq, m = 2, h.parallel = F, ...){
  p = length(Y[1,])
  n = length(Y[,1])
  Y = apply(Y, 2, function(x){x - mean(x)})

  help = function(h){
    w_h = local_polynomial_weights(p, h, p.eval = p, m, parallel = F, ...)
    max_diff = numeric(n)
    for(l in 1:n){
      Z_l         = tcrossprod(Y[l,], Y[l,])
      Z_minus_l   = observation_transformation(Y[-l,])
      estimate    = eval_weights(w_h, Z_minus_l)
      max_diff[l] = max( abs(Z_l - estimate)[!as.logical(diag(p))] )
    }
    mean(max_diff)
  }
  
  if (h.parallel) {
    mean_sup = future_sapply(h.seq, help, future.seed = T)
  } else {
    mean_sup = sapply(h.seq, help)
  }
  h.seq[which.min(mean_sup)]

}

lopocv_sim = function(N, n, p, h.seq, m = 1, w.parallel = F,
                      theta = 2, sigma = 3, sd = 0.75, ...){

  
  help = function(h){

    w_h = local_polynomial_weights(p = p, h = h, p.eval = p, parallel = w.parallel, m = m, parallel.environment = F,...)
    max_diff = numeric(n)

    future_replicate(N, {
      Y = FDA_observation(n, x.design = (1:p - 0.5)/p, f = biLocPol::mu,
                          r.process = OU, process.arg = list(alpha = theta, sigma = sigma, x0 = 0), eps.arg = list(sd = sd))
      Y = apply(Y, 2, function(x){x - mean(x)})
      for(l in 1:n){
        Z_l         = tcrossprod(Y[l,], Y[l,])
        Z_minus_l   = observation_transformation(Y[-l,])
        estimate    = eval_weights(w_h, Z_minus_l)
        max_diff[l] = max( abs(Z_l - estimate))
      }
      rm(Z_l);rm(Z_minus_l)
      mean(max_diff)
    }, future.seed = T)
  }
  mean_sup = sapply(h.seq, help)
  h.seq[apply(mean_sup, 1, which.min)]
}


k_fold_cv = function(Y, h.seq, K = 2, h.parallel = F, m = 1, ...){

  p = length(Y[1,])
  n = length(Y[,1])

  grp = sample(rep(1:K, ceiling(n/K)), n)

  help = function(h){
    w_h = local_polynomial_weights(p, h, p.eval = p, m, ...) # std.val
    max_diff = numeric(K)
    for(kk in 1:K){
      test_grp     = matrix(observation_transformation(Y[grp == kk,], grid.type = "full"), p, p)
      train_grp    = observation_transformation(Y[grp != kk,])
      pred         = eval_weights(w_h, train_grp)
      max_diff[kk] = max(abs(test_grp - pred)[!as.logical(diag(p))])
    }
    mean(max_diff)
  }

  if(h.parallel){
    cl = parallel::makeCluster(parallel::detectCores( ) - 1)
    future::plan(future::cluster)
    mean_sup = future_sapply(h.seq, help, future.seed = T)
    parallel::stopCluster(cl)
  }else{
    mean_sup = sapply(h.seq, help)
  }
  h.seq[which.min(mean_sup)]
}



k_fold_cv_simulation = function(N, n, p, h.seq, K = 5, m = 1, w.parallel = T, theta = 2, sigma = 3, sd = 0.75){

  grp = sample(rep(1:K, ceiling(n/K)), n)

  help = function(h){
    w_h = local_polynomial_weights(p, h, p.eval = p, m = m, parallel = w.parallel, parallel.environment = F)
    max_diff = numeric(K)

    future_replicate(N, {
      Y = FDA_observation(n, x.design = (1:p - 0.5)/p, f = biLocPol::mu,
                          r.process = OU, process.arg = list(alpha = theta, sigma = sigma, x0 = 0), eps.arg = list(sd = sd))
      for(kk in 1:K){
        test_grp     = matrix(observation_transformation(Y[grp == kk,], grid.type = "full"), p, p)
        train_grp    = observation_transformation(Y[grp != kk,], grid.type = "less")
        pred         = eval_weights(w_h, train_grp)
        max_diff[kk] = max(abs((test_grp - pred)[!as.logical(diag(p))]))
      }
      mean(max_diff)
    }, future.seed = T)
  }


  mean_sup = sapply(h.seq, help)
  
  h.seq[apply(mean_sup, 1, which.min)]
}


# k-fold simulation with process thats inspired by degras paper ##



k_fold_cv_simulation_2rv = function(N, n, p, h.seq, K = 5, m = 1, 
                                    w.parallel = T){
  
  grp = sample(rep(1:K, ceiling(n/K)), n)
  
  help = function(h){
    w_h = local_polynomial_weights(p, h, p.eval = p, m = m, 
                                   parallel = w.parallel, 
                                   parallel.environment = F)
    max_diff = numeric(K)
    
    future_replicate(N, {
      Y = FDA_observation(n, x.design = (1:p - 0.5)/p,  
                          f = biLocPol::mu, rprocess = z_two_rv, 
                          eps.arg = list(sd = 0.75))
      for(kk in 1:K){
        test_grp     = matrix(observation_transformation(Y[grp == kk,], grid.type = "full"), p, p)
        train_grp    = observation_transformation(Y[grp != kk,], grid.type = "less")
        pred         = eval_weights(w_h, train_grp)
        max_diff[kk] = max(abs((test_grp - pred)[!as.logical(diag(p))]))
      }
      mean(max_diff)
    }, future.seed = T)
  }
  
  
  mean_sup = sapply(h.seq, help)
  
  h.seq[apply(mean_sup, 1, which.min)]
}
