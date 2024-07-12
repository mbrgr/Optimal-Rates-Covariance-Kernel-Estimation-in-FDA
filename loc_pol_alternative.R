local_polynomial_weights_alternative = function (p, h, p.eval, parallel = F, m = 1, del = 0, x.design.grid = NULL, 
                                                 grid.type = "less", eval.type = "full", ...) 
{
  if (!(grid.type %in% c("less", "without diagonal"))) {
    stop("grid type is not feasible - choose -less- oder -without diagonal-.")
  }
  if (!is.null(x.design.grid) & is.vector(x.design.grid)) {
    x.design.grid = observation_grid(x = x.design.grid, 
                                     comp = grid.type)
  }
  else {
    x.design.grid = observation_grid(p = p, comp = grid.type)
  }
  if (!(eval.type %in% c("full", "diagonal"))) {
    stop("evaluation is only possible for lower diagonal -full- oder only on the -diagonal-.")
  }
  x.eval.grid = switch(eval.type, full = observation_grid(p = p.eval, 
                                                          comp = "lesseq"), diagonal = matrix((1:p.eval - 0.5)/p.eval, 
                                                                                              p.eval, 2))
  if (del > m) {
    m = 2
    del = 2
  }
  if (parallel) {

    w = future.apply::future_apply(x.eval.grid, 1, FUN = weights_point, 
                                   x.design.grid = x.design.grid, h = h, m = m, del = del, 
                                   ..., future.seed = T)
  }
  else {
    w = apply(x.eval.grid, 1, weights_point, x.design.grid = x.design.grid, 
              h = h, m = m, del = del, ...)
  }
  if (del == 0) {
    weights = w
  }
  else {
    if (grid.type == "less") {
      k = 2
    }
    else {
      k = 1
    }
    weights = slice_matrix(w, p * (p - 1)/k, switch(del, 
                                                    3, 6))
  }
  rm(w)
  L = list(design = x.design.grid, x.design = (1:p - 1/2)/p, 
           p = p, grid.type = grid.type, eval = x.eval.grid, x.eval = (1:p.eval - 
                                                                         1/2)/p.eval, p.eval = p.eval, eval.type = eval.type, 
           bandwidth = h, m = m, del = del, weights = weights)
}
