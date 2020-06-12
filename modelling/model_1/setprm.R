setprm <- function(model) {
  ## Init state
  model$setParameter(  X = c(init = 0.90, lb = 1E-5, ub = 100))
  
  model$setParameter(lambda0 = c(init = 0.04, lb = 1E-5, ub = 1E5))
  model$setParameter(lambda1 = c(init = 0.10, lb = 1E-5, ub = 1E5))
  model$setParameter(mu0 = c(init = 0.50, lb = 1E-5, ub = 100))
  model$setParameter(mu1 = c(init = 0.50, lb = 1E-5, ub = 100))
  model$setParameter(p11 = c(init = 1, lb = -10, ub = 1E5))
  model$setParameter(e11 = c(init = 1, lb = -10, ub = 1E5))
  ##-----------------------------------------------------------------
  invisible(model)
}
