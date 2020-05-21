setprm <- function(model) {
  ## Init state
  model$setParameter(  X = c(init = 0.90, lb = 1E-5, ub = 100))
  
  model$setParameter(lambda0 = c(init = 0.04, lb = 1E-5, ub = 1E5))
  model$setParameter(lambda1 = c(init = 0.10, lb = 1E-5, ub = 1E5))
  model$setParameter(mu0 = c(init = 0.50, lb = 1E-5, ub = 100))
  model$setParameter(mu1 = c(init = 0.50, lb = 1E-5, ub = 100))
  model$setParameter(p11 = c(init = 1, lb = -10, ub = 1E5))
  model$setParameter(e11 = c(init = 1, lb = -10, ub = 1E5))
  ## Constant, the thermostatic set temperature
  #model$setParameter(  Tset = c(init=5))
  ## From data report p. 23
  ##60.95+1862.09+713.56+732.93+424.31+451.85+425.92+573.94+778.77
  #model$setParameter(  sigmoidSlope = c(init=3, lb=0.1   ,ub=100.0)) # The slope parameter of the Sigmoid function tur
  ##-----------------------------------------------------------------
  invisible(model)
}
