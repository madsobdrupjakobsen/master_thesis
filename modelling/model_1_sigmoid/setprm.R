setprm <- function(model) {
  ## Init state
  model$setParameter(  Tin = c(init=6.5 ,lb=-10     ,ub=25 ) )
  model$setParameter(  T2 = c(init=0 ,lb=-10     ,ub=25 ) )
  model$setParameter(  Tout = c(init=0.2 ,lb=-20     ,ub=25 ) )
  ## Set the initial value for the optimization
  model$setParameter(  Ci = c(init=100   ,lb=1e-5  ,ub=500 ) )
  model$setParameter(  Cout = c(init=100   ,lb=1e-5  ,ub=500 ) )
  model$setParameter(  Cmax = c(init=0.1   ,lb=1E-5  ,ub=500 ) )
  model$setParameter(  setPoint = c(init=5.4   ,lb=-10  ,ub=10 ) )
  model$setParameter(  Tout = c(init=5.1   ,lb=0  ,ub=100 ) )
  model$setParameter( mu0 = c(init=20  ,lb=-10     ,ub=1E4) )
  model$setParameter( mu1 = c(init=-10  ,lb=-20     ,ub=20) )
  model$setParameter( lambda0 = c(init=0.09  ,lb=1E-4     ,ub=1E4) )
  model$setParameter( lambda1 = c(init=0.2 ,lb=1E-4     ,ub=1E4) )
  model$setParameter( p11 = c(init=1   ,lb=-30   ,ub=10 ) )
  model$setParameter( p22 = c(init=1   ,lb=-30   ,ub=10 ) )
  model$setParameter( e11 = c(init=-1  ,lb=-50   ,ub=10 ) )
  
  ## Constant, the thermostatic set temperature
  #model$setParameter(  Tset = c(init=5))
  ## From data report p. 23
  ##60.95+1862.09+713.56+732.93+424.31+451.85+425.92+573.94+778.77
  model$setParameter(  sigmoidSlope = c(init=3, lb=0.1   ,ub=100.0)) # The slope parameter of the Sigmoid function tur
  ##-----------------------------------------------------------------
  invisible(model)
}
