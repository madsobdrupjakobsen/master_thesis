setprm <- function(model) {
  ## Init state
  model$setParameter(  X1 = c(init=6.5 ,lb=-10     ,ub=25 ) )
  #model$setParameter(  T2 = c(init=600 ,lb=-10     ,ub=2500 ) )
  model$setParameter(  X0 = c(init=890 ,lb=-10     ,ub=2500 ) )
  model$setParameter(  T3 = c(init=0 ,lb=-10     ,ub=25 ) )
  model$setParameter(  Tout = c(init=0.2 ,lb=-20     ,ub=25 ) )
  
  
  ## Set the initial value for the optimization
  model$setParameter(  omega0 = c(init=0.04   ,lb=0.001  ,ub=0.2 ) )
  model$setParameter(  xi0 = c(init=0.2   ,lb=1e-5  ,ub=2 ) )
  model$setParameter( mu0 = c(init=950  ,lb=0     ,ub=1000) )
  
  
  model$setParameter(  omega1 = c(init=0.2   ,lb=1e-5  ,ub=0.6 ) )
  model$setParameter(  xi1 = c(init=0.3   ,lb=1e-5  ,ub=0.9 ) )
  model$setParameter( mu1 = c(init=650  ,lb=-0     ,ub=1000) )
  
  
  model$setParameter(  C12 = c(init=100   ,lb=1e-5  ,ub=2000 ) )
  model$setParameter(  C23 = c(init=1   ,lb=1e-5  ,ub=2000 ) )
  model$setParameter(  Cout = c(init=100   ,lb=1e-5  ,ub=500 ) )
  model$setParameter(  C_comp = c(init=10   ,lb=1E-5  ,ub=500 ) )
  model$setParameter(  setPoint = c(init=5.4   ,lb=-10  ,ub=10 ) )
  model$setParameter(  Tout = c(init=5.1   ,lb=0  ,ub=100 ) )
  
  
  model$setParameter( lambda0 = c(init=0.09  ,lb=1E-4     ,ub=1E4) )
  model$setParameter( lambda1 = c(init=0.2 ,lb=1E-4     ,ub=1E4) )
  model$setParameter( p11 = c(init=1   ,lb=-50   ,ub=10 ) )
  model$setParameter( p22 = c(init=1   ,lb=-30   ,ub=10 ) )
  model$setParameter( p33 = c(init=1   ,lb=-30   ,ub=10 ) )
  model$setParameter( e11 = c(init=-1  ,lb=-50   ,ub=10 ) )
  
  model$setParameter(  slope = c(init=1/90, lb=1/200   ,ub=1/30))
  model$setParameter(  offset = c(init=640, lb=400, ub=800))
  ## Constant, the thermostatic set temperature
  #model$setParameter(  Tset = c(init=5))
  ## From data report p. 23
  ##60.95+1862.09+713.56+732.93+424.31+451.85+425.92+573.94+778.77
  #model$setParameter(  sigmoidSlope = c(init=3, lb=0.1   ,ub=100.0)) # The slope parameter of the Sigmoid function tur
  ##-----------------------------------------------------------------
  invisible(model)
}
