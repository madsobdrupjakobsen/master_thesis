setprm <- function(model) {
  ## Init state
  model$setParameter(  T0 = c(init=0.9 ,lb=-10     ,ub=25 ) )
  #model$setParameter(  T2 = c(init=600 ,lb=-10     ,ub=2500 ) )
  model$setParameter(  T1 = c(init=0 ,lb=-10     ,ub=2500 ) )
  model$setParameter(  T2 = c(init=0 ,lb=-10     ,ub=25 ) )
  #model$setParameter(  Tout = c(init=0.2 ,lb=-20     ,ub=25 ) )
  
  ## Set the initial value for the optimization
  model$setParameter(  a0_1 = c(init=0.019370754   ,lb=-20  ,ub=50 ) )
  model$setParameter(  a0_0 = c(init=0.023590616   ,lb=-20  ,ub=50 ) )
  model$setParameter(  a1_0 = c(init=0.867249723     ,lb=-20  ,ub=50 ) )
  model$setParameter(  a1_1 = c(init=0.151822292   ,lb=-20  ,ub=50 ) )
  model$setParameter(  a2_0 = c(init=46.261177791     ,lb=-20  ,ub=50 ) )
  model$setParameter(  a2_1 = c(init=1.206710500   ,lb=-20  ,ub=50 ) )
  
  model$setParameter(  slow = c(init=0.318488737    ,lb=0.1  ,ub=30 ) )
  
  model$setParameter(  mu0 = c(init=0.929224978   ,lb=0  ,ub=1 ) )
  model$setParameter(  mu1 = c(init= 0.683920449    ,lb=0  ,ub=1 ) )

  model$setParameter( p11 = c(init=1   ,lb=-50   ,ub=10 ) )
  model$setParameter( p22 = c(init=1   ,lb=-30   ,ub=10 ) )
  model$setParameter( p33 = c(init=1   ,lb=-30   ,ub=10 ) )
  model$setParameter( e11 = c(init=-1  ,lb=-50   ,ub=10 ) )
  
  model$setParameter(  slope = c(init=4, lb=0.1   ,ub=100))
  model$setParameter(  offset = c(init=0.5, lb=0, ub=1))
  ## Constant, the thermostatic set temperature
  #model$setParameter(  Tset = c(init=5))
  ## From data report p. 23
  ##60.95+1862.09+713.56+732.93+424.31+451.85+425.92+573.94+778.77
  #model$setParameter(  sigmoidSlope = c(init=3, lb=0.1   ,ub=100.0)) # The slope parameter of the Sigmoid function tur
  ##-----------------------------------------------------------------
  invisible(model)
}
