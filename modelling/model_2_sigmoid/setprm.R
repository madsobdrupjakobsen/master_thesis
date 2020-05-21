setprm <- function(model) {
  ## Init state
  model$setParameter(  X1 = c(init=6.5 ,lb=-10     ,ub=25 ) )
  model$setParameter(  X0 = c(init=890 ,lb=-10     ,ub=2500 ) )
  
  ## Set the initial value for the optimization
  model$setParameter(  omega0 = c(init=0.0572797   ,lb=0.001  ,ub=10 ) )
  model$setParameter(  xi0 = c(init=0.1897310   ,lb=1e-5  ,ub=2 ) )
  model$setParameter( mu0 = c(init=962.5912549/700  ,lb=0     ,ub=2) )
  
  
  model$setParameter(  omega1 = c(init=0.3629961   ,lb=1e-5  ,ub=10 ) )
  model$setParameter(  xi1 = c(init=0.4967667   ,lb=1e-5  ,ub=0.9 ) )
  model$setParameter( mu1 = c(init=698.1883135/500  ,lb=1e-5     ,ub=2) )
  
  #model$setParameter( p11 = c(init=2.8219038   ,lb=-50   ,ub=10 ) )
  model$setParameter( p11 = c(init=-2   ,lb=-50   ,ub=10 ) )
  model$setParameter( p22 = c(init=-2   ,lb=-30   ,ub=10 ) )
  model$setParameter( e11 = c(init=-10  ,lb=-50   ,ub=10 ) )
  
  model$setParameter(  slope = c(init=1/140, lb=1/200   ,ub=1/30))
  model$setParameter(  offset = c(init=591.2794245/500, lb=0, ub=2))

  ##-----------------------------------------------------------------
  invisible(model)
}
