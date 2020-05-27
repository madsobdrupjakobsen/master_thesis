setprm <- function(model) {
  ## Init state
  model$setParameter(  X1 = c(init=-2.068630356 ,lb=-10     ,ub=25 ) )
  model$setParameter(  X0 = c(init=890.190900140/889.190900140 ,lb=-10/889.190900140     ,ub=2500/889.190900140 ) )
  
  ## Set the initial value for the optimization
  model$setParameter(  omega0 = c(init=0.054302070   ,lb=0.001  ,ub=10 ) )
  model$setParameter(  xi0 = c(init=0.131779060  ,lb=1e-5  ,ub=2 ) )
  model$setParameter( mu0 = c(init=1.0000  ,lb=1e-5     ,ub=1000/961.996347735) )
  
  
  model$setParameter(  omega1 = c(init=0.199631822    ,lb=1e-5  ,ub=10 ) )
  model$setParameter(  xi1 = c(init=0.613668062  ,lb=1e-5  ,ub=0.9 ) )
  model$setParameter( mu1 = c(init=1  ,lb=1e-5     ,ub=1000/698.861205845) )
  
  #model$setParameter( p11 = c(init=2.8219038   ,lb=-50   ,ub=10 ) )
  model$setParameter( p11 = c(init=2.720712912   ,lb=-50   ,ub=10 ) )
  model$setParameter( p22 = c(init=-0.501648509   ,lb=-30   ,ub=10 ) )
  model$setParameter( e11 = c(init=-9.935998459  ,lb=-50   ,ub=10 ) )
  
  model$setParameter(  slope = c(init=0.007309102, lb=1/200   ,ub=1/30))
  model$setParameter(  offset = c(init=1, lb=400/592.010123492, ub=800/592.010123492))
  
  ##-----------------------------------------------------------------
  invisible(model)
}
