setprm <- function(model) {
  ## Init state
  model$setParameter(  T0 = c(init=0.88999901      ,lb=-10     ,ub=25 ) )
  model$setParameter(  T1 = c(init=0.02016795 ,lb=-10     ,ub=2500 ) )
  model$setParameter(  T2 = c(init=-0.03208873 ,lb=-10     ,ub=25 ) )
  
  ## Set the initial value for the optimization
  model$setParameter(  a0_1 = c(init=2.537604   ,lb=-20  ,ub=50 ) )
  model$setParameter(  a0_0 = c(init=0.398642   ,lb=-20  ,ub=50 ) )
  model$setParameter(  a1_0 = c(init=0.762663   ,lb=-20  ,ub=50 ) )
  model$setParameter(  a1_1 = c(init=3.694011   ,lb=-20  ,ub=50 ) )
  model$setParameter(  a2_1 = c(init=1.383352   ,lb=-20  ,ub=50 ) )
  model$setParameter(  a2_0 = c(init=2.073256   ,lb=-20  ,ub=50 ) )
  
  model$setParameter(  slow = c(init=6.637094   ,lb=0.1  ,ub=30 ) )
  
  model$setParameter(  mu0 = c(init=0.930874   ,lb=0  ,ub=1 ) )
  model$setParameter(  mu1 = c(init=0.673448   ,lb=0  ,ub=1 ) )

  model$setParameter( p11 = c(init=-2   ,lb=-50   ,ub=10 ) )
  model$setParameter( p22 = c(init=-2   ,lb=-50   ,ub=10 ) )
  model$setParameter( p33 = c(init=-2   ,lb=-50   ,ub=10 ) )
  
  model$setParameter( e11 = c(init=-2  ,lb=-50   ,ub=10 ) )
  
  model$setParameter(  slope = c(init=5, lb=0.1   ,ub=100))
  model$setParameter(  offset = c(init=0.6, lb=0, ub=100))
  ##-----------------------------------------------------------------
  invisible(model)
}
