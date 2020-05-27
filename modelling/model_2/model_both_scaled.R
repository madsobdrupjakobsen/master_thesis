model <- function(){
  ## Generate a new object of class ctsm
  model <- ctsm()
  ## Add a system equation and thereby also a state
  #model$addSystem(dT1 ~ (status * (-2 * xi1 * omega1 * (T1) - omega1^2 * ((-T2)-(-mu1))) + 
  #                         (1-status)*(-2 * xi0 * omega0 * (T1) - omega0^2 * ((-T2)-(-mu0))) )  * dt + exp(p11)*dw1)
  #model$addSystem(dT2 ~ (-T1) * dt + exp(p22)*dw2)
  
  model$addSystem(dX0 ~ (X1) / 889.190900140 * dt + exp(p11)/889.190900140*dw1)
  model$addSystem(dX1 ~ (status * (-2 * xi1 * omega1 * (X1) - omega1**2 * (889.190900140*X0-698.861205845*mu1) ) + 
                           (1-status)*(-2 * xi0 * omega0 * (X1) - omega0**2 * (889.190900140*X0-961.996347735*mu0) ))  * dt + exp(p22)*dw2)
  
  
  ## Set the names of the inputs
  model$addInput(status)
  ## Set the observation equation: Ti is the state, Ph is now the output
  #model$addObs(CompCap ~ (1/1000 * T2) ) #(1/(1+exp((setPoint-T2)*sigmoidSlope)))) 
  model$addObs(CompCap ~ 1./(1. + exp(-(slope * (889.190900140*X0 - 592.010123492*offset)))) )  #(1/(1+exp((setPoint-T2)*sigmoidSlope)))) 
  ## Set the variance of the measurement error
  model$setVariance(CompCap ~ exp(e11))
  return(model)
}
