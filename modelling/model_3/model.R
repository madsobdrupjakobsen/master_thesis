model <- function(){
  ## Generate a new object of class ctsm
  model <- ctsm()
  ## Add a system equation and thereby also a state
  model$addSystem(dT0 ~ T1/slow   * dt + exp(p11)*dw1)
  model$addSystem(dT1 ~ T2/slow   * dt + exp(p22)*dw2)
  model$addSystem(dT2 ~  1/slow * (-a2_1 * status* T2 -a2_0 * (1-status)* T2 - 
                                     status * a1_1 * T1 - (1-status) * a1_0 * T1 - 
                                     status * a0_1*(T0-mu1) - a0_0*(1-status)*(T0-mu0) )*dt + 
                    exp(p33)*dw3)
  
  ## Set the names of the inputs
  model$addInput(status)
  ## Set the observation equation:
  model$addObs(CompCap ~ 1./(1. + exp(-(slope * (T0 - offset)))) )
  model$setVariance(CompCap ~ exp(e11))
  return(model)
}
