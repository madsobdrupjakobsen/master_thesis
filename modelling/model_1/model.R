model <- function(){
  ## Generate a new object of class ctsm
  model <- ctsm()
  
  ## Add a system equation and thereby also a state
  model$addSystem(dX ~ ((1 -status) * lambda0 * (mu0 - X) +  status * lambda1 * (mu1 - X))*dt + exp(p11)*dw1)

  ## Set the names of the inputs
  model$addInput(status)
  
  ## Set the observation equation:
  model$addObs(CompCap ~ X )
  
  ## Set the variance of the measurement error
  model$setVariance(CompCap ~ exp(e11))
  return(model)
}
