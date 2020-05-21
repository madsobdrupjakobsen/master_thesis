model <- function(){
  ## Generate a new object of class ctsm
  model <- ctsm()
  ## Add a system equation and thereby also a state
  #model$addSystem(dCompState ~ ( status*1/Ci*PhMax*(1/(1+exp(-1*(Tset-CompState)*sigmoidSlope))))*dt + exp(p11)*dw1 )
  #model$addSystem(dTin ~ (status*lambda1*(mu1-Tin) + 1/Ci * (Tout - Tin) )*dt + exp(p11)*dw1 )
  model$addSystem(dTin ~ (status*lambda1*(mu1-Tin) + (1-status)*lambda0*(mu0-Tin) )*dt + exp(p11)*dw1)
  #model$addSystem(dTout ~ ( 1/C_comp * (-1/(1+exp((5-Tout)*sigmoidSlope))+0.5) + 1/Ci * (Tin - Tout) + 0.001 * (15.0 - Tout))*dt + exp(p22)*dw2 )
  ## Set the names of the inputs
  model$addInput(status)
  ## Set the observation equation: Ti is the state, Ph is now the output
  model$addObs(CompCap ~ (1/(1+exp((5-Tin)*sigmoidSlope)))) 
  ## Set the variance of the measurement error
  model$setVariance(CompCap ~ exp(e11))
  return(model)
}