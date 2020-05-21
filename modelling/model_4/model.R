model <- function(){
  ## Generate a new object of class ctsm
  model <- ctsm()
  ## Add a system equation and thereby also a state
  #model$addSystem(dCompState ~ ( status*1/Ci*PhMax*(1/(1+exp(-1*(Tset-CompState)*sigmoidSlope))))*dt + exp(p11)*dw1 )
  #model$addSystem(dTin ~ (status*lambda1*(mu1-Tin) + 1/Ci * (Tout - Tin) )*dt + exp(p11)*dw1 )
  #model$addSystem(dTin ~ (status*lambda1*(mu1-Tin) + (1-status)*lambda0*(mu0-Tin) + 1/Ci * (T2 - Tin) )*dt + exp(p11)*dw1)
  #model$addSystem(dT2 ~ (-Cmax * 1/(1+exp((setPoint-T2)*sigmoidSlope)) + 1/Ci * (Tin - T2) + 1/Cout * (Tout - T2)  )*dt + exp(p22)*dw2)
  #model$addSystem(dTout ~ ( 1/C_comp * (-1/(1+exp((5-Tout)*sigmoidSlope))+0.5) + 1/Ci * (Tin - Tout) + 0.001 * (15.0 - Tout))*dt + exp(p22)*dw2 )
  
  model$addSystem(dT0 ~ T1/slow   * dt + exp(p11)*dw1)
  model$addSystem(dT1 ~ T2/slow   * dt + exp(p22)*dw2)
  #model$addSystem(dT2 ~ (-4 * T2 - 4 * T1 - 1 * (T0+3.8))  * dt + exp(p33)*dw3)
  model$addSystem(dT2 ~  1/slow * (-a2_1 * status* T2 -a2_0 * (1-status)* T2 - status * a1_1 * T1 - (1-status) * a1_0 * T1 + status * a0_1*(mu1-T0) + a0_0*(1-status)*(mu0-T0) )*dt + exp(p33)*dw3)
  #model$addSystem(dT2 ~ (8 - T2)  * dt + exp(p33)*dw3)
  #model$addSystem(dT1 ~ (-2 * xi1 * omega1 * (T1) - omega1^2 * ((-T2)-(-mu1))) * dt + exp(p11)*dw1) #(status*lambda1*(mu1-Tin)  )*dt + exp(p11)*dw1 )
  #model$addSystem(dT1 ~ (-2 * xi0 * omega0 * (T1) - omega0^2 * ((-T2)-(-mu0))) * dt + exp(p11)*dw1) #(status*lambda1*(mu1-Tin)  )*dt + exp(p11)*dw1 )
  
  #model$addSystem(dT1 ~ (-(T2-mu0)) * dt + exp(p11)*dw1) #(status*lambda1*(mu1-Tin)  )*dt + exp(p11)*dw1 )
  
  
  
  #model$addSystem(dT1 ~ (status*lambda1*(mu1-T1) + (1-status)*lambda0*(mu0-T1) - 1/C12 * (T1 - T2)) * dt + exp(p11)*dw1) #(status*lambda1*(mu1-Tin)  )*dt + exp(p11)*dw1 )
  #model$addSystem(dT2 ~ ( 1/C12 * (T1 - T2) - 1/C_comp * 1/(1+exp((setPoint-T2)*sigmoidSlope)) )*dt + exp(p22)*dw2 )
  #model$addSystem(dT3 ~ ( 1/C23 * (T2 - T3) )*dt + exp(p33)*dw3 )
  ## Set the names of the inputs
  model$addInput(status)
  ## Set the observation equation: Ti is the state, Ph is now the output
  model$addObs(CompCap ~ T0) #(1/(1+exp((offset-T0)*slope))) )
  #model$addObs(CompCap ~ T0) #1./(1. + exp(-(slope * (T3 - offset)))) )  #(1/(1+exp((setPoint-T2)*sigmoidSlope)))) 
  ## Set the variance of the measurement error
  model$setVariance(CompCap ~ exp(e11))
  return(model)
}
