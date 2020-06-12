library(MASS)
library(bdsmatrix)
createMeltSteps = function(changeTimes,obsTimes){
  status = 0
  
  melt = rep(NA,n)
  melt[obsTimes < changeTimes[1]] = status
  for(i in 1:(length(changeTimes)-1)){
    status = abs(status - 1)
    melt[(changeTimes[i] <= obsTimes) & (obsTimes < changeTimes[i+1])] = status
  }
  status = abs(status - 1)
  
  melt[obsTimes > changeTimes[length(changeTimes)]] = status
  return(melt)
}

diagtool <- function(residuals){
  par(mfrow=c(3,1), mar=c(3,3,1,1), mgp=c(2,0.7,0))
  plot(residuals, type="l" , cex.lab = 1.5 , cex.axis = 1.5)
  acf = acf(residuals,lag.max = 170 , plot = FALSE)
  #print(acf)
  pacf = pacf(residuals , lag.max = 170 , plot = FALSE)
  plot(acf, cex.lab = 1.5 , cex.axis = 1.5)
  plot(pacf , cex.lab = 1.5 , cex.axis = 1.5)
}

finiteDiffHessian2 = function(optPars,fit,r){
  #browser()
  npars = length(optPars)
  
  #h = optPars * h_rel #rep(h_rel,npars)
  h = (10^(-r/3))* (1. + abs(optPars))
  
  
  # All non-diagonal elements
  hess = matrix(0,npars,npars)
  for (i in seq(1,npars)){
    for (j in seq(1,npars)){
      par = optPars
      par[i] = par[i] + h[i]
      par[j] = par[j] + h[j]
      element1 = nllikelihood(par, fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,printit=FALSE)
      
      par = optPars
      par[i] = par[i] + h[i]
      par[j] = par[j] - h[j]
      element2 = nllikelihood(par, fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,printit=FALSE)
      
      par = optPars
      par[i] = par[i] - h[i]
      par[j] = par[j] + h[j]
      element3 = nllikelihood(par, fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,printit=FALSE)
      
      par = optPars
      par[i] = par[i] - h[i]
      par[j] = par[j] - h[j]
      element4 = nllikelihood(par, fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,printit=FALSE)
      hess[i,j] =  (element1 - element2- element3 +element4) /(4*h[i]*h[j])
      
    }
  }
  # Diagonal elements
  for (i in seq(1,npars)){
    par = optPars
    par[i] = par[i] + 2*h[i]
    element1 = nllikelihood(par, fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,printit=FALSE)
    
    par = optPars
    par[i] = par[i] + h[i]
    element2 = nllikelihood(par, fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,printit=FALSE)
    
    par = optPars
    element3 = nllikelihood(par, fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,printit=FALSE)
    
    par = optPars
    par[i] = par[i] - h[i]
    element4 = nllikelihood(par, fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,printit=FALSE)
    
    par = optPars
    par[i] = par[i] - 2*h[i]
    element5 = nllikelihood(par, fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,printit=FALSE)
    
    hess[i,i] =  (-element1 + 16*element2 - 30*element3 + 16*element4 - element5) /(12*h[i]^2)
  }
  #browser()
  
  
  H_sym = 1/2 * (hess + t(hess))
  
  #print(H_sym)
  
  sd = sqrt(diag(solve(H_sym)))
  print(sd)
  names(sd) = names(optPars)
  
  optPars/sd
  df = dim(ice_data)[1] - npars
  
  
  pvals = 1 - pt(abs(optPars/sd), df)
  names(pvals) = names(optPars)
  
  return(list(H=H_sym, sd = sd, df = df, pvals = 1 - pt(abs(optPars/sd), df)))
}

finiteDiffHessian = function(optPars,r){
  npars = length(optPars)
  
  h = optPars * h_rel
  print(h)
  
  f_zero = nllikelihood(optPars, fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,printit=FALSE)
  # Compute f(x+e_i) for all i
  f_singlel = rep(0,npars)
  for (i in seq(1,npars)){
    par = optPars
    par[i] = par[i] + h[i]
    f_singlel[i] =  nllikelihood(par, fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,printit=FALSE)
  }
  
  
  f_both = matrix(0,npars,npars)
  for (i in seq(1,npars)){
    for (j in seq(1,npars)){
      par = optPars
      par[i] = par[i] + h[i]
      par[j] = par[j] + h[j]
      f_both[i,j] =  nllikelihood(par, fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,printit=FALSE)
    }
  }
  
  hess = matrix(0,npars,npars)
  for (i in seq(1,npars)){
    for (j in seq(1,npars)){
      par = optPars
      par[i] = par[i] + h[i]
      par[j] = par[j] + h[j]
      #hess[i,j] =  (f_both[i,j] - f_singlel[i] - f_singlel[j] + f_zero) /(h[i]*h[j])
      
    }
  }
  
  H_sym = 1/2 * (hess + t(hess))
  
  sd = sqrt(diag(solve(H_sym)))
  print(sd)
  names(sd) = names(optPars)
  
  optPars/sd
  df = dim(ice_data)[1] - npars
  
  
  pvals = 1 - pt(abs(optPars/sd), df)
  names(pvals) = names(optPars)
  
  return(list(H=H_sym, sd = sd, df = df, pvals = 1 - pt(abs(optPars/sd), df)))
}


insert_xm <- function(fit, xm){
  fit$xm <- xm
  names(fit$xm) <- names(xm)
  fit$xm <- fit$xm[fit$model$pars]
  return(fit)
}

buildModels = function(model_func){
  # Construct a model with initial values
  mod_man = model_func
  setprm(mod_man)
  fit_manual <- list()
  fit_manual$model <- mod_man
  fit_manual$xm <- mod_man$ParameterValues$initial
  names(fit_manual$xm) = row.names(mod_man$ParameterValues)
  
  
  ##-----------------------------------------------------------------
  # Prepare model for nlminb
  model = model_func
  setprm(model)
  ## Fit with nlminb
  model$AnalyseModel()
  ## Make a fit object to pass to simulate.ctsmr
  fit <- list()
  fit$model <- model
  fit$xm <- model$ParameterValues$initial
  fit$lower <- model$ParameterValues$lower
  fit$upper <- model$ParameterValues$upper
  names(fit$xm) <- row.names(model$ParameterValues)
  names(fit$lower) <- row.names(model$ParameterValues)
  names(fit$upper) <- row.names(model$ParameterValues)
  ## Can be nessary such that xm has the right length
  class(fit) <- "ctsmr"
  ##-----------------------------------------------------------------
  
  return(list(fit=fit, model = model, mod_man=mod_man))
}

cols.blue = rgb(0, 0.4470, 0.7410)
cols.red = rgb(0.8500, 0.3250, 0.0980)
cols.orange = rgb(0.9290, 0.6940, 0.1250)
cols.purple = rgb(0.4940, 0.1840, 0.5560)
cols.green = rgb(0.4660, 0.6740, 0.1880)
cols.lightblue = rgb(0.3010, 0.7450, 0.9330)
cols.darkred = rgb(0.6350, 0.0780, 0.1840)
