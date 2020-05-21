setwd("~/Dropbox/Skole/DTU/Studie/MASTER/CODE/old/ice_tank_switching_times/model/peder")
library('ctsmr')
library(latex2exp)
library(boot)
source('../code/tools.R')
source("./yum_tools.R")


source('./prep.R')



##-----------------------------------------------------------------
# MODEL 2 - both
## -------------- -------------- -------------- --------------
model_dir = 'model_2_sigmoid'
source(paste0("./",model_dir,"/model_both.R"))
source(paste0("./",model_dir,"/setprm.R"))
#plotsim(simulate(makefit(setprm(model())), newdata=ice_data),ice_data)


sim = simulate(makefit(setprm(model())), newdata=ice_data)
par(mfrow = c(1,1))
plot(ice_data$CompCap, type = "b",ylim = c(0,1))
lines(sim$output$sim$CompCap, type = "l",col = "red")



sim = simulate(makefit(setprm(model())), newdata=ice_data)
plot(sim$state$sim$X0, type = "l")

par(mfrow = c(1,1))
plot(ice_data$CompCap, type = "b",ylim = c(0,1))
lines(sim$output$sim$CompCap, type = "l",col = "red")


mod_man = model()
setprm(mod_man)
fit_manual <- list()
fit_manual$model <- mod_man
fit_manual$xm <- mod_man$ParameterValues$initial
names(fit_manual$xm) = row.names(mod_man$ParameterValues)
##-----------------------------------------------------------------
model = model()
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


## The likelihood of the initial parameters set in fit
nllikelihood(fit$xm[model$pars], fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1)

## ------------------------------------------------------------------
## Fit with different optimizers

## Fit with ctsmr
model$options$maxNumberOfEval = 1000
fitctsmr_scaled <- model$estimate(ice_data, firstorder=TRUE)
#fitctsmr <- model$estimate(ice_data, firstorder=TRUE)

fitnlminb <- nlminb(fit$xm[model$pars], nllikelihood, fit=fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,
                    control=list(iter.max=400), lower = fit$lower[model$pars], upper = fit$upper[model$pars])

fitnlminb_scaled <- nlminb(fit$xm[model$pars], nllikelihood, fit=fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,
                    control=list(iter.max=400), lower = fit$lower[model$pars], upper = fit$upper[model$pars])

# Best in model 3 -165.3592

#### ESTIMATE UNCERTAINTY WITH FINITE DIFFERENCE
summary(fitctsmr)

hessObjNlminb = finiteDiffHessian2(fitnlminb$par,r = 4)
hessObjCtsmr2 = finiteDiffHessian2(fitctsmr$xm,r = 7)

cbind(fitnlminb$par, hessObjNlminb$sd, hessObjNlminb$pvals)
cbind(fitctsmr$xm, hessObjCtsmr$sd, hessObjCtsmr$pvals)

# Confidence intervals - ctmsr
#c(exp(fitctsmr$xm['p11']),exp(fitctsmr$xm['p11'] +  c(-1,1) * fitctsmr$sd['p11']))
c(exp(fitnlminb$par['p11']) ,exp(fitnlminb$par['p11'] +  c(-1,1) * hessObjNlminb$sd['p11']))

#c(exp(fitctsmr$xm['p22']),exp(fitctsmr$xm['p22'] +  c(-1,1) * fitctsmr$sd['p22']))
c(exp(fitnlminb$par['p22']) ,exp(fitnlminb$par['p22'] +  c(-1,1) * hessObjNlminb$sd['p22']))

#c(exp(fitctsmr$xm['e11']),exp(fitctsmr$xm['e11'] +  c(-1,1) * fitctsmr$sd['e11']))
c(exp(fitnlminb$par['e11']) ,exp(fitnlminb$par['e11'] +  c(-1,1) * hessObjNlminb$sd['e11']))




fitnlminb

## Plot the simulated output from the three estimated models
insert_xm <- function(fit, xm){
  fit$xm <- xm
  names(fit$xm) <- names(xm)
  fit$xm <- fit$xm[fit$model$pars]
  return(fit)
}

## ctsmr
(fitCtsmr <- insert_xm(fit, fitctsmr$xm))
simCtsmr <- simulate(fitCtsmr, newdata=ice_data)
## nlminb
(fitNlminb <- insert_xm(fit, fitnlminb$par))

# 864.6123 
#fitNlminb$xm['T10'] = 2.82
#fitNlminb$xm['T20'] = 864.6123 

summary(fitNlminb)
#fitNlminbMod = fitNlminb
#fitNlminbMod$xm['p22'] = 1.28
#fitNlminbMod$xm['xi0'] = 0.2
nllikelihood(fitNlminb$xm[model$pars], fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1)

plot(ice_data$t, ice_data$CompCap, type="l", ylab="Heat power (kW)")
simNlminb <- simulate(fitNlminbMod, newdata=ice_data)
lines(ice_data$t, simNlminb$output$sim$CompCap, type="l", col='blue')
lines(ice_data$t, simCtsmr$output$sim$CompCap, type="l", col='red')

summary(fitnlminb)

## manual
(fitManual <- insert_xm(fit, fit_manual$xm))
simManual <- simulate(fitManual, newdata=ice_data)

#simCtsmrStoch = stochastic.simulate(fitCtsmr, c(fitCtsmr$xm['X00'], fitCtsmr$xm['X10']), data=ice_data, dt = 0.00001)
simNlminbStoch = stochastic.simulate(fitNlminb, c(fitNlminb$xm['X00'], fitNlminb$xm['X10']), data=ice_data, dt = 0.0000001)

#plot(ice_data$CompCap, type = "b",ylim = c(0.4,1))
#lines(simCtsmr$output$sim,type = "l", col = 'red')
#lines(simCtsmrStoch$output$CompCap,type = "l", col = 'green')

plot(ice_data$CompCap, type = "b",ylim = c(0.4,1))
lines(simNlminb$output$sim,type = "l", col = 'red')
lines(simNlminb$state$sim$X0/1000,type = "l", col = 'blue')
lines(simNlminbStoch$output$CompCap,type = "l", col = 'green')

lines(simNlminbStoch$state$X0/1000,type = "l", col = 'green')
##-----------------------------------------------------------------

#source("./sourceFunctions.R")
#setpar("ts", mfrow=c(5,1))
par(mfrow = c(2,1))
Xplot <- ice_data

## One-step predictions
Xplot$residualCtsmr <- Xplot$CompCap- ctsmr:::predict.ctsmr(fitCtsmr, newdata=Xplot)$output$pred$CompCap
Xplot$residualNlminb <- Xplot$CompCap - ctsmr:::predict.ctsmr(fitNlminb, newdata=Xplot)$output$pred$CompCap
Xplot$residualManual <- Xplot$CompCap - ctsmr:::predict.ctsmr(fitManual, newdata=Xplot)$output$pred$CompCap



plot(Xplot$t, Xplot$CompCap, type="l", ylab="Heat power (kW)")
lines(Xplot$t, simCtsmr$output$sim$CompCap, type="l", col=2)
lines(Xplot$t, simNlminb$output$sim$CompCap, type="l", col=3)
lines(Xplot$t, simManual$output$sim$CompCap, type="l", col=4)
#legend("topright", c("Heat power","Predicted heat power"), lty=1, col=c(1,2), bg="white")
#plot(Xplot$t, Xplot$, type="l", ylab="Indoor temperature (C)")
#plot(Xplot$t, Xplot$Te, type="l", ylab="Outdoor temperature (C)")
#plot(Xplot$t, Xplot$Ps, type="l", ylab="Solar radiation (W/m2)")
#plotTSXAxis(Xplot$t, format="%m-%d %H:%M")
plot(Xplot$t, Xplot$residualCtsmr, type="l", ylab="one-step residuals", col=2)
lines(Xplot$t, Xplot$residualNlminb, type="l", ylab="one-step residuals", col=3)
#lines(Xplot$t, Xplot$residualManual, type="l", ylab="one-step residuals", col=4)



##-----------------------------------------------------------------
# MODEL 2 - first
## -------------- -------------- -------------- --------------
model_dir = 'model_2_sigmoid'
source(paste0("./",model_dir,"/model_first.R"))
source(paste0("./",model_dir,"/setprm.R"))
#plotsim(simulate(makefit(setprm(model())), newdata=ice_data),ice_data)


sim = simulate(makefit(setprm(model())), newdata=ice_data)
par(mfrow = c(1,1))
plot(ice_data$CompCap, type = "b",ylim = c(0,1))
lines(sim$output$sim$CompCap, type = "l",col = "red")



sim = simulate(makefit(setprm(model())), newdata=ice_data)
plot(sim$state$sim$X0, type = "l")

par(mfrow = c(1,1))
plot(ice_data$CompCap, type = "b",ylim = c(0,1))
lines(sim$output$sim$CompCap, type = "l",col = "red")


mod_man = model()
setprm(mod_man)
fit_manual <- list()
fit_manual$model <- mod_man
fit_manual$xm <- mod_man$ParameterValues$initial
names(fit_manual$xm) = row.names(mod_man$ParameterValues)
##-----------------------------------------------------------------
model = model()
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


## The likelihood of the initial parameters set in fit
nllikelihood(fit$xm[model$pars], fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1)

## ------------------------------------------------------------------
## Fit with different optimizers

## Fit with ctsmr
model$options$maxNumberOfEval = 1000
fitctsmr <- model$estimate(ice_data, firstorder=TRUE)

fitnlminb <- nlminb(fit$xm[model$pars], nllikelihood, fit=fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,
                    control=list(iter.max=400), lower = fit$lower[model$pars], upper = fit$upper[model$pars])

fitnlminb

## Plot the simulated output from the three estimated models
insert_xm <- function(fit, xm){
  fit$xm <- xm
  names(fit$xm) <- names(xm)
  fit$xm <- fit$xm[fit$model$pars]
  return(fit)
}

## ctsmr
(fitCtsmr <- insert_xm(fit, fitctsmr$xm))
simCtsmr <- simulate(fitCtsmr, newdata=ice_data)
## nlminb
(fitNlminb <- insert_xm(fit, fitnlminb$par))

plot(ice_data$t, ice_data$CompCap, type="l", ylab="Heat power (kW)")
simCtsmr <- simulate(fitCtsmr, newdata=ice_data)
lines(ice_data$t, simCtsmr$output$sim$CompCap, type="l", col=3)

# 864.6123 
#fitNlminb$xm['T10'] = 2.82
#fitNlminb$xm['T20'] = 864.6123 
plot(ice_data$t, ice_data$CompCap, type="l", ylab="Heat power (kW)")
simNlminb <- simulate(fitNlminb, newdata=ice_data)
lines(ice_data$t, simNlminb$output$sim$CompCap, type="l", col=3)

summary(fitnlminb)

## manual
(fitManual <- insert_xm(fit, fit_manual$xm))
simManual <- simulate(fitManual, newdata=ice_data)

#simCtsmrStoch = stochastic.simulate(fitCtsmr, c(fitCtsmr$xm['X00'], fitCtsmr$xm['X10']), data=ice_data, dt = 0.00001)
simNlminbStoch = stochastic.simulate(fitNlminb, c(fitNlminb$xm['X00'], fitNlminb$xm['X10']), data=ice_data, dt = 0.0000001)

#plot(ice_data$CompCap, type = "b",ylim = c(0.4,1))
#lines(simCtsmr$output$sim,type = "l", col = 'red')
#lines(simCtsmrStoch$output$CompCap,type = "l", col = 'green')

plot(ice_data$CompCap, type = "b",ylim = c(0.4,1))
lines(simNlminb$output$sim,type = "l", col = 'red')
lines(simNlminb$state$sim$X0/1000,type = "l", col = 'blue')
lines(simNlminbStoch$output$CompCap,type = "l", col = 'green')

lines(simNlminbStoch$state$X0/1000,type = "l", col = 'green')
##-----------------------------------------------------------------

#source("./sourceFunctions.R")
#setpar("ts", mfrow=c(5,1))
par(mfrow = c(2,1))
Xplot <- ice_data

## One-step predictions
Xplot$residualCtsmr <- Xplot$CompCap- ctsmr:::predict.ctsmr(fitCtsmr, newdata=Xplot)$output$pred$CompCap
Xplot$residualNlminb <- Xplot$CompCap - ctsmr:::predict.ctsmr(fitNlminb, newdata=Xplot)$output$pred$CompCap
Xplot$residualManual <- Xplot$CompCap - ctsmr:::predict.ctsmr(fitManual, newdata=Xplot)$output$pred$CompCap



plot(Xplot$t, Xplot$CompCap, type="l", ylab="Heat power (kW)")
lines(Xplot$t, simCtsmr$output$sim$CompCap, type="l", col=2)
lines(Xplot$t, simNlminb$output$sim$CompCap, type="l", col=3)
lines(Xplot$t, simManual$output$sim$CompCap, type="l", col=4)
#legend("topright", c("Heat power","Predicted heat power"), lty=1, col=c(1,2), bg="white")
#plot(Xplot$t, Xplot$, type="l", ylab="Indoor temperature (C)")
#plot(Xplot$t, Xplot$Te, type="l", ylab="Outdoor temperature (C)")
#plot(Xplot$t, Xplot$Ps, type="l", ylab="Solar radiation (W/m2)")
#plotTSXAxis(Xplot$t, format="%m-%d %H:%M")
plot(Xplot$t, Xplot$residualCtsmr, type="l", ylab="one-step residuals", col=2)
lines(Xplot$t, Xplot$residualNlminb, type="l", ylab="one-step residuals", col=3)
#lines(Xplot$t, Xplot$residualManual, type="l", ylab="one-step residuals", col=4)


##-----------------------------------------------------------------
# MODEL 2 - second
## -------------- -------------- -------------- --------------
model_dir = 'model_2_sigmoid'
source(paste0("./",model_dir,"/model_second.R"))
source(paste0("./",model_dir,"/setprm.R"))
#plotsim(simulate(makefit(setprm(model())), newdata=ice_data),ice_data)


sim = simulate(makefit(setprm(model())), newdata=ice_data)
par(mfrow = c(1,1))
plot(ice_data$CompCap, type = "b",ylim = c(0,1))
lines(sim$output$sim$CompCap, type = "l",col = "red")



sim = simulate(makefit(setprm(model())), newdata=ice_data)
plot(sim$state$sim$X0, type = "l")

par(mfrow = c(1,1))
plot(ice_data$CompCap, type = "b",ylim = c(0,1))
lines(sim$output$sim$CompCap, type = "l",col = "red")


mod_man = model()
setprm(mod_man)
fit_manual <- list()
fit_manual$model <- mod_man
fit_manual$xm <- mod_man$ParameterValues$initial
names(fit_manual$xm) = row.names(mod_man$ParameterValues)
##-----------------------------------------------------------------
model = model()
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


## The likelihood of the initial parameters set in fit
nllikelihood(fit$xm[model$pars], fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1)

## ------------------------------------------------------------------
## Fit with different optimizers

## Fit with ctsmr
model$options$maxNumberOfEval = 1000
fitctsmr <- model$estimate(ice_data, firstorder=TRUE)

fitnlminb <- nlminb(fit$xm[model$pars], nllikelihood, fit=fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,
                    control=list(iter.max=400), lower = fit$lower[model$pars], upper = fit$upper[model$pars])

# Best in model 3 -165.3592


fitnlminb

## Plot the simulated output from the three estimated models
insert_xm <- function(fit, xm){
  fit$xm <- xm
  names(fit$xm) <- names(xm)
  fit$xm <- fit$xm[fit$model$pars]
  return(fit)
}

## ctsmr
(fitCtsmr <- insert_xm(fit, fitctsmr$xm))
simCtsmr <- simulate(fitCtsmr, newdata=ice_data)
## nlminb
(fitNlminb <- insert_xm(fit, fitnlminb$par))

# 864.6123 
plot(ice_data$t, ice_data$CompCap, type="l", ylab="Heat power (kW)")
lines(ice_data$t, simCtsmr$output$sim$CompCap, type="l", col='red')
#fitNlminb$xm['T10'] = 2.82
#fitNlminb$xm['T20'] = 864.6123 
plot(ice_data$t, ice_data$CompCap, type="l", ylab="Heat power (kW)")
simNlminb <- simulate(fitNlminb, newdata=ice_data)
lines(ice_data$t, simNlminb$output$sim$CompCap, type="l", col=3)

summary(fitnlminb)

## manual
(fitManual <- insert_xm(fit, fit_manual$xm))
simManual <- simulate(fitManual, newdata=ice_data)

#simCtsmrStoch = stochastic.simulate(fitCtsmr, c(fitCtsmr$xm['X00'], fitCtsmr$xm['X10']), data=ice_data, dt = 0.00001)
simNlminbStoch = stochastic.simulate(fitNlminb, c(fitNlminb$xm['X00'], fitNlminb$xm['X10']), data=ice_data, dt = 0.0000001)

#plot(ice_data$CompCap, type = "b",ylim = c(0.4,1))
#lines(simCtsmr$output$sim,type = "l", col = 'red')
#lines(simCtsmrStoch$output$CompCap,type = "l", col = 'green')

plot(ice_data$CompCap, type = "b",ylim = c(0.4,1))
lines(simNlminb$output$sim,type = "l", col = 'red')
lines(simNlminb$state$sim$X0/1000,type = "l", col = 'blue')
lines(simNlminbStoch$output$CompCap,type = "l", col = 'green')

lines(simNlminbStoch$state$X0/1000,type = "l", col = 'green')
##-----------------------------------------------------------------

#source("./sourceFunctions.R")
#setpar("ts", mfrow=c(5,1))
par(mfrow = c(2,1))
Xplot <- ice_data

## One-step predictions
Xplot$residualCtsmr <- Xplot$CompCap- ctsmr:::predict.ctsmr(fitCtsmr, newdata=Xplot)$output$pred$CompCap
Xplot$residualNlminb <- Xplot$CompCap - ctsmr:::predict.ctsmr(fitNlminb, newdata=Xplot)$output$pred$CompCap
Xplot$residualManual <- Xplot$CompCap - ctsmr:::predict.ctsmr(fitManual, newdata=Xplot)$output$pred$CompCap



plot(Xplot$t, Xplot$CompCap, type="l", ylab="Heat power (kW)")
lines(Xplot$t, simCtsmr$output$sim$CompCap, type="l", col=2)
lines(Xplot$t, simNlminb$output$sim$CompCap, type="l", col=3)
lines(Xplot$t, simManual$output$sim$CompCap, type="l", col=4)
#legend("topright", c("Heat power","Predicted heat power"), lty=1, col=c(1,2), bg="white")
#plot(Xplot$t, Xplot$, type="l", ylab="Indoor temperature (C)")
#plot(Xplot$t, Xplot$Te, type="l", ylab="Outdoor temperature (C)")
#plot(Xplot$t, Xplot$Ps, type="l", ylab="Solar radiation (W/m2)")
#plotTSXAxis(Xplot$t, format="%m-%d %H:%M")
plot(Xplot$t, Xplot$residualCtsmr, type="l", ylab="one-step residuals", col=2)
lines(Xplot$t, Xplot$residualNlminb, type="l", ylab="one-step residuals", col=3)
#lines(Xplot$t, Xplot$residualManual, type="l", ylab="one-step residuals", col=4)

