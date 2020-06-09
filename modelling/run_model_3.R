setwd("~/Dropbox/Skole/DTU/Studie/MASTER/master_thesis/modelling")
library('ctsmr')
library(latex2exp)
library(boot)
source('./src/tools.R')
source('./src/ldf.R')
source('./src/leaveOneOut.R')
source("./src/yum_tools.R")
source('./src/prep.R')



# New subsetting
# Subsampling
nSkip = 50
Idx = seq(1,N,nSkip)
ice_data = ice_data_full[Idx,]

##-----------------------------------------------------------------
# MODEL 3 - Results changes with p11 = -12 and p11 = =-2 as initial. Best with -2
## -------------- -------------- -------------- --------------
model_dir = 'model_3'
source(paste0("./",model_dir,"/model.R"))
source(paste0("./",model_dir,"/setprm.R"))

# Simulate based on initial parameter values
sim = simulate(makefit(setprm(model())), newdata=ice_data)
par(mfrow = c(1,1))
plot(ice_data$CompCap, type = "b",ylim = c(0,1))
lines(sim$output$sim$CompCap, type = "l",col = "red")

##-----------------------------------------------------------------
# Prepare models
fit_and_models = buildModels(model())
fit = fit_and_models$fit
model = fit_and_models$model
mod_man = fit_and_models$mod_man
##-----------------------------------------------------------------

## The likelihood of the initial parameters set in fit
nllikelihood(fit$xm[model$pars], fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,printit = FALSE)


##-----------------------------------------------------------------
## ESTIMATE PARAMETERS
## -------------- -------------- -------------- --------------
## Fit with different optimizers

## Fit with ctsmr
model$options$maxNumberOfEval = 1000
fitctsmr <- model$estimate(ice_data, firstorder=TRUE)

## Fit with nlminb
fitnlminb <- nlminb(fit$xm[model$pars], nllikelihood, fit=fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,
                    control=list(iter.max=400), lower = fit$lower[model$pars], upper = fit$upper[model$pars])

# UNCERTAINTY
hessObjNlminb = finiteDiffHessian2(fitnlminb$par,fit,r = 40)
hessObjCtsmr = finiteDiffHessian2(fitctsmr$xm,r = 16)

fitnlminb

# Confidence intervals - ctmsr
c(exp(fitctsmr$xm['p11']),exp(fitctsmr$xm['p11'] +  c(-1,1) * fitctsmr$sd['p11']))
c(exp(fitctsmr$xm['p22']),exp(fitctsmr$xm['p22'] +  c(-1,1) * fitctsmr$sd['p22']))
c(exp(fitctsmr$xm['p33']),exp(fitctsmr$xm['p33'] +  c(-1,1) * fitctsmr$sd['p33']))
c(exp(fitctsmr$xm['e11']),exp(fitctsmr$xm['e11'] +  c(-1,1) * fitctsmr$sd['e11']))
sqrt(c(exp(fitctsmr$xm['e11']),exp(fitctsmr$xm['e11'] +  c(-1,1) * fitctsmr$sd['e11'])))

exp(fitnlminb$par['p11'])
exp(fitnlminb$par['p22'])
exp(fitnlminb$par['p33'])
exp(fitnlminb$par['e11'])
sqrt(exp(fitnlminb$par['e11']))
##-----------------------------------------------------------------
## EVALUATE MODELS
## -------------- -------------- -------------- --------------

## Plot the simulated output from the two estimated models
## ctsmr
(fitCtsmr <- insert_xm(fit, fitctsmr$xm))
simCtsmr <- simulate(fitCtsmr, newdata=ice_data)
## nlminb
(fitNlminb <- insert_xm(fit, fitnlminb$par))
simNlminb <- simulate(fitNlminb, newdata=ice_data)


plot(ice_data$t, ice_data$CompCap, type="l", ylab="Heat power (kW)")
lines(ice_data$t, simCtsmr$output$sim$CompCap, type="l", col=2)
lines(ice_data$t, simNlminb$output$sim$CompCap, type="l", col=3)

## Plot the stochastic simulated output from the two estimated models
simCtsmrStoch = stochastic.simulate(fitCtsmr, c(fitCtsmr$xm['X0']), data=ice_data, dt = 0.01)
simNlminStoch = stochastic.simulate(fitNlminb, c(fitNlminb$xm['X10'], fitNlminb$xm['X10']), data=ice_data, dt = 0.01)

plot(ice_data$CompCap, type = "b",ylim = c(0.4,1))
lines(simCtsmr$output$sim,type = "l", col = 'red')
lines(simCtsmr$output$sim + simCtsmr$output$sd,type = "l", col = 'red')
lines(simCtsmr$output$sim - simCtsmr$output$sd,type = "l", col = 'red')
lines(simCtsmrStoch$output$CompCap,type = "l", col = 'green')

plot(ice_data$CompCap, type = "b",ylim = c(0.4,1))
lines(simNlminb$output$sim,type = "l", col = 'red')
lines(simNlminStoch$output$CompCap,type = "l", col = 'green')

#### RESIDUAL ANALYSIS
res = ice_data$CompCap - simNlminb$output$sim$CompCap
plot(res, type = "l")
diagtool(res)

plot(ice_data$CompCap,res)

lags = seq(1,20)
ldf_vals = ldf(res,lags,nBoot=30,plotIt=FALSE,plotFits=TRUE)
val = ldf_vals$val
iidVal = ldf_vals$iidVal

plot(c(0,lags), c(1,val), type="n", ylim=c(-1,1), ylab="LDF", main="Lag Dependence Functions", xaxt="n", xlab="lag")
axis(1,c(0,lags))
abline(0,0,col="gray")
lines(c(0,lags), c(1,val), type="h")
## Draw the approximate 95% confidence interval
abline(h=quantile(iidVal,0.95), col="blue", lty=2)














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
plot(ice_data$t, ice_data$CompCap, type="l", ylab="Heat power (kW)")
lines(ice_data$t, simCtsmr$output$sim$CompCap, type="l", col=3)

simNlminb <- simulate(fitNlminb, newdata=ice_data, dt = 0.00001)
lines(ice_data$t, simNlminb$output$sim$CompCap, type="l", col='red')
#lines(ice_data$t, sim$output$sim$CompCap, type="l",col='red')

summary(fitNlminb)

## manual
(fitManual <- insert_xm(fit, fit_manual$xm))
simManual <- simulate(fitManual, newdata=ice_data)

simCtsmrStoch = stochastic.simulate(fitCtsmr, c(fitCtsmr$xm['T00'], fitCtsmr$xm['T10'], fitCtsmr$xm['T20']), data=ice_data, dt = 0.00001)
simNlminbStoch = stochastic.simulate(fitNlminb, c(fitNlminb$xm['T00'], fitNlminb$xm['T10'], fitNlminb$xm['T20']), data=ice_data, dt = 0.00001)

plot(simCtsmr$state$sim$T2,type = "l", col = 'red')
plot(simCtsmrStoch$state$T2,type = "l", col = 'green')

#plot(simNlminb$state$sim$T1,type = "l", col = 'red')
#plot(simNlminbStoch$state$T2,type = "l", col = 'green')

plot(ice_data$CompCap, type = "b",ylim = c(0.4,1))
lines(simCtsmr$output$sim,type = "l", col = 'red')
lines(simCtsmrStoch$output$CompCap,type = "l", col = 'green')

plot(ice_data$CompCap, type = "b",ylim = c(0.4,1))
lines(simNlminb$output$sim,type = "l", col = 'red')
lines(simNlminbStoch$output$CompCap,type = "l", col = 'green')
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

