# Set the working directory to the location of this file
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
ldf_vals = ldf(res,lags,nBoot=30,plotIt=FALSE,plotFits=FALSE)
val = ldf_vals$val
iidVal = ldf_vals$iidVal

plot(c(0,lags), c(1,val), type="n", ylim=c(-1,1), ylab="LDF", main="Lag Dependence Functions", xaxt="n", xlab="lag")
axis(1,c(0,lags))
abline(0,0,col="gray")
lines(c(0,lags), c(1,val), type="h")
## Draw the approximate 95% confidence interval
abline(h=quantile(iidVal,0.95), col="blue", lty=2)
