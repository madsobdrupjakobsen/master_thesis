setwd("~/Dropbox/Skole/DTU/Studie/MASTER/master_thesis/modelling")
library('ctsmr')
library(latex2exp)
library(boot)
source('./src/tools.R')
source("./src/yum_tools.R")
source('./src/prep.R')

##-----------------------------------------------------------------
# MODEL 1
## -------------- -------------- -------------- --------------
model_dir = 'model_1'
source(paste0("./",model_dir,"/model.R"))
source(paste0("./",model_dir,"/setprm.R"))

# Simulate based on initial parameter values
sim = simulate(makefit(setprm(model())), newdata=ice_data)
par(mfrow = c(1,1))
plot(ice_data$CompCap, type = "b",ylim = c(0,1))
lines(sim$output$sim$CompCap, type = "l",col = "red")

# Construct a model with initial values
mod_man = model()
setprm(mod_man)
fit_manual <- list()
fit_manual$model <- mod_man
fit_manual$xm <- mod_man$ParameterValues$initial
names(fit_manual$xm) = row.names(mod_man$ParameterValues)


##-----------------------------------------------------------------
# Prepare model for nlminb
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
##-----------------------------------------------------------------

## The likelihood of the initial parameters set in fit
nllikelihood(fit$xm[model$pars], fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,printit = FALSE)


##-----------------------------------------------------------------
## ESTIMATE PARAMETERS
## -------------- -------------- -------------- --------------
## Fit with different optimizers

## Fit with ctsmr
fitctsmr <- model$estimate(ice_data, firstorder=TRUE)

## Fit with nlminb
fitnlminb <- nlminb(fit$xm[model$pars], nllikelihood, fit=fit, D=ice_data, firstorder=TRUE, c=3, n.ahead=1,
                    control=list(iter.max=400), lower = fit$lower[model$pars], upper = fit$upper[model$pars])


## Plot the simulated output from the two estimated models
## ctsmr
(fitCtsmr <- insert_xm(fit, fitctsmr$xm))
simCtsmr <- simulate(fitCtsmr, newdata=ice_data)
## nlminb
(fitNlminb <- insert_xm(fit, fitnlminb$par))
simNlminb <- simulate(fitCtsmr, newdata=ice_data)
## manual
(fitManual <- insert_xm(fit, fit_manual$xm))
simManual <- simulate(fitManual, newdata=ice_data)




##-----------------------------------------------------------------
## EVALUATE MODELS
## -------------- -------------- -------------- --------------
plot(ice_data$t, ice_data$CompCap, type="l", ylab="Heat power (kW)")
lines(ice_data$t, simCtsmr$output$sim$CompCap, type="l", col=2)
lines(ice_data$t, simNlminb$output$sim$CompCap, type="l", col=3)


## Plot the stochastic simulated output from the two estimated models
simCtsmrStoch = stochastic.simulate(fitCtsmr, c(fitCtsmr$xm['X0']), data=ice_data, dt = 0.01)
simNlminStoch = stochastic.simulate(fitNlminb, c(fitNlminb$xm['X0']), data=ice_data, dt = 0.01)

plot(ice_data$CompCap, type = "b",ylim = c(0.4,1))
lines(simCtsmr$output$sim,type = "l", col = 'red')
lines(simCtsmr$output$sim + simCtsmr$output$sd,type = "l", col = 'red')
lines(simCtsmr$output$sim - simCtsmr$output$sd,type = "l", col = 'red')
lines(simCtsmrStoch$output$CompCap,type = "l", col = 'green')

plot(ice_data$CompCap, type = "b",ylim = c(0.4,1))
lines(simNlminb$output$sim,type = "l", col = 'red')
lines(simNlminStoch$output$CompCap,type = "l", col = 'green')
##-----------------------------------------------------------------

#### ESTIMATE UNCERTAINTY WITH FINITE DIFFERENCE
summary(fitctsmr)

hessObjNlminb = finiteDiffHessian2(fitnlminb$par,r = 15)
hessObjCtsmr = finiteDiffHessian2(fitctsmr$xm,r=15)

cbind(fitnlminb$par, hessObjNlminb$sd, hessObjNlminb$pvals)

# Confidence intervals - ctmsr
c(exp(fitctsmr$xm['p11']),exp(fitctsmr$xm['p11'] +  c(-1,1) * fitctsmr$sd['p11']))
c(exp(fitctsmr$xm['e11']),exp(fitctsmr$xm['e11'] +  c(-1,1) * fitctsmr$sd['e11']))
sqrt(c(exp(fitctsmr$xm['e11']),exp(fitctsmr$xm['e11'] +  c(-1,1) * fitctsmr$sd['e11'])))

# Confidence intervals - nlminb
c(exp(fitnlminb$par['p11']) ,exp(fitnlminb$par['p11'] +  c(-1,1) * hessObjNlminb$sd['p11']))
c(exp(fitnlminb$par['e11']) ,exp(fitnlminb$par['e11'] +  c(-1,1) * hessObjNlminb$sd['e11']))
sqrt(c(exp(fitnlminb$par['e11']) ,exp(fitnlminb$par['e11'] +  c(-1,1) * hessObjNlminb$sd['e11'])))


#### RESIDUAL ANALYSIS
res = ice_data$CompCap - simNlminb$output$sim$CompCap
plot(res)
diagtool(res)


