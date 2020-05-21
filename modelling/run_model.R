setwd("~/Dropbox/Skole/DTU/Studie/MASTER/CODE/old/ice_tank_switching_times/model/peder")
library('ctsmr')
library(latex2exp)
library(boot)
source('../code/tools.R')
source("./yum_tools.R")



##-----------------------------------------------------------------
## Read data
dataFull = read.csv2('../data/high_tc_20190913.csv', fileEncoding="latin1")
N = dim(dataFull)[1]
set.seed(144257)


# Subsampling
nSkip = 50
Idx = seq(1,N,nSkip)
data = dataFull[Idx,]
n = dim(data)[1]


# Define provided status changes
changeTimes <- c(as.POSIXct("2019-09-17 13:00:00"),
                 as.POSIXct("2019-09-17 13:30:00"),
                 as.POSIXct("2019-09-17 14:00:00"),
                 as.POSIXct("2019-09-17 14:45:00"),
                 as.POSIXct("2019-09-17 15:30:00"),
                 as.POSIXct("2019-09-17 16:30:00"),
                 as.POSIXct("2019-09-17 17:30:00"),
                 as.POSIXct("2019-09-17 19:00:00"))


# Data
time = as.POSIXct(data$Time)
hour_index = which(c(1,diff(as.numeric(strftime(time, format="%H"))))==1)

# Convert to step function
status = createMeltSteps(changeTimes,time)
ice_data = data.frame("status" = status, "CompCap" = data$A40_A_CompCap/100)
ice_data$t = seq(1,n)*6

plot.ts(ice_data[,c("status","CompCap")])



##-----------------------------------------------------------------
## Plot the sigmoid function
Ti <- seq(0,1, by=0.001)
Tset <- 0.5
sigmoidSlope <- 4
plot(Ti, (1/(1+exp((Tset-Ti)*sigmoidSlope))),type = "l")


##-----------------------------------------------------------------
# MODEL 2
## -------------- -------------- -------------- --------------
model_dir = 'model_2_sigmoid'
source(paste0("./",model_dir,"/model.R"))
source(paste0("./",model_dir,"/setprm.R"))
#plotsim(simulate(makefit(setprm(model())), newdata=ice_data),ice_data)


sim = simulate(makefit(setprm(model())), newdata=ice_data)
plot(sim$state$sim$X0, type = "l")

par(mfrow = c(1,1))
plot(ice_data$CompCap, type = "b",ylim = c(0,1))
lines(sim$output$sim$CompCap, type = "l",col = "red")

 
plot(sim$state$sim$T1, type = "l")
plot(sim$state$sim$T2, type = "l",col = "red")
lines(sim$state$sim$T3, type = "b",col = "green")

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

simCtsmrStoch = stochastic.simulate(fitCtsmr, c(fitCtsmr$xm['X00'], fitCtsmr$xm['X10']), data=ice_data, dt = 0.00001)
simNlminbStoch = stochastic.simulate(fitNlminb, c(fitNlminb$xm['X00'], fitNlminb$xm['X10']), data=ice_data, dt = 0.00001)

plot(ice_data$CompCap, type = "b",ylim = c(0.4,1))
lines(simCtsmr$output$sim,type = "l", col = 'red')
lines(simCtsmrStoch$output$CompCap,type = "l", col = 'green')

plot(ice_data$CompCap, type = "b",ylim = c(0.4,1))
lines(simNlminb$output$sim,type = "l", col = 'red')
lines(simNlminb$state$sim$X0/1000,type = "l", col = 'blue')
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

##-----------------------------------------------------------------


##-----------------------------------------------------------------
# MODEL 3
## -------------- -------------- -------------- --------------
model_dir = 'model_3'
source(paste0("./",model_dir,"/model.R"))
source(paste0("./",model_dir,"/setprm.R"))
#plotsim(simulate(makefit(setprm(model())), newdata=ice_data),ice_data)


sim = simulate(makefit(setprm(model())), newdata=ice_data)
par(mfrow = c(1,1))
plot(ice_data$CompCap, type = "b",ylim=c(0,1))
lines(sim$output$sim$CompCap, type = "l",col = "red")
T2 = sim$state$sim$T2
T2

T0 = sim$state$sim$T0
T1 = sim$state$sim$T1
T2 = sim$state$sim$T2
(-4 * T2 - 4 * T1 - 1 * (T0+3.8))


plot(sim$state$sim$T1, type = "l")
plot(sim$state$sim$T2, type = "l",col = "red")
lines(sim$state$sim$T3, type = "b",col = "green")

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

# 864.6123 
#fitNlminb$xm['T10'] = 2.82
#fitNlminb$xm['T20'] = 864.6123 
plot(ice_data$t, ice_data$CompCap, type="l", ylab="Heat power (kW)")
simNlminb <- simulate(fitNlminb, newdata=ice_data, dt = 0.00001)
lines(ice_data$t, simNlminb$output$sim$CompCap, type="l", col=3)

summary(fitnlminb)

## manual
(fitManual <- insert_xm(fit, fit_manual$xm))
simManual <- simulate(fitManual, newdata=ice_data)

simCtsmrStoch = stochastic.simulate(fitCtsmr, c(fitCtsmr$xm['T00'], fitCtsmr$xm['T10'], fitCtsmr$xm['T20']), data=ice_data, dt = 0.00001)
simNlminbStoch = stochastic.simulate(fitNlminb, c(fitNlminb$xm['T00'], fitNlminb$xm['T10'], fitNlminb$xm['T20']), data=ice_data, dt = 0.00001)

plot(simCtsmr$state$sim$T2,type = "l", col = 'red')
plot(simCtsmrStoch$state$T2,type = "l", col = 'green')

plot(simNlminb$state$sim$T1,type = "l", col = 'red')
plot(simNlminbStoch$state$T2,type = "l", col = 'green')

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



##-----------------------------------------------------------------
# MODEL 4
## -------------- -------------- -------------- --------------
model_dir = 'model_4'
source(paste0("./",model_dir,"/model.R"))
source(paste0("./",model_dir,"/setprm.R"))
#plotsim(simulate(makefit(setprm(model())), newdata=ice_data),ice_data)


sim = simulate(makefit(setprm(model())), newdata=ice_data)
par(mfrow = c(1,1))
plot(ice_data$CompCap, type = "b",ylim=c(0,1))
lines(sim$output$sim$CompCap, type = "l",col = "red")
T2 = sim$state$sim$T2
T2

T0 = sim$state$sim$T0
T1 = sim$state$sim$T1
T2 = sim$state$sim$T2
(-4 * T2 - 4 * T1 - 1 * (T0+3.8))


plot(sim$state$sim$T1, type = "l")
plot(sim$state$sim$T2, type = "l",col = "red")
lines(sim$state$sim$T3, type = "b",col = "green")

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

# 864.6123 
#fitNlminb$xm['T10'] = 2.82
#fitNlminb$xm['T20'] = 864.6123 
plot(ice_data$t, ice_data$CompCap, type="l", ylab="Heat power (kW)")
simNlminb <- simulate(fitNlminb, newdata=ice_data, dt = 0.00001)
lines(ice_data$t, simNlminb$output$sim$CompCap, type="l", col=3)

summary(fitnlminb)

## manual
(fitManual <- insert_xm(fit, fit_manual$xm))
simManual <- simulate(fitManual, newdata=ice_data)

simCtsmrStoch = stochastic.simulate(fitCtsmr, c(fitCtsmr$xm['T00'], fitCtsmr$xm['T10'], fitCtsmr$xm['T20']), data=ice_data, dt = 0.00001)
simNlminbStoch = stochastic.simulate(fitNlminb, c(fitNlminb$xm['T00'], fitNlminb$xm['T10'], fitNlminb$xm['T20']), data=ice_data, dt = 0.00001)

plot(simCtsmr$state$sim$T2,type = "l", col = 'red')
plot(simCtsmrStoch$state$T2,type = "l", col = 'green')

plot(simNlminb$state$sim$T1,type = "l", col = 'red')
plot(simNlminbStoch$state$T2,type = "l", col = 'green')

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

