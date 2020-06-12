import numpy as np
from models import firstordermodel, secondordermodel, secondordermodel_scaled, thirdordermodel

# First order model
lambda_MELT = 0.077234
lambda_IDLE = 0.044628
mu_MELT = 0.654915
mu_IDLE = 0.947308
logsigma = -4.444617
logR = -5.979276
pars = np.array([lambda_MELT, mu_MELT, lambda_IDLE, mu_IDLE, np.exp(logsigma), np.exp(logR)])
m1 = firstordermodel(pars)
x0_1 = np.array([0.879327])

# Second order - both - NLMINB
mu_IDLE = 961.9963477
mu_MELT = 698.8612058
omega_IDLE = 0.0543021
omega_MELT = 0.1996318 
xi_IDLE = 0.1317791 
xi_MELT = 0.6136681
slope = 0.0073091
offset = 592.0101235 
logsigma = np.array([2.7207129 , -0.5016485])
logR =-9.9359985

pars = np.array([mu_IDLE, mu_MELT, omega_IDLE, omega_MELT, xi_IDLE, xi_MELT, slope, offset, np.exp(logsigma), np.exp(logR)])
m2_both = secondordermodel(pars)
x0_2_both = np.array([889.1909001, -2.0686304])

# Second order - both - CTSMR
mu_IDLE = 692.80388228
mu_MELT = 571.21642534
omega_IDLE = 0.06247610 
omega_MELT =0.20686358 
xi_IDLE = 0.50249192
xi_MELT = 0.51699198
slope = 0.01849967
offset = 524.75178335
logsigma = np.array([-7.52515120 , 0.35376150])
logR =-10.22442105

pars = np.array([mu_IDLE, mu_MELT, omega_IDLE, omega_MELT, xi_IDLE, xi_MELT, slope, offset, np.exp(logsigma), np.exp(logR)])
m2_both_ctsmr = secondordermodel(pars)
x0_2_both_ctsmr = np.array([642.57154115, 2.47069021])


# Second order - NLMIB - scaled
mu_IDLE = 0.999999999
mu_MELT = 1.000000000
omega_IDLE = 0.054302066
omega_MELT = 0.199631822
xi_IDLE = 0.131779060
xi_MELT = 0.613668062 
slope = 0.007308264
offset = 1.000000013
logsigma = np.array([2.720712912 * (1-np.log(889.190900140)/2.720712912) , -0.501648509])
logR =-9.935998459

pars = np.array([mu_IDLE, mu_MELT, omega_IDLE, omega_MELT, xi_IDLE, xi_MELT, slope, offset, np.exp(logsigma), np.exp(logR)])
m2_both_scaled = secondordermodel_scaled(pars)
x0_2_both_scaled = np.array([1.001124600, -2.068630356])




# Second order model - used
m2 = m2_both_scaled
x0_2 = x0_2_both_scaled


# Third order model

# Define model - nlminb
mu_IDLE = 0.978603
mu_MELT = 0.374562 
a0_IDLE = 0.137103
a0_MELT = 11.053499
a1_IDLE = 0.294471
a1_MELT = 8.627593
a2_IDLE = 1.243450
a2_MELT = 3.532068
slow = 4.909405
slope = 3.468666
offset = 0.152946

logsigma = np.array([-3.635138, -11.418277,  -4.228776])
logR = -10.977195

pars = np.array([mu_IDLE, mu_MELT, a0_IDLE, a0_MELT, a1_IDLE, a1_MELT, a2_IDLE, a2_MELT, slow, offset, slope, np.exp(logsigma), np.exp(logR)])
m3_nlminb = thirdordermodel(pars)
x0_3_nlminb = np.array([0.764273,  0.065197,  -0.132599])

# Same model butwith 0 diffusion
sigma = np.array([0.,0.,0.])
pars = np.array([mu_IDLE, mu_MELT, a0_IDLE, a0_MELT, a1_IDLE, a1_MELT, a2_IDLE, a2_MELT, slow, offset, slope, sigma, np.exp(logR)])
m3_zero_diff = thirdordermodel(pars)

# Define model - ctsmr
mu_IDLE = 0.6183479
mu_MELT = 0.2234170
a0_IDLE = 0.0167358
a0_MELT = 0.0463627
a1_IDLE = 0.3411097
a1_MELT = 0.3264330
a2_IDLE = 49.7535156
a2_MELT = 2.4016786
slow = 0.2836255
slope = 4.9465409
offset = 0.0642175

logsigma = np.array([-4.2524933, -7.9162707,  -4.4415517])
logR = -10.977195

pars = np.array([mu_IDLE, mu_MELT, a0_IDLE, a0_MELT, a1_IDLE, a1_MELT, a2_IDLE, a2_MELT, slow, offset, slope, np.exp(logsigma), np.exp(logR)])
m3_ctsmr = thirdordermodel(pars)
x0_3_ctsmr = np.array([0.4914763,  0.1570753,  -7.8600795])


m3 = m3_nlminb
x0_3 = x0_3_nlminb
