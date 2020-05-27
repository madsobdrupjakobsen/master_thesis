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




# Second order - first
mu_IDLE = 689.980557
mu_MELT = 579.448947
omega_IDLE = 0.050398
omega_MELT = 0.196751 
xi_IDLE = 0.087644
xi_MELT = 0.596060 
slope = 0.016102
offset = 528.908275
logsigma = np.array([1.990430 , -999999])
logR =-9.921732

pars = np.array([mu_IDLE, mu_MELT, omega_IDLE, omega_MELT, xi_IDLE, xi_MELT, slope, offset, np.exp(logsigma), np.exp(logR)])
m2_first = secondordermodel(pars)
x0_2_first = np.array([664.112676, -1.143137])

# Second order - second ctsmr
mu_IDLE = 695.444864
mu_MELT = 573.927181
omega_IDLE = 0.062476
omega_MELT = 0.206864
xi_IDLE = 0.502492
xi_MELT =  0.516992 
slope = 0.018510 
offset = 527.489203
logsigma = np.array([-99999, 0.353187])
logR =-10.224420

pars = np.array([mu_IDLE, mu_MELT, omega_IDLE, omega_MELT, xi_IDLE, xi_MELT, slope, offset, np.exp(logsigma), np.exp(logR)])
m2_second_ctsmr = secondordermodel(pars)
x0_2_second_ctsmr= np.array([645.241334, 2.469267])

# Second order - second nlminb
mu_IDLE = 981.1224368
mu_MELT = 701.2119410
omega_IDLE =  0.0624759 
omega_MELT = 0.2068639 
xi_IDLE = 0.5024971 
xi_MELT = 0.5169938
slope = 0.0080359
offset =594.2449282
logsigma = np.array([-99999, 1.1875941])
logR =-10.2244098

pars = np.array([mu_IDLE, mu_MELT, omega_IDLE, omega_MELT, xi_IDLE, xi_MELT, slope, offset, np.exp(logsigma), np.exp(logR)])
m2_second_nlminb = secondordermodel(pars)
x0_2_second_nlminb = np.array([865.4797506, 5.6877789])

# Second order model - used
m2 = m2_both
x0_2 = x0_2_both

'''
T00      0.764273         NA      NA       NA
T10      0.065197         NA      NA       NA
T20     -0.132599         NA      NA       NA
a0_0     0.137103         NA      NA       NA
a0_1    11.053499         NA      NA       NA
a1_0     0.294471         NA      NA       NA
a1_1     8.627593         NA      NA       NA
a2_0     1.243450         NA      NA       NA
a2_1     3.532068         NA      NA       NA
e11    -10.977195         NA      NA       NA
mu0      0.978603         NA      NA       NA
mu1      0.374562         NA      NA       NA
offset   0.152946         NA      NA       NA
p11     -3.635138         NA      NA       NA
p22    -11.418277         NA      NA       NA
p33     -4.228776         NA      NA       NA
slope    3.468666         NA      NA       NA
slow     4.909405         NA      NA       NA
'''

# Third order model
'''
mu_IDLE = 0.929224978
mu_MELT = 0.683920449 
a0_IDLE = 0.023590616
a0_MELT = 0.019370754
a1_IDLE = 0.867249723  
a1_MELT = 0.151822292
a2_IDLE = 46.261177791   
a2_MELT = 1.206710500
slow = 0.318488737 

logsigma = np.array([-5.110433517, -18.493395485,  -5.116712113])
logR = -8.444535962
'''

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
'''
a0_0     0.0167358   0.0061738   2.7108 0.0075726 ** 
a0_1     0.0463627   0.0220815   2.0996 0.0375916 *  
a1_0     0.3411097   0.1697582   2.0094 0.0464550 *  
a1_1     0.3264330   0.1642066   1.9879 0.0488058 *  
a2_0    49.7535156   0.1904151 261.2898 < 2.2e-16 ***
a2_1     2.4016786   1.3464798   1.7837 0.0766826 .  
e11    -10.8960154   0.3771985 -28.8867 < 2.2e-16 ***
mu0      0.6183479   0.2160731   2.8618 0.0048762 ** 
mu1      0.2234170   0.0585640   3.8149 0.0002062 ***
offset   0.0642175   0.0114054   5.6304 1.024e-07 ***
p11     -4.2524933   0.4268121  -9.9634 < 2.2e-16 ***
p22     -7.9162707   0.7001128 -11.3071 < 2.2e-16 ***
p33     -4.4415517   0.6851678  -6.4824 1.666e-09 ***
slope    4.9465409   2.0322216   2.4341 0.0162176 *  
slow     0.2836255   0.0423370   6.6992 5.603e-10 ***
'''


# pars from jupyter
""" mu_IDLE = 0.929224978
mu_MELT = 0.683920449 
a0_IDLE = 0.023590616
a0_MELT = 0.019370754
a1_IDLE = 0.867249723  
a1_MELT = 0.151822292
a2_IDLE = 40. #46.261177791   
a2_MELT = 1.206710500
slow = 0.318488737 
logsigma = np.array([-5.110433517, -18.493395485,  -5.116712113])
logR = -8.444535962 """