import numpy as np

# MATH
def sigmoid(x,slope,offset):
    return 1./(1. + np.exp(-(slope * (x - offset))))

def derive_regimes(switches,tf,endpoints):
    n_switches = switches.size
    n_cycles = int(n_switches/2)
    tau_MELT = switches[np.arange(0,n_switches,2)]
    tau_IDLE = switches[np.arange(1,n_switches,2)]
    
    if endpoints:
        tau_IDLE = np.insert(tau_IDLE,0,0)
        tau_MELT = np.append(tau_MELT,tf)
    
    return tau_MELT, tau_IDLE 

# ENERGY COST
def smooth_dap(t,dap,slope = 1):
    dat = 60 * np.arange(dap.size + 1); dat[0] = dat[0] - 60; dat[-1] = dat[-1] + 60
    
    _dap = 0
    for k in range(24):
            _dap += dap[k] / ((1. + np.exp(np.minimum(-slope * (t - dat[k]), 15.0 ))) *
                             (1. + np.exp( np.minimum(slope * (t - dat[k + 1]), 15.))))
    return _dap


def smooth_regime(T,switches,slope=10.):
    n_s = int(len(switches)/2)
    tau_MELT, tau_IDLE  = switches[:n_s] , switches[n_s:] #derive_regimes(switches,T[-1],0)
    regime = 0
    for k in range(len(tau_MELT)):
            regime += 1/ ((1 + np.exp(np.minimum(-slope* (T - tau_MELT[k]), 15.0 ))) *
                             (1 + np.exp( np.minimum(slope* (T - tau_IDLE[k]), 15.))))
            
    return regime

class firstordermodel: 
    def __init__(self, pars,regime_slope=20.,smooth=False):
        self.lambda_MELT = pars[0]
        self.mu_MELT = pars[1]
        self.lambda_IDLE = pars[2]
        self.mu_IDLE = pars[3]
        self.sigma = pars[4]
        self.R = pars[5]
        
        self.drift_pars = np.array([pars[i] for i in range(4)])

        self.pars = pars
        self.regime_slope = regime_slope
        self.smooth=smooth
        
        
    
        # Consider adding noise term
        
    def f_MELT(self, t,x):
        x = np.array(x)
        return(self.lambda_MELT * (self.mu_MELT - x))

    def f_IDLE(self, t,x):
        x = np.array(x)
        return(self.lambda_IDLE * (self.mu_IDLE - x))

    def f(self,t,x,switches):
        regime = smooth_regime(t,switches,self.regime_slope)

        if self.smooth:
            f = regime * self.f_MELT(t,x) + (1-regime) * self.f_IDLE(t,x)
        else:
            if regime > 0.5:
                f = self.f_MELT(t,x)
            else:
                f = self.f_IDLE(t,x)
        
        return(f)
    
    def f_wrt_x(self,t,x,switches):
        regime = smooth_regime(t,switches)
        
        fdx = regime * (-self.lambda_IDLE) + (1-regime) * -self.lambda_MELT
        #for i in range(tau_MELT.size-1):
        #    if tau_MELT[i] <= t and t < tau_IDLE[i+1]: # Check if MELT
        #        return(-self.lambda_IDLE)

        #    elif tau_IDLE[i] <= t and t < tau_MELT[i]: # Check if IDLE
        #        return(-self.lambda_IDLE)

        # If last interval
        return(fdx)
    
    def h(self,x):
        y = np.array(x)
        return y
    
    
    
class secondordermodel: 
    def __init__(self, pars,regime_slope=20.,smooth=False):
        self.mu_IDLE = pars[0]
        self.mu_MELT = pars[1]
        self.omega_IDLE = pars[2]
        self.omega_MELT = pars[3]
        self.xi_IDLE = pars[4]
        self.xi_MELT = pars[5]
        self.slope = pars[6]
        self.offset = pars[7]
        self.sigma = pars[8]
        self.R = pars[9]

        self.drift_pars = np.array([pars[i] for i in range(8)])

        self.pars = pars
        self.regime_slope = regime_slope
        self.smooth=smooth


    def f_MELT(self,t,x):
        dxdt = np.zeros(2)
        dxdt[0] = x[1]
        dxdt[1] = -2 * self.xi_MELT * self.omega_MELT * (x[1]) - self.omega_MELT**2 * (x[0]-self.mu_MELT)

        return dxdt

    def f_IDLE(self,t,x):
        dxdt = np.zeros(2)
        dxdt[0] = x[1]
        dxdt[1] = -2 * self.xi_IDLE * self.omega_IDLE * (x[1]) - self.omega_IDLE**2 * (x[0]-self.mu_IDLE)

        return dxdt


    def f(self,t,x,switches):
        regime = smooth_regime(t,switches,self.regime_slope)

        if self.smooth:
            f = regime * self.f_MELT(t,x) + (1-regime) * self.f_IDLE(t,x)
        else:
            if regime > 0.5:
                f = self.f_MELT(t,x)
            else:
                f = self.f_IDLE(t,x)
        
        return(f)
    
    def h(self,x):
        #y = x[1]
        #y = np.expand_dims(sigmoid(x[1],0.00778941,599.12794983), axis = 0)
        y = np.expand_dims(sigmoid(x[0],self.slope,self.offset), axis = 0)
        return y
    
    
class thirdordermodel: 
    def __init__(self, pars,regime_slope=20.,smooth=False):
        self.mu_IDLE = pars[0]
        self.mu_MELT = pars[1]
        self.a0_IDLE = pars[2]
        self.a0_MELT = pars[3]
        self.a1_IDLE = pars[4]
        self.a1_MELT = pars[5]
        self.a2_IDLE = pars[6]
        self.a2_MELT = pars[7]
        self.slow = pars[8]
        self.offset = pars[9]
        self.slope = pars[10]
        self.sigma = pars[11]
        self.R = pars[12]
        
        self.drift_pars = np.array([pars[i] for i in range(11)])

        self.pars = pars
        self.regime_slope = regime_slope
        self.smooth=smooth

    def f_MELT(self,t,x):
        dxdt = np.zeros(3)
        dxdt[0] = x[1]/self.slow
        dxdt[1] = x[2]/self.slow
        dxdt[2] = 1/self.slow * (-self.a2_MELT * x[2] - self.a1_MELT * x[1] - self.a0_MELT * (x[0] - self.mu_MELT) )

        return dxdt

    def f_IDLE(self,t,x):
        dxdt = np.zeros(3)
        dxdt[0] = x[1]/self.slow
        dxdt[1] = x[2]/self.slow
        dxdt[2] = 1/self.slow * (-self.a2_IDLE * x[2] -self.a1_IDLE * x[1] - self.a0_IDLE * (x[0] - self.mu_IDLE ))

        return dxdt


    def f(self,t,x,switches):
        regime = smooth_regime(t,switches,self.regime_slope)

        if self.smooth:
            f = regime * self.f_MELT(t,x) + (1-regime) * self.f_IDLE(t,x)
        else:
            if regime > 0.5:
                f = self.f_MELT(t,x)
            else:
                f = self.f_IDLE(t,x)
        
        return(f)
    
    def h(self,x):
        #y = x[1]
        #y = np.expand_dims(sigmoid(x[1],0.00778941,599.12794983), axis = 0)
        y = np.expand_dims(1./(1. + np.exp(-(self.slope * (x[0] - self.offset)))) , axis = 0)
        return y