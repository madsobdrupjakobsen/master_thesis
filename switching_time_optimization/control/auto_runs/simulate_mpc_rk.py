#%%

import os
import numpy as np
import random
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patch
from scipy.integrate import solve_ivp
import scipy.optimize as sc_opt
import math
from scipy.integrate import quad
import time
from IPython.core.debugger import set_trace
import datetime as dt

from numpy import genfromtxt

import sys
sys.path.append('../../models')
from models import firstordermodel, secondordermodel, thirdordermodel

sys.path.append('../')
from tools import stochasticSimulation, derive_regimes, discrete_ivp_solver, \
                    smooth_dap, sol_ivp_wrapper, sol_ivp_wrapper_discrete,\
                    smooth_regime, cost, sigmoid,\
                    simulate_MPC, mergeSwitchCont,\
                    build_initial_ipopt_object,\
                    consMatrix, f_with_consumption, plotSwitches, f_with_cost, \
                    consumption, stochasticSimulation_IDLE, removeRedundantSwitches



#%load_ext autoreload
#%autoreload 2

# PATHS
#FIGS = '/Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/THESIS/FIGS'

import os
import numpy as np
import random
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patch
from scipy.integrate import solve_ivp
import scipy.optimize as sc_opt
import math
from scipy.integrate import quad
import time
from IPython.core.debugger import set_trace
import datetime as dt

from numpy import genfromtxt


os.system("../model_loader.py")
os.system("../price_loader.py")
os.system("../../models/models.py")

import sys
sys.path.append('../../models')
from models import firstordermodel, secondordermodel, thirdordermodel


sys.path.append('../')
from tools import stochasticSimulation, derive_regimes, discrete_ivp_solver, \
                    smooth_dap, sol_ivp_wrapper, sol_ivp_wrapper_discrete,\
                    smooth_regime, cost, sigmoid,\
                    simulate_MPC, mergeSwitchCont,\
                    build_initial_ipopt_object,\
                    consMatrix, f_with_consumption, plotSwitches, f_with_cost, \
                    consumption, stochasticSimulation_IDLE
from model_loader import *
from price_loader import *

# PATHS
#FIGS = '/Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/THESIS/FIGS'

# Switching time optimization modules
import switching_times_1st as st1
import switching_times_2nd as st2
import switching_times_3rd as st3

import switching_times_1st_rk as st1_rk
import switching_times_2nd_rk as st2_rk
import switching_times_3rd_rk as st3_rk
import lxml


k_baseline = 12400.
k_MELT = 250.
k_IDLE = 9.

# Define model
#mu_IDLE = 0.929224978
#mu_MELT = 0.683920449 
#a0_IDLE = 0.023590616
#a0_MELT = 0.019370754
#a1_IDLE = 0.867249723  
#a1_MELT = 0.151822292
#a2_IDLE = 46.261177791   
#a2_MELT = 1.206710500
#slow = 0.318488737 
#logsigma = np.array([-5.110433517, -18.493395485,  -5.116712113])
#logR = -8.444535962

#pars = np.array([mu_IDLE, mu_MELT, a0_IDLE, a0_MELT, a1_IDLE, a1_MELT, a2_IDLE, a2_MELT, slow, logsigma, logR])
#m3 = thirdordermodel(pars)


model_sys = int(sys.argv[1])
n_days = int(sys.argv[2])
n_s = int(sys.argv[3])
true_sys = int(sys.argv[4])
price_slope = float(sys.argv[5])
regime_slope = float(sys.argv[6])

logfile = open("../results/current_run_rk_"+str(model_sys)+".txt","w")
logfile.write("CURRENT RUN""\n")
logfile.write("--------------""\n")
logfile.write("Model: " + str(model_sys) +  "\n")
logfile.write("Number of days: " + str(n_days) +  "\n")
logfile.write("n_s: " + str(n_s) +  "\n")
logfile.write("True: " + str(true_sys) +  "\n")
logfile.write("--------------""\n")
logfile.close()


# BUILD OPTIMIZATION FUNCTION
# Define parameters
#n_s = 6
max_melt = 16. * 60.
dt_opt = 0.1

# ------ Define models -------
# Expectation and optimization object model
# -------------
if model_sys == 1:
    print(1)
    system_model = m1; x0_model = x0_1; tank = st1_rk.plant(); tank_pars = m1.drift_pars; sys_mod = 'm1'
elif model_sys == 2:
    print(2)
    system_model = m2; x0_model = x0_2; tank = st2_rk.plant(); tank_pars = m2.drift_pars; sys_mod = 'm2'
elif model_sys == 3:
    print(3)
    system_model = m3; x0_model = x0_3; tank = st3_rk.plant(); tank_pars = m3.drift_pars; sys_mod = 'm3'
# -------------

# "True" model
# -------------
if true_sys == 1:
    print(1)
    system_true = m1; x0_true = x0_1; sys_true = 'm1'
elif true_sys == 2:
    print(2)
    system_true = m2; x0_true = x0_2; sys_true = 'm2'
elif true_sys == 3:
    print(3)
    system_true = m3; x0_true = x0_3; sys_true = 'm3'
elif true_sys == 4: # With diffusion 0
    print(4)
    system_true = m3_zero_diff; x0_true = x0_3; sys_true = 'm4'
# -------------


#%%

#%%
# Build object
t0 = 0.
tf_sim = 60 * 24
tf_ph = 60 * 48
opt_rk = True
tank = build_initial_ipopt_object(tank, tank_pars, dt_opt, k_baseline, k_MELT, k_IDLE, t0, tf_ph, max_melt, n_s, price_slope, regime_slope, opt_rk)


# Build initial switches
# Build initial values - Other methods could be considered
idle = tf_ph * np.linspace(0.1,0.9,n_s) # Spread out on the whole interval
melt = idle - max_melt/n_s * 0.866 # Assign melt period to a little before idle
switch0_dap = np.concatenate((melt,idle)) # put together

idle = tf_ph * np.linspace(0.15,0.95,n_s) # Spread out on the whole interval
melt = idle - max_melt/n_s * 0.9 # Assign melt period to a little before idle
switch0_rk = np.concatenate((melt,idle)) # put together

switch0 = np.concatenate((switch0_dap,switch0_rk))
switch0

tank.set_p_optimize(switch0)



#%%
# SIMUALTE
#n_days = 2
n_skip = 2 #int(1/ och_sim)
start_date = '2019-01-01 12:00:00'
dt_stoch_sim = 0.1
seed = 1235

# Derive convenient parameters
np.random.seed(1234)
n_per_day =  int(1/n_skip * 1/dt_stoch_sim *60 * 24)
t0_sim = t0
t0_ph = t0
total_sim_time = tf_sim - t0_sim

# Initialize processes
nx_true = x0_true.size
nx_model = x0_model.size

# Initialzie history
history = {}
history['X_true'] = np.zeros([n_days, nx_true, n_per_day])
history['X_model_promised'] = np.zeros([n_days, nx_model, n_per_day])
history['X_model_rk'] = np.zeros([n_days, nx_model, n_per_day])
history['X_true_init_promised'] = np.zeros([n_days, nx_true, n_per_day])
history['X_true_init_rk'] = np.zeros([n_days, nx_true, n_per_day])
history['X_true_idle'] = np.zeros([n_days, nx_true, n_per_day])

history['Z_true'] = np.zeros([n_days, 1, n_per_day])
history['Z_model_promised'] = np.zeros([n_days, 1, n_per_day])
history['Z_model_rk'] = np.zeros([n_days, 1, n_per_day])
history['Z_true_init_promised'] = np.zeros([n_days, 1, n_per_day])
history['Z_true_init_rk'] = np.zeros([n_days, 1, n_per_day])
history['Z_true_idle'] = np.zeros([n_days, 1, n_per_day])


history['T'] = np.zeros([n_days, n_per_day])

history['SWITCHES_promised'] = np.zeros([n_days, 2 * n_s])
history['SWITCHES_rk'] = np.zeros([n_days, 2 * n_s])
history['SWITCHES_init_promised'] = np.zeros([n_days, 2 * n_s])
history['SWITCHES_init_rk'] = np.zeros([n_days, 2 * n_s])

history['PRICES_dap'] = np.zeros([n_days, 48])
history['PRICES_down'] = np.zeros([n_days, 48])
history['PRICES_up'] = np.zeros([n_days, 48])
history['PRICES_rk'] = np.zeros([n_days, 48])

history['price_model_promised'] = np.zeros([n_days,n_per_day])
history['price_model_rk'] = np.zeros([n_days,n_per_day])
history['price_true'] = np.zeros([n_days,n_per_day])
history['price_true_init_switch'] = np.zeros([n_days,n_per_day])
history['price_true_idle'] = np.zeros([n_days,n_per_day])

history['price_model_promised_each_hour'] = np.zeros([n_days,n_per_day])
history['price_model_rk_each_hour'] = np.zeros([n_days,n_per_day])
history['price_true_each_hour'] = np.zeros([n_days,n_per_day])
history['price_true_each_hour_init_switch'] = np.zeros([n_days,n_per_day])
history['price_true_each_hour_idle'] = np.zeros([n_days,n_per_day])

history['cons_model_promised'] = np.zeros([n_days,n_per_day])
history['cons_model_rk'] = np.zeros([n_days,n_per_day])
history['cons_true'] = np.zeros([n_days,n_per_day])
history['cons_true_init_promised'] = np.zeros([n_days,n_per_day])
history['cons_true_init_rk'] = np.zeros([n_days,n_per_day])
history['cons_true_idle'] = np.zeros([n_days,n_per_day])

history['cons_model_promised_each_hour'] = np.zeros([n_days,n_per_day])
history['cons_model_rk_each_hour'] = np.zeros([n_days,n_per_day])
history['cons_true_each_hour'] = np.zeros([n_days,n_per_day])
history['cons_true_each_hour_init_promised'] = np.zeros([n_days,n_per_day])
history['cons_true_each_hour_init_rk'] = np.zeros([n_days,n_per_day])
history['cons_true_each_hour_idle'] = np.zeros([n_days,n_per_day])

history['start_date'] = prices.index[0]

history['opt_time'] = np.zeros([n_days,1])
history['sim_time'] = np.zeros([n_days,1])


price_true_lasthour = 0.
price_true_lasthour_init_switch = 0.
price_true_lasthour_idle = 0.
price_promised_lasthour = 0.
price_rk_lasthour = 0.


cons_true_lasthour = 0.
cons_true_lasthour_init_promised = 0.
cons_true_lasthour_init_rk = 0.
cons_true_lasthour_idle = 0.
cons_promised_lasthour = 0.
cons_rk_lasthour = 0.

start_true = np.array(x0_true)
start_true_init_promised = np.array(x0_true)
start_true_init_rk = np.array(x0_true)
start_true_idle = np.array(x0_true)
start_promised = np.array(x0_model)
start_expected = np.array(x0_model)

idx_offset = np.where(prices.index == start_date)[0][0]


for day in range(n_days):
    print('  ----- Simulating day ' + str(day) + ' ----  ')
    logfile = open("../results/current_run_rk_"+str(model_sys)+".txt","a")
    logfile.write("Simulating day " + str(day) + '\n')
    logfile.write("--------------\n")
    logfile.close()
    
    t0_sub = t0_sim/24.
    tf_sub = tf_sim/24.
    
    # Extract prices for optimization
    idx_offset = np.where(prices.index == start_date)[0][0]
    idx = np.arange(48) + day * 24 + idx_offset
    q_dap = np.array(prices['spot'][idx]).astype(np.float) # * 1/1000000
    q_down = np.array(prices['down'][idx]).astype(np.float) # * 1/1000000
    q_up = np.array(prices['up'][idx]).astype(np.float) # * 1/1000000
    q_rk = np.array(prices['RK'][idx]).astype(np.float) # * 1/1000000
    
    # Insert new prices in optimization object
    p_dynamic = tank.get_p_dynamic()
    for k in range(0, 48):
        p_dynamic[k] = q_dap[k]
        p_dynamic[48 + k] = q_rk[k]
        tank.set_p_dynamic(p_dynamic)
        
    # Initial state of the system
    start_promised = start_expected # Should be estimated with Kalman Filter
    start_expected = start_expected # Should be estimated with Kalman Filter
    tank.set_x0(np.append(start_promised,np.array([0,0]))) 

    # Find optimal switches
    init_switch = tank.get_p_optimize()
    init_promised = init_switch[:2*n_s]
    init_rk = init_switch[2*n_s:]
    print(' ... Started optimizing')
    time_start = time.time()
    tank.solve()
    history['opt_time'][day] = time.time() - time_start
    print(' ... Finished optimizing')
    opt_switch = tank.get_p_optimize_ipopt()
    switch_promised = opt_switch[:2*n_s]
    switch_used = opt_switch[2*n_s:]
    
    
    #switch_opt_dap = np.sort(np.random.uniform(0,tf_ph,2*n_s))
    #switch_promised = np.concatenate(derive_regimes(switch_opt_dap,0,0))
    #print(switch_opt_dap)

    #switch_opt_rk = np.sort(np.random.uniform(0,tf_ph,2*n_s))
    #switch_used = np.concatenate(derive_regimes(switch_opt_rk,0,0))
    #print(switch_opt_rk)
    
    
    # Simulate process and compute prices
    print(' ... Simulating processes')
    time_start = time.time()
    for hour in range(24):

        # ---------------------------
        # Simulate 1 day
        # ---------------------------
        # Remove switches closer to each other than 10 minutes
        switch_used_filtered = removeRedundantSwitches(switch_used,10)
        switch_promised_filtered = removeRedundantSwitches(switch_promised,10)

        # Solve SDEs

        # True - Optimal RK switches
        if (switch_used_filtered.size > 0):
            T_tmp, X_true, Y_true, Z_true = stochasticSimulation(system_true,switch_used_filtered,start_true,t0_sub,tf_sub,dt_stoch_sim)
        else:
            T_tmp, X_true, Y_true, Z_true = stochasticSimulation_IDLE(system_true,start_true,t0_sub,tf_sub,dt_stoch_sim)

        # Expected Promised  - Optimal switches
        T_tmp1, X_promised, Z_promised = sol_ivp_wrapper_discrete(system_model,start_promised,switch_promised_filtered,t0_sub,tf_sub,T_tmp)

        # Exptected RK  - Optimal switches
        T_tmp2, X_rk_expected, Z_rk_expected = sol_ivp_wrapper_discrete(system_model,start_expected,switch_used_filtered,t0_sub,tf_sub,T_tmp)

        # True - Initial switches - Promised
        T_tmp_init_switch, X_true_init_switch_promised, Y_true_init_switch_promised, Z_true_init_switch_promised = stochasticSimulation(system_true,init_promised,start_true_init_promised,t0_sub,tf_sub,dt_stoch_sim)
        
        # True - Initial switches - RK
        T_tmp_init_switch, X_true_init_switch_rk, Y_true_init_switch_rk, Z_true_init_switch_rk = stochasticSimulation(system_true,init_rk,start_true_init_rk,t0_sub,tf_sub,dt_stoch_sim)
        
        # True - Always idle
        T_tmp_idle, X_true_idle, Y_true_idle, Z_true_idle = stochasticSimulation_IDLE(system_true,start_true_idle,t0_sub,tf_sub,dt_stoch_sim)


        # -----------------------------
        # Compute total consumption
        # -----------------------------
        # True
        cons_momental_true, cons_acc_true = consumption(Z_true[0],T_tmp,switch_used,k_baseline,k_MELT,k_IDLE)
        
        # Promised
        cons_momental_promised, cons_acc_promised = consumption(Z_promised[0],T_tmp,switch_promised,k_baseline,k_MELT,k_IDLE)

        # Expected true
        cons_momental_exptrue, cons_acc_exptrue = consumption(Z_rk_expected[0],T_tmp,switch_used,k_baseline,k_MELT,k_IDLE)

        # True - init switch - Promised
        cons_momental_true_init_promised, cons_acc_true_init_promised = consumption(Z_true_init_switch_promised[0],T_tmp,init_promised,k_baseline,k_MELT,k_IDLE)
        
        # True - init switch - rk
        cons_momental_true_init_rk, cons_acc_true_init_rk = consumption(Z_true_init_switch_rk[0],T_tmp,init_rk,k_baseline,k_MELT,k_IDLE)
        
        # True IDLE - zero price of running the tank
        cons_momental_true_idle, cons_acc_true_idle = consumption(Z_true_idle[0],T_tmp,switch_promised,k_baseline,0,0)

        

        # Compute payed prices
        price_payed =             1/1000000 * ((cons_acc_true - cons_acc_promised) * q_rk[hour] + cons_acc_promised * q_dap[hour])
        price_exptrue =           1/1000000 * ((cons_acc_exptrue - cons_acc_promised) * q_rk[hour] + cons_acc_promised * q_dap[hour])
        price_promised =          1/1000000 * ( cons_acc_promised * q_dap[hour])
        price_payed_init =        1/1000000 * ((cons_acc_true_init_rk - cons_acc_true_init_promised) * q_rk[hour] + cons_acc_true_init_promised * q_dap[hour])
        price_payed_idle =        1/1000000 * ( cons_acc_true_idle * q_dap[hour])



        # -----------------------------
        # Save monitored variables
        # -----------------------------
        
        start = int(1/n_skip*1/dt_stoch_sim*hour*60)
        end = int(1/n_skip*1/dt_stoch_sim*(hour+1)*60)
        history['Z_true'][day][0][start:end] = Z_true[:,::n_skip]
        history['Z_model_promised'][day][0][start:end] = Z_promised[:,::n_skip]
        history['Z_model_rk'][day][0][start:end] = Z_rk_expected[:,::n_skip]
        history['Z_true_init_promised'][day][0][start:end] = Z_true_init_switch_promised[:,::n_skip]
        history['Z_true_init_rk'][day][0][start:end] = Z_true_init_switch_rk[:,::n_skip]
        history['Z_true_idle'][day][0][start:end] = Z_true_idle[:,::n_skip]


        history['X_true'][day][:,start:end] = X_true[:,::n_skip]
        history['X_model_promised'][day][:,start:end] = X_promised[:,::n_skip]
        history['X_model_rk'][day][:,start:end] = X_rk_expected[:,::n_skip]
        history['X_true_init_promised'][day][:,start:end] = X_true_init_switch_promised[:,::n_skip]
        history['X_true_init_rk'][day][:,start:end] = X_true_init_switch_rk[:,::n_skip]
        history['X_true_idle'][day][:,start:end] = X_true_idle[:,::n_skip]


        history['T'][day][start:end] = T_tmp[::n_skip] # Same in each day
        history['SWITCHES_promised'][day,:] = switch_promised
        history['SWITCHES_rk'][day,:] = switch_used
        history['SWITCHES_init_promised'][day,:] = init_promised
        history['SWITCHES_init_rk'][day,:] = init_rk

        history['PRICES_dap'][day,:] = q_dap
        history['PRICES_up'][day,:] = q_up
        history['PRICES_down'][day,:] = q_down
        history['PRICES_rk'][day,:] = q_rk

        # Prices payed
        history['price_true'][day][start:end] = price_payed[::n_skip] + price_true_lasthour
        history['price_model_promised'][day][start:end] = price_promised[::n_skip] + price_promised_lasthour
        history['price_model_rk'][day][start:end] = price_exptrue[::n_skip] + price_rk_lasthour
        history['price_true_init_switch'][day][start:end] = price_payed_init[::n_skip] + price_true_lasthour_init_switch
        history['price_true_idle'][day][start:end] = price_payed_idle[::n_skip] + price_true_lasthour_idle

        history['price_true_each_hour'][day][start:end] = price_payed[::n_skip]
        history['price_model_promised_each_hour'][day][start:end] = price_promised[::n_skip]
        history['price_model_rk_each_hour'][day][start:end] = price_exptrue[::n_skip]
        history['price_true_each_hour_init_switch'][day][start:end] = price_payed_init[::n_skip]
        history['price_true_each_hour_idle'][day][start:end] = price_payed_idle[::n_skip]


        # Consumption
        history['cons_true'][day][start:end] = cons_acc_true[::n_skip] + cons_true_lasthour
        history['cons_model_promised'][day][start:end] = cons_acc_promised[::n_skip] + cons_promised_lasthour
        history['cons_model_rk'][day][start:end] = cons_acc_exptrue[::n_skip] + cons_rk_lasthour
        history['cons_true_init_promised'][day][start:end] = cons_acc_true_init_promised[::n_skip] + cons_true_lasthour_init_promised
        history['cons_true_init_rk'][day][start:end] = cons_acc_true_init_rk[::n_skip] + cons_true_lasthour_init_rk
        history['cons_true_idle'][day][start:end] = cons_acc_true_idle[::n_skip] + cons_true_lasthour_idle


        history['cons_true_each_hour'][day][start:end] = cons_acc_true[::n_skip]
        history['cons_model_promised_each_hour'][day][start:end] = cons_acc_promised[::n_skip]
        history['cons_model_rk_each_hour'][day][start:end] = cons_acc_exptrue[::n_skip]
        history['cons_true_each_hour_init_promised'][day][start:end] = cons_acc_true_init_promised[::n_skip]
        history['cons_true_each_hour_init_rk'][day][start:end] = cons_acc_true_init_promised[::n_skip]
        history['cons_true_each_hour_idle'][day][start:end] = cons_acc_true_idle[::n_skip]



        # Update variables
        start_true = X_true[:,-1]
        start_promised = X_promised[:,-1] 
        start_expected = X_rk_expected[:,-1] 
        start_true_init_promised = X_true_init_switch_promised[:,-1]
        start_true_init_rk = X_true_init_switch_rk[:,-1]
        start_true_idle = X_true_idle[:,-1]

        t0_sub += total_sim_time/24.
        tf_sub += total_sim_time/24.

        price_true_lasthour = history['price_true'][day][end-1]
        price_promised_lasthour = history['price_model_promised'][day][end-1]
        price_rk_lasthour = history['price_model_rk'][day][end-1]
        price_true_lasthour_init_switch = history['price_true_init_switch'][day][end-1]
        price_true_lasthour_idle = history['price_true_idle'][day][end-1]

        cons_true_lasthour = history['cons_true'][day][end-1]
        cons_promised_lasthour = history['cons_model_promised'][day][end-1]
        cons_rk_lasthour = history['cons_model_rk'][day][end-1]
        cons_true_lasthour_init_promised = history['cons_true_init_promised'][day][end-1]
        cons_true_lasthour_init_rk = history['cons_true_init_rk'][day][end-1]
        cons_true_lasthour_idle = history['cons_true_idle'][day][end-1]


        

    history['sim_time'][day] = time.time() - time_start

    history['end_date'] = prices.index[idx[0]]

    filenname = '../results/sim_history/rk_history_(' + start_date + ')_(' + str(n_days) + '_days)' + '_(price_slope_' + str(price_slope) +  ')_(regime_slope_' + str(regime_slope) +  ')_(seed_' + str(seed) +  ')_(n_s_' + str(n_s) + ')_(sys_model_' +  sys_mod + ')_(sys_true_' + sys_true + ')' +  '.npy'
    np.save(filenname,history)
