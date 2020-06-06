import numpy as np
import scipy
from scipy.integrate import solve_ivp
import torch
import matplotlib.pyplot as plt

# Expand f to also compute the consumptions
def f_with_consumption(t,x,switches,model,price,k_baseline,k_IDLE,k_MELT):
    
    # first nx-2 values is the state of the process
    # -2 value is the total power consumption
    # -1 is the total cost
    
    model_regime = smooth_regime(t,switches,1.)

    nx = x.size
    x_no_cost = x[:-1]
    dxdt_no_cost = model.f(t,x_no_cost,switches)
    dxdt = np.zeros(len(dxdt_no_cost)+1)

    # Normal process
    for i in range(nx-1):
        dxdt[i] = dxdt_no_cost[i]

    # Consumption
    dxdt[nx-1] = (k_baseline * model.h(x_no_cost) +  model_regime * k_MELT + (1-model_regime) * k_IDLE) * 1./60;   
    return(dxdt)

# Expand f to also compute the cost
def f_with_cost(t,x,switches,model,price,k_baseline,k_IDLE,k_MELT,hour):
    
    # first nx-2 values is the state of the process
    # -2 value is the total power consumption
    # -1 is the total cost
    _price = smooth_dap(t,price)
    model_regime = smooth_regime(t,switches,1.)

    nx = x.size
    x_no_cost = x[:-2]
    dxdt_no_cost = model.f(t,x_no_cost,switches)
    dxdt = np.zeros(len(dxdt_no_cost)+2)

    # Normal process
    for i in range(nx-2):
        dxdt[i] = dxdt_no_cost[i]

    # Consumption
    dxdt[nx-2] = (k_baseline * model.h(x_no_cost) +  model_regime * k_MELT + (1-model_regime) * k_IDLE) * 1./60.;
    # Cost
    dxdt[nx-1] = _price * (k_baseline * model.h(x_no_cost) +  model_regime * k_MELT + (1-model_regime) * k_IDLE) * 1./60.;   
    return(dxdt)

def stochasticSimulation_IDLE(model,x0,t0,tf,dt):
    
    # Compute regimes
    #tau_MELT_all, tau_IDLE_all = derive_regimes(switches,tf,1)
    
    # Model functions
    f_IDLE = model.f_IDLE
    h = model.h
    sigma = model.sigma
    R = model.R
    
    # Initialize vairables
    x = x0
    #print(tf)
    #print(dt)
    N = int((tf-t0)/dt)
    nx = x.size
    nz = h(x).size
    X = np.zeros((nx,N))
    Y = np.zeros((nz,N))
    Z = np.zeros((nz,N))
    T = np.zeros((N))


    # Draw noise
    dW = np.sqrt(dt) * np.random.multivariate_normal(np.zeros(nx), np.eye(nx), N)
    dW = np.sqrt(dt) * np.random.multivariate_normal(np.zeros(nx), np.eye(nx), N)
    v = np.sqrt(R) * np.random.multivariate_normal(np.zeros(nz), np.eye(nz), N)


    # Run simulation
    t = t0
    for k in range(N):
        dx = f_IDLE(t,x)
        dw = sigma * dW[k]

        x = x + dx * dt + dw
        z = h(x)
        t += dt

        # Save current state
        X[:,k] = x
        Z[:,k] = z
        Y[:,k] = z + v[k]
        T[k] = t
        
    return np.squeeze(T), X, Y, Z

def stochasticSimulation(model,switches,x0,t0,tf,dt):
    
    # Compute regimes
    #tau_MELT_all, tau_IDLE_all = derive_regimes(switches,tf,1)
    
    # Model functions
    f = model.f
    h = model.h
    sigma = model.sigma
    R = model.R
    
    # Initialize vairables
    x = x0
    N = int((tf-t0)/dt)
    nx = x.size
    nz = h(x).size
    X = np.zeros((nx,N))
    Y = np.zeros((nz,N))
    Z = np.zeros((nz,N))
    T = np.zeros((N))


    # Draw noise
    dW = np.sqrt(dt) * np.random.multivariate_normal(np.zeros(nx), np.eye(nx), N)
    v = np.sqrt(R) * np.random.multivariate_normal(np.zeros(nz), np.eye(nz), N)


    # Run simulation
    t = t0
    for k in range(N):
        dx = f(t,x,switches)
        dw = sigma * dW[k]

        x = x + dx * dt + dw
        z = h(x)
        t += dt

        # Save current state
        X[:,k] = x
        Z[:,k] = z
        Y[:,k] = z + v[k]
        T[k] = round(t,10)
        
    return np.squeeze(T), X, Y, Z


def simulate_MPC(system_true,system_model,tank,
                 x0_true,x0_model,n_days,prices,k,k_MELT,k_IDLE,dt,
                 start_date,seed,save_to_file):   
    
    np.random.seed(seed)
    
    n_s = int(len(tank.get_p_optimize())/2)
    print(n_s)

    start_true = x0_true
    start_model = x0_model
    # Simulate parameters
    #n_days = 30
    #dt = 0.01
    #start_date = '01-08-2013'

    # Initialize process
    #x0 = np.array([2.6659426, 899.8884004])
    nx_true = start_true.shape[0]
    nx_model = start_true.shape[0]
    
    tf_ph = 48 * 60
    tf_sim = 24 * 60
    dt = 0.01

    # Initialzie history
    history = {}
    history['X_true'] = np.zeros([n_days, nx_true, 60 * 24])
    history['X_model'] = np.zeros([n_days, nx_model, 60 * 24])
    history['Z_true'] = np.zeros([n_days, 1, 60 * 24])
    history['Z_model'] = np.zeros([n_days, 1, 60 * 24])
    history['T'] = np.zeros([n_days, 60 * 24])

    history['SWITCHES'] = np.zeros([n_days, 2 * n_s])
    history['PRICES'] = np.zeros([n_days, 48])

    history['price_model'] = np.zeros(n_days)
    history['price_true'] = np.zeros(n_days)


    idx_offset = np.where(prices.index == start_date)[0][0]
    for day in range(n_days):
        #print(' --- Simulating day ' + str(day), end = '')
        print('  ----- Simulating day ' + str(day) + ' ----  ')

        # Extract variables related to the day
        idx = np.arange(48) + day * 48 + idx_offset
        future_days = np.array(prices.index[idx])
        #future_hours = np.array(prices['Hours'][idx])
        future_price = np.array(prices['spot'][idx]).astype(np.float) * 1/1000000


        # ---------------------------
        # Compute optimal switches over 2 days
        # ---------------------------
        
        # -- Set new optimization object
        
        # Initial state of the system
        tank.set_x0(np.append(start_model,0)) 
        
        # Update prices
        p_dynamic = tank.get_p_dynamic()
        for k in range(0, 48):
            p_dynamic[k] = future_price[k]
        tank.set_p_dynamic(p_dynamic)

        
        # -- Solve the problem - Uncommented until number of iterations is fixed
        print(' ... Optimizing')
        tank.solve()
        switch_opt = np.array(tank.get_p_optimize_ipopt())
        

        # To be changed to C++ optimizer
        #switch_opt = np.sort(np.random.uniform(0,tf_ph,2*n_s))
        #switch_opt = np.concatenate(derive_regimes(switch_opt,0,0))


        # ---------------------------
        # Simulate 1 day
        # ---------------------------
        switch_sim = switch_opt # Could be corrupted by delay is in reality

        # Consider only saving every 1 minutes of the simulation
        print(' ... Simulating true process')
        T_tmp, X_tmp, Y_tmp, Z_tmp = stochasticSimulation(system_true,switch_sim,start_true,tf_sim,dt)
        #T_tmp = np.linspace(0,tf_sim,1000)
        
        print(' ... Simulating model process')
        T_tmp_ode, X_tmp_ode, Z_tmp_ode = sol_ivp_wrapper_discrete(system_model,start_model,switch_sim,t0_sim,tf_sim,T_tmp)
        #T_tmp_ode, X_tmp_ode, Z_tmp_ode =  sol_ivp_wrapper(system_model,start_model,switch_sim,tf_sim,T_tmp)


        # ---------------------------
        # Compute Cost
        # ---------------------------
        dap = future_price[:24]
        cost_true, cost_acc_true = cost(Z_tmp[0],T_tmp,switch_sim,dap,k,k_MELT,k_IDLE)
        cost_model, cost_acc_model = cost(Z_tmp_ode[0],T_tmp_ode,switch_sim,dap,k,k_MELT,k_IDLE)



        # ---------------------------
        # Update variables
        # ---------------------------

        # Process
        #if nx > 1:
        #    start = X_tmp[i][-1] for i in range(nx)]
        #else:
        #    start = X_tmp[-1]
        #x0 = np.array([X_tmp[i][-1] for i in range(nx)])
        start_true = X_tmp[:,-1]
        start_model = X_tmp_ode[:,-1] # Should be estimated with Kalman filter.

        
        # Update monitored variables
        history['price_true'][day] = cost_true[-1]
        history['price_model'][day] = cost_model[-1]

        history['Z_true'][day] = Z_tmp[:,::int(1/dt)]
        history['Z_model'][day] = Z_tmp_ode[:,::int(1/dt)]
        history['X_true'][day,:] = X_tmp[:,::int(1/dt)]
        history['X_model'][day,:] = X_tmp_ode[:,::int(1/dt)]
        history['T'][day] = T_tmp[::int(1/dt)] # Same in each day
        history['SWITCHES'][day,:] = switch_opt
        history['PRICES'][day,:] = future_price

    filenname = '../results/sim_history/history_' + start_date + '_' + str(n_days) + '_days' + '_seed_' + str(seed) + '.npy'
    np.save(filenname,history)
    
    return history




def discrete_ivp_solver(f,x0,switch,T,tf,extra_args=()): 
    # Notice that the very last value in T is not evaluated
    
    # Extract eventual extra arguments to f
    args = (switch,) + extra_args
    
    # There must not be duplicates in switches
    nx = x0.size
    T_out = np.array([T[0]])
    X_out = np.expand_dims(x0, axis=1) #np.zeros((nx,1))
    
    # Derive the melt and idle starts
    n_s = int(len(switch)/2)
    tau_MELT_all = np.append(switch[:int(n_s)],tf)
    tau_IDLE_all = np.insert(switch[int(n_s):],0,0)
  


    # Find initial regime
    #t0_subint = T[0]
    #if smooth_regime(t0_subint,switch,2000000.) > 0.5:
    #    current_state = 1
    #else:
    #    current_state = 0
    
    t0_subint = T[0]
    current_state = np.sum(t0_subint >= tau_IDLE_all) == np.sum(t0_subint >= tau_MELT_all)

    


    x0_subint = x0
    tf_subint = T[-2] # Dummy value
    
    while (tf_subint < T[-1]):

        if current_state == 1:
            tf_subint = tau_IDLE_all[t0_subint < tau_MELT_all][0]
        else:
            tf_subint = tau_MELT_all[t0_subint >= tau_IDLE_all][-1]
            
        #print(t0_subint, tf_subint)


        # Derive t_eval
        #print('Tsize', T.size)
        T_subint = T[(t0_subint < T) & (T <= tf_subint)]
        #print(T_subint.size)
        if len(T_subint) > 0:
            #print('Tsubsize',T_subint.size)
            # Solve in the interval anf save
            #print(args)
            sol = solve_ivp(f, [t0_subint, T_subint[-1]], x0_subint, args=args, t_eval = T_subint)
            #print(sol.t.size)
            #print(X_out)
            #print(X_out.shape)
            #print(sol.y)
            #print(sol.y.shape)
            T_out = np.hstack((T_out,sol.t))
            X_out = np.hstack([X_out,sol.y])

            # Update x only if an iteragration was executed
            x0_subint = np.array([sol.y[j][-1] for j in range(nx) ])
            
        else:
            # Make sure to keep integration though we do not want to return any points in the region
            sol = solve_ivp(f, [t0_subint, tf_subint], x0_subint, args=args, t_eval = np.expand_dims(tf_subint, axis=0))

            # Update x only if an iteragration was executed
            x0_subint = np.array([sol.y[j][-1] for j in range(nx) ])

        # Update status of the system
        current_state = 1 - current_state
        t0_subint = tf_subint
        
    # Evaluate end point
    #sol = solve_ivp(f, [tf-1, tf], x0_subint, args=args, t_eval = np.array(tf))
    #T_out = np.concatenate((T_out, sol.t))
    #X_out = np.hstack([X_out,sol.y])

    
    # Restructure solution
    #X_out = np.array([X_out[i][1:] for i in range(nx)])
    #Z_out = model.h(X_out)
    
    return T, X_out #, Z_out


def sol_ivp_wrapper_discrete(model,x0,switches,t0,tf,T):
    if switches.size > 0:
        T, X = discrete_ivp_solver(model.f,x0,switches,T,tf)
    else:
        T[-1] = round(T[-1],10)
        #print(T[0],T[-1],t0,tf)
        sol = solve_ivp(model.f_IDLE, [t0, tf], x0,t_eval = T)
        T = sol.t
        X = sol.y
    
    Z = model.h(X)
    return T, X, Z
    

def sol_ivp_wrapper(model,x0,switches,tf,t_plot):
    x0 = [x0[i] for i in range(len(x0))]
    
    sol = solve_ivp(model.f, [0, tf], x0, args=(switches,),t_eval = t_plot)
    
    X = sol.y
    Z = model.h(X)
    T = sol.t
    return np.squeeze(T), X, Z

def derive_regimes(switches,tf,endpoints):
    n_switches = switches.size
    n_cycles = int(n_switches/2)
    tau_MELT = switches[np.arange(0,n_switches,2)]
    tau_IDLE = switches[np.arange(1,n_switches,2)]
    
    if endpoints:
        tau_IDLE = np.insert(tau_IDLE,0,0)
        tau_MELT = np.append(tau_MELT,tf)
    
    return tau_MELT, tau_IDLE 

def mergeSwitchCont(switches):
    n_s = int(len(switches)/2)
    melt = switches[:n_s]
    idle = switches[n_s:]
    return(np.insert(idle, np.arange(len(melt)), melt))

# ENERGY COST

def smooth_dap(t,dap,slope = 1,n_hours=24):
    dat = 60 * np.arange(dap.size + 1); dat[0] = dat[0] - 60; dat[-1] = dat[-1] + 60
    
    _dap = 0
    for k in range(n_hours):
            _dap += dap[k] / ((1. + np.exp(np.minimum(-slope * (t - dat[k]), 15.0 ))) *
                             (1. + np.exp( np.minimum(slope * (t - dat[k + 1]), 15.))))
    return _dap


def smooth_regime(T,switches,slope=20.):
    n_s = int(len(switches)/2)
    tau_MELT, tau_IDLE  = switches[:n_s] , switches[n_s:] #derive_regimes(switches,T[-1],0)
    regime = 0
    for k in range(len(tau_MELT)):
            regime += 1/ ((1 + np.exp(np.minimum(-slope* (T - tau_MELT[k]), 15.0 ))) *
                             (1 + np.exp( np.minimum(slope* (T - tau_IDLE[k]), 15.))))
            
    return regime

def consumption(traj,T,switches,k,k_MELT,k_IDLE):
    
    
    dt = np.diff(T)
    regime = smooth_regime(T,switches,20000000.)
    cost_momental =  (k * traj + k_IDLE * (regime <= 0.5) + k_MELT * (regime > 0.5)) * 1/60 # 1/60 to Compensate for unit change
    cost_acc = np.cumsum(cost_momental[1:] * dt)

    #plt.plot(cost_momental)
    return cost_momental, cost_acc


def cost(traj,T,switches,dap,k,k_MELT,k_IDLE):
    
    
    dt = np.diff(T)
    regime = smooth_regime(T,switches,1.)
    _dap = smooth_dap(T,dap)
    cost_momental = _dap * (k * traj + k_IDLE * (regime <= 0.5) + k_MELT * (regime > 0.5)) * 1/60 # 1/60 to Compensate for unit change
    cost_acc = np.cumsum(cost_momental[1:] * dt)

    #plt.plot(cost_momental)
    return cost_momental, cost_acc


# MATH
def sigmoid(x,slope,offset):
    return 1./(1. + np.exp(-(slope * (x - offset))))

def plotSwitches(switches,t0,c1,c2,ax):
    n_s = int(len(switches)/2)
    for xc in switches[:n_s]:
        ax.axvline(x=xc+t0,color = c1, alpha = 0.2)
    
    for xc in switches[n_s:]:
        ax.axvline(x=xc+t0,color = c2, alpha = 0.2)

def removeRedundantSwitches(switches,tolerance):
    current_n_s = int(len(switches)/2)
    curr_melt = switches[:current_n_s]
    curr_idle = switches[current_n_s:]

    rm = []
    for i in range(current_n_s):
        if abs(curr_melt[i] - curr_idle[i]) < tolerance:
            rm = np.append(rm,[i,current_n_s + i])

    switches = np.delete(switches, np.array(rm).astype(int))
    rm = []
    current_n_s = int(len(switches)/2)
    curr_melt = switches[:current_n_s]
    curr_idle = switches[current_n_s:]

    for i in range(current_n_s-1):
        if abs(curr_idle[i] - curr_melt[i+1]) < tolerance:
            rm = np.append(rm,[current_n_s + i, i+1])

    switches = np.delete(switches, np.array(rm).astype(int))
    
    return(switches)

def removeWrongOrder(switches):
    current_n_s = int(len(switches)/2)
    curr_melt = switches[:current_n_s]
    curr_idle = switches[current_n_s:]

    rm = []
    for i in range(current_n_s):
        if (curr_melt[i] - curr_idle[i] > 0):
            rm = np.append(rm,[i,current_n_s + i])

    switches = np.delete(switches, np.array(rm).astype(int))
    rm = []
    current_n_s = int(len(switches)/2)
    curr_melt = switches[:current_n_s]
    curr_idle = switches[current_n_s:]

    for i in range(current_n_s-1):
        if (curr_idle[i] - curr_melt[i+1] > 0):
            rm = np.append(rm,[current_n_s + i, i+1])

    switches = np.delete(switches, np.array(rm).astype(int))
    
    return(switches)

def makeANNFeasible(pred_switch,redundant_tolerance = 10):
        correct_order = removeWrongOrder(pred_switch)
        filtered_switch = removeRedundantSwitches(correct_order,10)
        new_ns = int(len(filtered_switch)/2)

        melt = filtered_switch[new_ns:] - filtered_switch[:new_ns]
        
        exceed = np.sum(melt) - 16*60
        adjust = melt / np.sum(melt) * exceed / 2

        if exceed > 0:
            filtered_switch[new_ns:] -= adjust
            filtered_switch[:new_ns] += adjust
        
        return filtered_switch

def predict_switch_ann_rk(model,input,switch_type):
    
    model.load_state_dict(torch.load('./trained_networks/best_model_rk_switchtype_' + switch_type))
    model.eval()

    # Convert predictions to actual switches
    output = model(torch.tensor(input).float())
    switches_conological = np.array(output.data[0])

    switches = np.concatenate(derive_regimes(switches_conological,2*1440,0))
    
    return switches

def predict_switch_ann_og(model,input):
    
    #model.load_state_dict(torch.load('./trained_networks/best_model_rk_switchtype_' + switch_type))
    #model.eval()

    # Convert predictions to actual switches
    output = model(torch.tensor(input).float())
    switches_conological = np.array(output.data[0])

    switches = np.concatenate(derive_regimes(switches_conological,2*1440,0))
    
    return switches



def build_initial_ipopt_object(tank, tank_pars, dt, k, k_MELT, k_IDLE, t0, tf_ph, max_melt, n_s, price_slope, regime_slope,opt_rk):
    
    #p_x0 = np.append(x0,0)
    p_const = np.concatenate((tank_pars, np.array([k, k_MELT, k_IDLE, price_slope, regime_slope])))
    #p_const = np.concatenate((tank_pars, np.array([k, k_MELT, k_IDLE])))
    
    
    if opt_rk:
        p_dynamic = np.zeros(48 + 48 + 49) # day-ahead prices and day-ahead times -> vary in each new optimization -> therefore dynamic
        for k in range(0, 48):
            p_dynamic[k] = -1 #future_price[k] * 1/1000000 Price noit set untill MPC loop
            p_dynamic[k + 48] = -1 #future_price[k] * 1/1000000 Price noit set untill MPC loop
        for k in range(0, 49):
            p_dynamic[48 + 48 + k] = k * 60.
            
        n_bounds = 4 * n_s
            
    else:
        p_dynamic = np.zeros(48 + 49) # day-ahead prices and day-ahead times -> vary in each new optimization -> therefore dynamic
        for k in range(0, 48):
            p_dynamic[k] = -1 #future_price[k] * 1/1000000 Price noit set untill MPC loop
        for k in range(0, 49):
            p_dynamic[48 + k] = k * 60.
            
            
        n_bounds = 2 * n_s
        

    lower_bound = np.zeros(n_bounds) # lower bounds for optimization variables
    upper_bound = np.zeros(n_bounds) # upper bounds for optimization variables
    on_bound = np.array([0., tf_ph, max_melt])
    off_bound = np.array([0., tf_ph])
    for k in range(0, n_bounds):
        lower_bound[k], upper_bound[k] = t0, tf_ph


    # Build initial values - Other methods could be considered
    #idle = tf_ph * np.linspace(0.1,0.9,n_s) # Spread out on the whole interval
    #melt = idle - max_melt/n_s * 0.9 # Assign melt period to a little before idle
    #p_optimize = np.concatenate((melt,idle)) # put together

    #print(upper_bound)
    tank.set_p_const(p_const)
    tank.set_p_dynamic(p_dynamic)
    #tank.set_p_optimize(p_optimize)
    tank.set_t0(t0)
    tank.set_tf(tf_ph)
    tank.set_dt(dt)
    #tank.set_x0(p_x0)
    tank.set_lower_bound(lower_bound)
    tank.set_upper_bound(upper_bound)
    tank.set_on_bound(on_bound)
    tank.set_off_bound(off_bound)
    
    return tank

def consMatrix(n_s,tf,max_melt):
    # Consider to rebuild since t0 > 0 and tN < tf is implciit since all variables must be greater than 0 and less than tf
    #nswitch = len(tau)
    n_s_all = 2*n_s
    A = np.zeros((n_s_all+2,n_s_all))
    A[:n_s_all,:n_s_all] += + -1 * np.eye(n_s_all)
    A[:n_s,n_s:] += + np.eye(n_s)
    A[(n_s+1):(n_s_all+1),:n_s] += np.eye(n_s)
    A[-1,:] = np.concatenate((np.ones(n_s),-1 * np.ones(n_s)))
    
    A = np.hstack((A[:,n_s:],A[:,:n_s])) # Let [tau_idle , tau_melt]' become [tau_melt, tau_idle]' 
    
    ub = np.zeros(n_s_all + 2)
    ub[-2] = tf
    ub[-1] = max_melt
    
    return(A,ub)


def constraintASparse(n_s):

    row = np.zeros( 2 * 6 * n_s - 2 )
    col = np.zeros( 2 * 6 * n_s - 2 )
    val = np.zeros( 2 * 6 * n_s - 2 )

    # Build coordinates

    # Block 1
    index = 0
    for k in range(n_s):
        row[index] = k; col[index] = k;
        index += 1

        row[index] = k; col[index] = n_s + k;
        index += 1

    for k in range(n_s-1):
        row[index] = n_s + k; col[index] = 1 + k;
        index += 1

        row[index] = n_s + k; col[index] = n_s + k;
        index += 1

    for k in range(n_s):
        row[index] = 2*n_s-1; col[index] = k;
        index += 1

        row[index] = 2*n_s-1; col[index] = n_s + k;
        index += 1

    # Block 2
    for k in range(n_s):
        row[index] = 2*n_s + k; col[index] = 2*n_s + k;
        index += 1

        row[index] = 2*n_s + k; col[index] = 2*n_s + n_s + k;
        index += 1

    for k in range(n_s-1):
        row[index] = 2*n_s + n_s + k; col[index] = 2*n_s + 1 + k;
        index += 1

        row[index] = 2*n_s + n_s + k; col[index] = 2*n_s + n_s + k;
        index += 1

    for k in range(n_s):
        row[index] = 2*n_s + 2*n_s-1; col[index] = 2*n_s + k;
        index += 1

        row[index] = 2*n_s + 2*n_s-1; col[index] = 2*n_s + n_s + k;
        index += 1


    # Insert values

    # Block 1
    index = 0
    for k in range(n_s):
        val[index] = 1
        index += 1

        val[index] = -1
        index += 1

    for k in range(n_s-1):
        val[index] = -1
        index += 1

        val[index] = 1
        index += 1

    for k in range(n_s):
        val[index] = -1
        index += 1

        val[index] = 1
        index += 1

    # Block 2
    for k in range(n_s):
        val[index] = 1
        index += 1

        val[index] = -1
        index += 1

    for k in range(n_s-1):
        val[index] = -1
        index += 1

        val[index] = 1
        index += 1

    for k in range(n_s):
        val[index] = -1
        index += 1

        val[index] = 1
        index += 1

    A = scipy.sparse.coo_matrix((val, (row, col)), shape=(2 * 2*n_s, 2 * 2*n_s))
    
    return A

def setupAxis(ax,scale=1,ncol_legend=1,nolegend=False):
    
    #ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #ax.rc('font', **{'size':str(int(30/scale))})

    # Show the major grid lines with dark grey lines
    if not nolegend:
        ax.legend(prop={'size': 30/scale},ncol=ncol_legend)
    ax.grid()
    ax.grid(b=True, which='major', color='#666666', linestyle='-',alpha=0.4)

    # Show the minor grid lines with very faint and almost transparent grey lines
    ax.minorticks_on()
    ax.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    #ax.set_ylim(0,1000)
    ax.set_xlabel(r'$123$', fontsize=40/scale)

    ax.set_ylabel(r'$123$', fontsize=40/scale)

    ax.tick_params(axis="x", labelsize=30)
    ax.tick_params(axis="y", labelsize=30)
    
    