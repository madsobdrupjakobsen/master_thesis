import numpy as np
from tools import discrete_ivp_solver, f_with_consumption



def obj_rk(switches,number_of_hours,dt,model,x0,q_dap,q_rk,full_return):
    
    k_baseline = 12400.
    k_MELT = 250.
    k_IDLE = 9.
    
    n_s = int(len(switches)/4)
    switch_dap = switches[:(2*n_s)]
    switch_rk = switches[(2*n_s):]
    #print('dap', switch_dap)
    #print('rk', switch_rk)


    assert (number_of_hours <= int(q_dap.size)), 'Too short price array'

    # Extract variables from input
    n_obs = int(number_of_hours * 60/dt)
    #print(n_obs)
    tf = 60 * number_of_hours

    # Prepare variables
    x0_consumption_dap = np.concatenate((x0,np.array([0.])))
    x0_consumption_rk = np.concatenate((x0,np.array([0.])))

    nx = x0_consumption_rk.size
    T_out = np.zeros((n_obs))
    total_price_regulated = np.zeros((n_obs))
    total_price_spot = np.zeros((n_obs))
    total_price_regulated_sum = np.zeros((n_obs))
    total_price_spot_sum = np.zeros((n_obs))

    X_out_dap = np.zeros((nx,n_obs))
    X_out_rk = np.zeros((nx,n_obs))

    current_acc_spot_price = 0.
    current_acc_rk_price = 0.


    for hour in range(number_of_hours):
        #print(current_acc_rk_price)
        #print(current_acc_spot_price)


        #print(hour)
        # Takes 0.2 second to execute the solvers for dt = 0.001
        # Run for 60 minutes in hour x. Add dt * 1/900 to avoid collision with switching times
        # which currently cannot be handled by discrete_ivp_solver()
        T = dt * 1/900 + hour * 60 + dt * np.arange(0,int(60/dt))
        #print(T)
        T_hour, X_dap = discrete_ivp_solver(f_with_consumption,x0_consumption_dap,switch_dap,T,tf,extra_args=(model,q_dap,k_baseline,k_IDLE,k_MELT))
        T_hour, X_rk = discrete_ivp_solver(f_with_consumption,x0_consumption_rk,switch_rk,T,tf,extra_args=(model,q_rk,k_baseline,k_IDLE,k_MELT))
        #sol_dap = solve_ivp(f_with_cost, [T[0], T[-1]], x0_cost, args=(model,q_dap,switch_dap,k,k_IDLE,k_MELT),t_eval=T)
        #sol_rk = solve_ivp(f_with_cost, [T[0], T[-1]], x0_cost, args=(model,q_rk,switch_rk,k,k_IDLE,k_MELT),t_eval=T)

        x0_consumption_dap = np.concatenate(([X_dap[i][-1] for i in range(X_dap.shape[0]-1)], np.array([0.])))
        x0_consumption_rk = np.concatenate(([X_rk[i][-1] for i in range(X_rk.shape[0]-1)], np.array([0.])))



        #print(X_dap[-1][-1])
        #print(X_rk[-1][-1])
        
        # Expected use
        #exp_use = X_dap[-2]

        # Difference from expected
        #diff_use = X_rk[-2] - X_dap[-2]

        total_rk_price = (X_rk[-1] - X_dap[-1]) * q_rk[hour] +  (X_dap[-1]) * q_dap[hour]
        total_spot_price = X_dap[-1] * q_dap[hour]

        # Save values
        idx = np.array(range(0,int(60/dt))) + hour * int(60/dt)
        #print('Index',idx)
        T_out[idx] = T_hour
        X_out_dap[:,idx] = X_dap
        X_out_rk[:,idx] = X_rk

        # Save all accumulated prices for this very hour
        total_price_regulated[idx] = total_rk_price
        total_price_spot[idx] = total_spot_price

        # Save all accumulated prices for after hour
        total_price_regulated_sum[idx] = total_rk_price + current_acc_rk_price
        total_price_spot_sum[idx] = total_spot_price + current_acc_spot_price

        # Update current accumulated price
        current_acc_rk_price = total_price_regulated_sum[idx[-1]]
        current_acc_spot_price = total_price_spot_sum[idx[-1]]  

        #print(100 * ((total_rk_price - total_spot_price)[-1]) / total_spot_price[-1])

    # Choose to return data structures or objective function
    if full_return:
        out_dict = {}
        out_dict['T'] = T_out
        out_dict['X_dap'] = X_out_dap
        out_dict['X_rk'] = X_out_rk
        out_dict['total_price_regulated'] = total_price_regulated
        out_dict['total_price_spot'] = total_price_spot
        out_dict['total_price_regulated_sum'] = total_price_regulated_sum
        out_dict['total_price_spot_sum'] = total_price_spot_sum
        
        print('Done')
        
        return out_dict
    #T, X_out_dap, X_out_rk, total_price_regulated, total_price_spot, total_price_regulated_sum, total_price_spot_sum
    else:
        return total_price_regulated_sum[-1]

