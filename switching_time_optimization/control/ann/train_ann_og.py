#%%

##
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from matplotlib import rc
from sklearn.preprocessing import MinMaxScaler
from torch import nn, optim
from torch.autograd import Variable
from torch.optim import lr_scheduler
from torch.utils.data import DataLoader
from torch.utils.data.dataset import Dataset


from dense_nets import simpleNet

sys.path.append('..')
from data import SwitchData
from tools import (cost, derive_regimes, removeRedundantSwitches,
                   removeWrongOrder, sigmoid, smooth_dap, smooth_regime,
                   sol_ivp_wrapper, stochasticSimulation)


#%%
def my_loss_sum(output, target):
    x = torch.linspace(0, 2*24*60, 1000)

    loss = 0
    for k in range(output.shape[0]):
        pred_switch_torch =output[k]
        gt_switch_torch =target[k]
        tau_MELT = pred_switch_torch[0::2] 
        tau_IDLE = pred_switch_torch[1::2]
        tau_MELT_gt = gt_switch_torch[0::2] 
        tau_IDLE_gt = gt_switch_torch[1::2]

        #tau_MELT, tau_IDLE  = pred_switch_torch[:n_s] , pred_switch_torch[n_s:] #derive_regimes(switches,T[-1],0)
        regime = 0
        regime_gt = 0
        for k in range(6):
                regime += 1/ ((1 + torch.exp(torch.min(-0.2* (x - tau_MELT[k]), torch.tensor(15.0) ))) *
                                 (1 + torch.exp( torch.min(0.2* (x - tau_IDLE[k]),  torch.tensor(15.0)))))

                regime_gt += 1/ ((1 + torch.exp(torch.min(-0.2* (x - tau_MELT_gt[k]), torch.tensor(15.0) ))) *
                                 (1 + torch.exp( torch.min(0.2* (x - tau_IDLE_gt[k]),  torch.tensor(15.0)))))

    #loss = torch.sum(torch.pow(output - target,2))
        loss += torch.sum(torch.pow(regime - regime_gt,2))
    
    return loss

lossFunction = sys.argv[1]

#%%
model_sys = 'm1'
n_s = 6
filename = '../results/sim_history/og_history_(2018-01-01 12:00:00)_(100_days)' + '_(price_slope_' + str(0.2) +  ')_(regime_slope_' + str(0.2) +  ')_(seed_1235)_(n_s_' + str(n_s) + ')_(sys_model_' + model_sys + ')_(sys_true_m3).npy'
history = np.load(filename,allow_pickle=True).item()

n_s_all = history['SWITCHES_dap'].shape[1]


scaler = MinMaxScaler()


# Create training data set
batch_size = 4
train_range = range(0,80)
dset = SwitchData(history,train_range[0],train_range[-1]+1,transform=scaler.fit_transform )
train_loader = DataLoader(dset,
                          batch_size=batch_size,
                          shuffle=True,
                          num_workers=0,
                          pin_memory=True,# CUDA only
                          )


# Create test data set
batch_size_pred = 20
test_range = range(80,100)
dtest = SwitchData(history,test_range[0],test_range[-1]+1,transform=scaler.transform)


test_loader = DataLoader(dtest,
                              batch_size=batch_size_pred,
                              shuffle=False,
                              num_workers=0,
                              pin_memory=True)  # CUDA only)



#%% Network settings
learning_rate = 0.001
model = simpleNet(n_s_all=n_s_all)
optimizer = optim.RMSprop(model.parameters(), lr=learning_rate, weight_decay=0.0)  # n
#assert lossFunction != '-f'

print(lossFunction)
if lossFunction == 'mse':
    criterion = nn.MSELoss(reduction='sum')
elif lossFunction == 'l2':
    criterion = my_loss_sum


#criterion = my_loss_sum
# Train network
best_val_loss = None
max_epochs = 10000
losses_sgd = []
losses_train = []
losses_test = []
it = 0
for epoch in range(max_epochs):
    
    # Train
    model.train()
    loss_ = 0.
    predicted = []
    gt = []
    for batch_idx, (data, target) in enumerate(train_loader):
        #data = Variable(data.permute(0, 2, 1)).contiguous()
        target = Variable(target.unsqueeze_(0))
        optimizer.zero_grad()
        if target.data.size()[1] == batch_size:
            output = model(data)
            loss = criterion(output, target[0])
            losses_sgd.append(loss.data.item()/batch_size)
            loss_ += loss.data
            loss.backward()
            optimizer.step()
            
            for k in range(batch_size):
                predicted.append(output.data[k, 0])
                gt.append(target.data[:, k, 0])
        it += 1
    losses_train.append(loss_.item()/(len(train_loader) * batch_size))
    
    
    # Test
    model.eval()
    loss_ = 0.
    with torch.no_grad():
        for batch_idx, (data_test, target_test) in enumerate(test_loader):
        #data = Variable(data.permute(0, 2, 1)).contiguous()
            target_test = Variable(target_test.unsqueeze_(1))
            if target_test.data.size()[0] == batch_size_pred:
                output_test = model(data_test)
                loss = criterion(output_test, target_test[:,0,:])
                loss_ += loss.data
        

    losses_test.append(loss_.item()/(len(test_loader) * batch_size_pred))
    
    if not best_val_loss or loss_.item() < best_val_loss:
        torch.save(model.state_dict(), './trained_networks/best_model_og_loss' + lossFunction)
        optimal_number_of_epochs = epoch
        best_val_loss = loss_.item()
    
    if epoch % 10 == 9:    # print every 100th mini-batches
            print(f'Epoch {epoch+1}, Training loss: {losses_train[-1]:.2e}, Testing loss: {losses_test[-1]:.2e}, Learning rate {float(optimizer.param_groups[0]["lr"]):.2e}')
    #writer.add_scalar("loss_epoch", loss_, i)

    #scheduler_model.step()

#%%
np.save('./train_info/og_loss_epoch_loss' + lossFunction + '.npy',np.vstack([losses_train, losses_test]))
np.save('./train_info/og_loss_iteration_loss' + lossFunction + '.npy',losses_sgd)
np.save('./train_info/og_opt_epoch_loss' + lossFunction + '.npy',optimal_number_of_epochs)


# %%
