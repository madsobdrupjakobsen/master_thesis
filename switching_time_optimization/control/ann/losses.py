import torch

def L2diff(output, target):
    x = torch.linspace(0, 2*24*60, 1000)

    loss = 0
    for k in range(output.shape[0]):
        pred_switch_torch =(2*1440)*output[k]
        gt_switch_torch =(2*1440)*target[k]
        tau_MELT = pred_switch_torch[0::2] 
        tau_IDLE = pred_switch_torch[1::2]
        tau_MELT_gt = gt_switch_torch[0::2] 
        tau_IDLE_gt = gt_switch_torch[1::2]

        #tau_MELT, tau_IDLE  = pred_switch_torch[:n_s] , pred_switch_torch[n_s:] #derive_regimes(switches,T[-1],0)
        regime = 0
        regime_gt = 0
        for k in range(6):
                regime += 1/ ((1 + torch.exp(torch.min(-1* (x - tau_MELT[k]), torch.tensor(15.0) ))) *
                                 (1 + torch.exp( torch.min(1* (x - tau_IDLE[k]),  torch.tensor(15.0)))))

                regime_gt += 1/ ((1 + torch.exp(torch.min(-1* (x - tau_MELT_gt[k]), torch.tensor(15.0) ))) *
                                 (1 + torch.exp( torch.min(1* (x - tau_IDLE_gt[k]),  torch.tensor(15.0)))))

    #loss = torch.sum(torch.pow(output - target,2))
        loss += torch.sum(torch.pow(regime - regime_gt,2))
    
    return loss


def MSECons(output,target,mu_pos = 0):
    pos_cons = torch.mean(torch.exp(torch.min(-1 * output,torch.tensor(15.))))
    #pos_cons = torch.mean(1/(1 + torch.exp(torch.min(3*output,torch.tensor(5.)))))
    #melt = output[:,1::2] - output[:,0::2]
    #idle = output[:,2:-1:2] - output[:,1:-1:2]
    #order_violation = torch.mean(torch.exp(-100 * idle)) + torch.mean(torch.exp(-100 * melt))
    #positive_violation = torch.mean(torch.exp(-10 * output))
    
    return 1/10000 * torch.mean(torch.pow(output - target,2)) + mu_pos*pos_cons #+ positive_violation