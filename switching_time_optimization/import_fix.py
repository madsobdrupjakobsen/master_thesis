#%%
import datetime as dt
import math
##
import os
import random
import sys
import time

import lxml
import matplotlib.patches as patch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize as sc_opt

# Switching time optimization modules
import switching_times_1st as st1
import switching_times_1st_rk as st1_rk
import switching_times_2nd as st2
import switching_times_2nd_rk as st2_rk
import switching_times_3rd as st3
import switching_times_3rd_rk as st3_rk

import torch
from IPython.core.debugger import set_trace
from numpy import genfromtxt
from scipy.integrate import quad, solve_ivp
from sklearn.preprocessing import MinMaxScaler
from torch import nn, optim
from torch.autograd import Variable
from torch.optim import lr_scheduler
from torch.utils.data import DataLoader
from torch.utils.data.dataset import Dataset

from ann.dense_nets import simpleNet
from data import SwitchData, SwitchData_rk
from model_loader import *
from models import firstordermodel, secondordermodel, thirdordermodel
from price_loader import *
from tools import (build_initial_ipopt_object, consMatrix, consumption, cost,
                   derive_regimes, discrete_ivp_solver, f_with_consumption,
                   f_with_cost, makeANNFeasible, mergeSwitchCont, plotSwitches,
                   predict_switch_ann_og, removeRedundantSwitches,
                   removeWrongOrder, sigmoid, simulate_MPC, smooth_dap,
                   smooth_regime, sol_ivp_wrapper, sol_ivp_wrapper_discrete,
                   stochasticSimulation, stochasticSimulation_IDLE)

sys.path.append('..')
%load_ext autoreload
%autoreload 2


sys.path.append('../../models')
sys.path.append('../')

os.system("../model_loader.py")
os.system("../price_loader.py")
os.system("../../models/models.py")
