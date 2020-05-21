import os

import torch
from torch.utils.data.dataset import Dataset
from sklearn.preprocessing import MinMaxScaler
import pandas as pd
import numpy as np

class SwitchData(Dataset):
    def __init__(self, history, day_start, day_end,transform=None):
        """
        :param folder_dataset: str
        :param T: int
        :param symbols: list of str
        :param use_columns: bool
        :param start_date: str, date format YYY-MM-DD
        :param end_date: str, date format YYY-MM-DD
        """

        assert transform is not None
        
        self.history=history
        self.day_start=day_start
        self.day_end=day_end

        self.z0 = self.history['Z_true'][self.day_start:self.day_end,0,0]
        self.price = self.history['PRICES_dap'][self.day_start:self.day_end,:]
        self.switch = self.history['SWITCHES_dap'][self.day_start:self.day_end]
        self.n_s = int(self.switch.shape[1]/2)
        #print(self.n_s)

        self.numpy_data = np.column_stack((self.z0, self.price))
        #self.train_data = torch.FloatTensor(self.scaler.fit_transform(self.numpy_data))
        
        self.x_all = torch.FloatTensor(transform(self.numpy_data))

        # Un order to make i t easier to learn the cronology
        self.switch_new = np.zeros((self.x_all.shape[0],2*self.n_s))
        self.switch_new[:,0::2] = self.switch[:,:self.n_s]
        self.switch_new[:,1::2] = self.switch[:,self.n_s:]
        self.y_all = torch.FloatTensor(self.switch_new)


    def __getitem__(self, index):

        x = self.x_all[index,:]
        y = self.y_all[index,:]
        return x, y

    def __len__(self):
        return self.x_all.shape[0]

class SwitchData_rk(Dataset):
    def __init__(self, history, day_start, day_end,type,transform=None):
        """
        :param folder_dataset: str
        :param T: int
        :param symbols: list of str
        :param use_columns: bool
        :param start_date: str, date format YYY-MM-DD
        :param end_date: str, date format YYY-MM-DD
        """
        
        self.history=history
        self.day_start=day_start
        self.day_end=day_end
        self.type = type

        assert transform is not None


        self.z0 = self.history['Z_true'][self.day_start:self.day_end,0,0]
        self.price_rk = self.history['PRICES_rk'][self.day_start:self.day_end,:]
        self.price_dap = self.history['PRICES_dap'][self.day_start:self.day_end,:]

        if self.type  == 'rk':
            self.switch = self.history['SWITCHES_rk'][self.day_start:self.day_end]
        else:
            self.switch = self.history['SWITCHES_promised'][self.day_start:self.day_end]
        self.n_s = int(self.switch.shape[1]/2)
        #print(self.n_s)

        self.numpy_data = np.column_stack((self.z0, self.price_rk, self.price_dap))
        #self.train_data = torch.FloatTensor(self.scaler.fit_transform(self.numpy_data))
        
        self.x_all = torch.FloatTensor(transform(self.numpy_data))

        # Un order to make i t easier to learn the cronology
        self.switch_new = np.zeros((self.x_all.shape[0],2*self.n_s))
        self.switch_new[:,0::2] = self.switch[:,:self.n_s]
        self.switch_new[:,1::2] = self.switch[:,self.n_s:]
        self.y_all = torch.FloatTensor(self.switch_new)


    def __getitem__(self, index):

        x = self.x_all[index,:]
        y = self.y_all[index,:]
        return x, y

    def __len__(self):
        return self.x_all.shape[0]

    
class MarketClass(Dataset):
    def __init__(self, priceData):
        """
        :param Pandas dataframe with price information
        """
        
        #history_np=np.load(filename + '.npy',allow_pickle=True)
        self.prices = priceData

        self.scaler_train = MinMaxScaler()
        #self.scaler_test = MinMaxScaler()

        self.numpy_data = np.column_stack((self.prices['spot'], self.prices['t_sin']))
        #self.train_data = torch.FloatTensor(self.scaler.fit_transform(self.numpy_data))
        
        self.x_all = torch.FloatTensor(self.scaler_train.fit_transform(self.numpy_data))
        self.y_all = torch.LongTensor(self.prices['market_class_all'])


    def __getitem__(self, index):

        x = self.x_all[index,:]
        y = self.y_all[index]
        return x, y

    def __len__(self):
        return self.x_all.shape[0]