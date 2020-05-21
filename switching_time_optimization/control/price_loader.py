import pandas as pd
import numpy as np

regulation = pd.read_csv('../../data/realtimemarket.csv',header=0, index_col=0,skiprows = 0)
regulation = regulation.loc[regulation['PriceArea'] == 'DK1']
regulation = regulation[{'BalancingPowerPriceDownEUR','BalancingPowerPriceUpEUR'}]
regulation.index = pd.to_datetime(regulation.index)
regulation = regulation.reindex(index=regulation.index[::-1])
spot = pd.read_csv('../../data/elspotprices.csv',header=0, index_col=0,skiprows = 0)
spot = spot.loc[spot['PriceArea'] == 'DK1']
spot = spot[{'SpotPriceEUR'}]
spot.index = pd.to_datetime(spot.index)
spot = spot.reindex(index=spot.index[::-1])

prices = pd.concat([spot, regulation['BalancingPowerPriceDownEUR'], regulation['BalancingPowerPriceUpEUR']], axis=1, join='inner')
prices.columns = ['spot', 'down', 'up']
#prices = prices * 1e-6


# Make correct round off
num_deci = 2
prices_og = prices
prices = np.round(np.round(prices_og * 10**num_deci,num_deci) / 10**num_deci,num_deci) #= np.round(prices,3)

# Derive markets
tol = 1e-2

prices['market_spot'] = 0
prices['market_down'] = 0
prices['market_up'] = 0

# When both up and down is less than tol away from the spot price
prices.loc[(np.abs(np.round(prices['spot'] - prices['down'],2)) <= tol) & (np.abs(np.round(prices['spot'] - prices['up'],2)) <= tol),'market_spot'] = 1

# When both are different, but even the smallest is bigger than 0.1 we have an error. Put it to spot market
prices.loc[np.minimum(np.abs(prices['spot'] - prices['up']),np.abs(prices['spot'] - prices['down'])) > 0.1,'market_spot'] = 1


prices.loc[(prices['market_spot'] != 1) & (np.abs(np.round(prices['spot'] - prices['down'],2)) > np.abs(np.round(prices['spot'] - prices['up'],2))),'market_down'] = 1
prices.loc[(prices['market_spot'] != 1) & (np.abs(np.round(prices['spot'] - prices['down'],2)) <= np.abs(np.round(prices['spot'] - prices['up'],2))),'market_up'] = 1


frac_down = np.mean(prices['market_down'])
frac_spot = np.mean(prices['market_spot'])
frac_up = np.mean(prices['market_up'])

prices['RK'] = prices['spot'] * prices['market_spot'] + prices['up'] * prices['market_up'] + prices['down'] * prices['market_down']
#prices['RK'] = -50 

#[frac_down, frac_spot, frac_up], np.sum([frac_down, frac_spot, frac_up])

# Repeat missing hours and add hour variable
idx = pd.date_range(prices.index[0], prices.index[-1], freq='H')
prices.index = pd.DatetimeIndex(prices.index)
prices = prices.reindex(idx, fill_value='h')

idx_nan = np.where(prices['spot'] == 'h')
n_missing = len(idx_nan[0])
for i in range(n_missing):
    prices.iloc[idx_nan[0][i]] = prices.iloc[np.array(idx_nan[0][i])-1]
    
prices['t_sin'] = np.sin(prices.index.hour/(24/np.pi))
prices['t_sin'] = np.sin(2*np.pi*prices.index.hour/24)

prices['market_class_all'] = -1
prices['market_class_all'].loc[prices['market_down']==1] = 0
prices['market_class_all'].loc[prices['market_spot']==1] = 1
prices['market_class_all'].loc[prices['market_up']==1] = 2