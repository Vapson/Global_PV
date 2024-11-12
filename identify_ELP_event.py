# -*- coding: utf-8 -*-
"""
Created on Sat Aug 13 17:30:16 2022

@author: wangx
"""
import os
os.chdir(r'G:\global_photovoltatic')


import numpy as np
import pandas as pd
import datetime
import xarray as xr


# In[event without backup energy]

def main_0():
 
    for year in range(1986, 2022):
        file_path = r'G:\global_photovoltatic\daily_power_by_ERA5_pvlib'
        file_name = 'power_daily_' + str(year) + '.nc'
        daily_data = xr.open_dataset( os.path.join(file_path, file_name) )['band_data']
        days = len(daily_data['time'])
        lat_length, lon_length = len(daily_data.lat), len(daily_data.lon)
        
        if days == 366:
            daily_data = daily_data[~((daily_data['time'].dt.month==2) & (daily_data['time'].dt.day==29))]
            days = 365
        
        # generate time, lon_index, lat_index
        doy = np.arange(0,days,1)
        lat = np.arange(0,lat_length,1)
        lon = np.arange(0,lon_length,1)
        lats, doys, lons  = np.meshgrid(lat,doy, lon)         
        
        
        # read threshold
        thresholds = xr.open_dataset(r'G:\global_photovoltatic\threshold\threshold2.nc')
        
        
        threshold = thresholds.sel(thrs='10th')['__xarray_dataarray_variable__'].values
        threshold1 = thresholds.sel(thrs='5th')['__xarray_dataarray_variable__'].values
        threshold2 = thresholds.sel(thrs='1th')['__xarray_dataarray_variable__'].values
        
    
        power = daily_data.values
        judge = np.where((power<threshold)|(power<0), 1, 0)
        lon_index = lons[judge==1]
        lat_index = lats[judge==1]
        doy_index = doys[judge==1]
        power_event = power[judge==1]
        thr = threshold[judge==1]
        
        thr1 = threshold1[judge==1]
        thr2 = threshold2[judge==1]
        
        # columns=['lon_index','lat_index','threshold','country_FID','doy','power']
        results = np.concatenate((lon_index.reshape(-1,1), lat_index.reshape(-1,1), thr.reshape(-1,1), \
                                  thr1.reshape(-1,1), thr2.reshape(-1,1),\
                                  doy_index.reshape(-1,1), power_event.reshape(-1,1)), axis=1, dtype=np.float16) 
        save_path = os.path.join(r'G:\global_photovoltatic\low_output') 
        save_name = 'event_' + str(year) + '.npy'
        np.save(os.path.join(save_path,save_name), results)     
    




# In[event with backup energy]

from scipy.sparse import coo_matrix
def reconstruct_2d_array(index_value_pairs, num_rows, num_cols):
    rows, cols, values = zip(*index_value_pairs)
    return coo_matrix((values, (rows, cols)), shape=(num_rows, num_cols)).toarray()
    

import calendar
def get_annual_mean_potential():
    
    yearly_data = np.full((36, 600, 1440), np.nan, dtype=np.float32)
    
    for year in range(1986, 2022):  
        pathToRaster = 'G:\\global_photovoltatic\\yearly_daily_power_by_ERA5_pvlib\\power_sum_' + str(year) + '.tif'
        raster = xr.open_dataset(pathToRaster)
        
        raster.coords['x'] = (raster.coords['x'] + 180) % 360 - 180
        raster = raster.sortby(raster.x) 
        raster = raster.sortby(raster.y, ascending = False) 
        
        if calendar.isleap(year):
            days = 366
        else:
            days = 365
            
        yearly_data[year-1986,:,:] = raster['band_data'].values[0,:600,:] / days
        
    yearly_data[np.isnan(yearly_data)] = 0
    mean_pv = np.mean(yearly_data, axis=0)
    return mean_pv

    

def get_backup_level():
    mean_pv = get_annual_mean_potential()
    
    #
    # Read the low production record and calculate the backup 
    for year in range(1986, 2022):
        
        save_path = os.path.join(r'G:\global_photovoltatic\low_output') 
        save_name = 'event_' + str(year) + '.npy'
        df = np.load(os.path.join(save_path,save_name))         
        df = pd.DataFrame(df, columns = ['lon_index','lat_index',
                                         'thrs_10th','thrs_5th','thrs_1th','doy','power'])   
        
        # Exclude the part where the threshold is negative, which is usually due to the season, such as winter
        df = df[df['thrs_10th'] > 0]
        df['power'][df['power'] < 0] = 0
        df['diff'] = df['thrs_10th'] - df['power']
        cum_loss = df.groupby(['lon_index', 'lat_index']).sum()[['diff', 'thrs_10th']].reset_index()
        num_days = df.groupby(['lon_index', 'lat_index']).count()['doy'].reset_index()
        cum_loss = cum_loss.merge(num_days, on = ['lon_index', 'lat_index'])
 
        
        if year == 1986:
            cum_diff = cum_loss.copy()
        else:
            cum_diff = pd.concat((cum_diff, cum_loss))
        print(year)
        
    cum_loss = cum_diff.groupby(['lon_index', 'lat_index']).sum().reset_index()
    cum_loss['diff_daily_ave'] = cum_loss['diff'] / cum_loss['doy']
    cum_loss['daily_potential'] = mean_pv[np.array(cum_loss['lat_index']).astype(int), np.array(cum_loss['lon_index']).astype(int)]  
    
    cum_loss.to_csv(r'G:\global_photovoltatic\adaptation\backup\cum_loss.csv')
    

        
        
def main_1():    
    cum_loss = pd.read_csv(r'G:\global_photovoltatic\adaptation\backup\cum_loss.csv', index_col=0)
    cum_loss = cum_loss.rename({'doy':'lp_days'}, axis = 'columns')
    #
    # Read the low production record and recalculate the low-production event     
    for year in range(1986, 2022):    
        save_path = os.path.join(r'G:\global_photovoltatic\low_output') 
        save_name = 'event_' + str(year) + '.npy'
        df = np.load(os.path.join(save_path,save_name))         
        df = pd.DataFrame(df, columns = ['lon_index','lat_index', 
                                         'thrs_10th','thrs_5th','thrs_1th','doy','power'])   
        
        df = df.merge(cum_loss, on = ['lon_index', 'lat_index'])
        
        for factor in [0.5, 0.8, 1, 1.2, 1.5]:
            df['new_power_' + str(factor)] = df['power'] + df['diff_daily_ave'] * factor

  
        df = df[ (df['new_power_0.5']<df['thrs_10th']) | (df['new_power_0.5']<0) ]
        
        # ['lon_index', 'lat_index', 'thrs_10th', 'thrs_5th', 'thrs_1th', 'doy_x',
        #  'power', 'diff', 'doy_y'(low_output_days), 'diff_daily_ave', 'daily_potential',
        #  'new_power_0.5', 'new_power_0.8', 'new_power_1', 'new_power_1.2',
        #  'new_power_1.5']
        results = np.array(df, dtype=np.float16) 
            
        save_path = os.path.join(r'G:\global_photovoltatic\adaptation\backup\low_output') 
        save_name = 'event_' + str(year) + '.npy'
        np.save(os.path.join(save_path,save_name), results)     






