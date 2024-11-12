# -*- coding: utf-8 -*-
"""
Created on Sat May 18 18:17:10 2024

@author: wangx
"""


import os
os.chdir(r'G:\global_photovoltatic')

import numpy as np
import xarray as xr
import datetime
import warnings



def main():
    #
    # a 15-day rolling window    
    file_path = os.path.join(r'G:\global_photovoltatic\daily_power_by_ERA5_pvlib') 
    
    lats = np.arange(90, -60, -0.25) - 0.25/2
    lons = np.arange(-180, 180, 0.25) + 0.25/2
    
    thresholds = xr.DataArray(np.zeros((3, 365, 600, 1440),dtype=np.float32),
                                dims = ['thrs', 'doy', 'lat', 'lon'],
                                coords={'thrs':['1th', '5th', '10th'], 'doy': range(1, 366), 'lat': lats, 'lon': lons })    
    
    for block in range(2):
        
        daily_data = []
        t = datetime.datetime.now()
        
        for year in range(1986,2022):
            file_name = 'power_daily_' + str(year) + '.nc'
            data = xr.open_dataset( os.path.join(file_path, file_name) )
            
            if block == 0:
                data = data.sel(lon=slice(-180,0))
            if block == 1:
                data = data.sel(lon=slice(0,180))

            data['time'] = data['time'].dt.dayofyear
            daily_data.append(data)
            
        daily_data = xr.concat(daily_data, dim='time')       
        print('finish: concat', datetime.datetime.now() - t)    
        
 
        for doy in range(1, 366):
            t = datetime.datetime.now()
            #
            # Get the window of 7-8 days before and after the current dayofyear, a total of 15 days
            window_days = [(doy - i) % 365 +1 for i in range(-7, 8)]             
            window_data = daily_data.sel(time = daily_data['time'].isin(window_days))['band_data'].values
            
            if block == 0:
                thresholds[0, doy-1, :, :720] = np.percentile(window_data, 1, axis = 0)
                thresholds[1, doy-1, :, :720] = np.percentile(window_data, 5, axis = 0)
                thresholds[2, doy-1, :, :720] = np.percentile(window_data, 10, axis = 0)
                print(doy, datetime.datetime.now() - t)  
                
            if block == 1:
                thresholds[0, doy-1, :, 720:] = np.percentile(window_data, 1, axis = 0)
                thresholds[1, doy-1, :, 720:] = np.percentile(window_data, 5, axis = 0)
                thresholds[2, doy-1, :, 720:] = np.percentile(window_data, 10, axis = 0)
                print(doy, datetime.datetime.now() - t)   
                
        del window_data, daily_data
                
    save_path = r'G:\global_photovoltatic\threshold'                
    file_name =  'threshold.nc'
    thresholds.to_netcdf( os.path.join(save_path, file_name)
                         ,encoding={'__xarray_dataarray_variable__': {'zlib': True, 'complevel': 6}}) 



        