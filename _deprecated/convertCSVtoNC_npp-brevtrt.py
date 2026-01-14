# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:44:45 2022

@author: arik
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 11:48:20 2022

@author: arik
"""

import numpy as np
import numpy.ma as ma
import xarray as xr
import pandas as pd
from netCDF4 import Dataset
import netCDF4 as nc

 
# reading in original .nc for dimensions
   # reading in the dataset
data_nc = Dataset('J:/Cai_data/TCFD/BrdEvrgrnTrTropical/caraib_gfdl-esm2m_ewembi_rcp60_2005soc_co2_npp-brevtrt_global_annual_2006_2099.nc4')
    # accessing the variables
origLon = data_nc.variables['lon']
origLat = data_nc.variables['lat']
fireint = data_nc.variables['npp-brevtrt']       # note time is in 3 dimensions, lon lat and time



fn = 'j:/Cai_data/TCFD/BrdEvrgrnTrTropical/npp-brevtrt_RCP26et60_2020sto2090s_CARAIB_4May2022.nc'

ds = nc.Dataset(fn,'w',format='NETCDF4')

stYr = ds.createDimension('stYr', 8)
#HSP = ds.createDimension('HSP', 6)
RCP = ds.createDimension('RCP', 2)
lon = ds.createDimension('lon', len(origLon[:]))
lat = ds.createDimension('lat', len(origLat[:]))

stYrs = ds.createVariable('stYr', 'f4', ('stYr',))
#HSPs = ds.createVariable('HSP', 'f4', ('HSP',))
RCPs = ds.createVariable('RCP', 'f4', ('RCP',))
lons = ds.createVariable('lon', 'f4', ('lon',))
lats = ds.createVariable('lat', 'f4', ('lat',))

fireint = ds.createVariable('npp-brevtrt', 'f4', ('stYr', 'RCP', 'lat','lon',))
fireint.units = 'kg m-2 s-1'

stYrs[:] = np.array([2020,2030,2040,2050,2060,2070,2080,2090])
#HSPs[:] = np.array([5,10,20,30,50,100])
RCPs[:] = np.array([26,60])
lons[:] = np.array(origLon[:])    # np.arange(-179.75,179.76,0.5)#
lats[:] = np.array(origLat[:]) # this reverses the array, which appears high to low np.arange(-89.75,89.76,0.5)# 


##########################################
### now reading in the data

  ### RCP 26
# reading in csvs for compilation
data_cs_20 = pd.read_csv('J:/Cai_data/TCFD/BrdEvrgrnTrTropical/nppchng_CARAIB_rcp26_2029.csv')
data_dr_20 = data_cs_20.iloc[:,1:]
data_cs_30 = pd.read_csv('J:/Cai_data/TCFD/BrdEvrgrnTrTropical/nppchng_CARAIB_rcp26_3039.csv')
data_dr_30 = data_cs_30.iloc[:,1:]
data_cs_40 = pd.read_csv('J:/Cai_data/TCFD/BrdEvrgrnTrTropical/nppchng_CARAIB_rcp26_4049.csv')
data_dr_40 = data_cs_40.iloc[:,1:]
data_cs_50 = pd.read_csv('J:/Cai_data/TCFD/BrdEvrgrnTrTropical/nppchng_CARAIB_rcp26_5059.csv')
data_dr_50 = data_cs_50.iloc[:,1:]
data_cs_60 = pd.read_csv('J:/Cai_data/TCFD/BrdEvrgrnTrTropical/nppchng_CARAIB_rcp26_6069.csv')
data_dr_60 = data_cs_60.iloc[:,1:]
data_cs_70 = pd.read_csv('J:/Cai_data/TCFD/BrdEvrgrnTrTropical/nppchng_CARAIB_rcp26_7079.csv')
data_dr_70 = data_cs_70.iloc[:,1:]
data_cs_80 = pd.read_csv('J:/Cai_data/TCFD/BrdEvrgrnTrTropical/nppchng_CARAIB_rcp26_8089.csv')
data_dr_80 = data_cs_80.iloc[:,1:]
data_cs_90 = pd.read_csv('J:/Cai_data/TCFD/BrdEvrgrnTrTropical/nppchng_CARAIB_rcp26_9099.csv')
data_dr_90 = data_cs_90.iloc[:,1:]

csvDataAr_20 = pd.DataFrame(data_dr_20).to_numpy()
csvDataAr_30 = pd.DataFrame(data_dr_30).to_numpy()
csvDataAr_40 = pd.DataFrame(data_dr_40).to_numpy()
csvDataAr_50 = pd.DataFrame(data_dr_50).to_numpy()
csvDataAr_60 = pd.DataFrame(data_dr_60).to_numpy()
csvDataAr_70 = pd.DataFrame(data_dr_70).to_numpy()
csvDataAr_80 = pd.DataFrame(data_dr_80).to_numpy()
csvDataAr_90 = pd.DataFrame(data_dr_90).to_numpy()
#csvDataArRshp = np.reshape(csvDataAr, (len(origLat),len(origLon)))

fireint[0, 0, : , :] =csvDataAr_20     # ('stYr','RCP', 'lat','lon',)
fireint[1, 0, : , :] =csvDataAr_30
fireint[2, 0, : , :] =csvDataAr_40
fireint[3, 0, : , :] =csvDataAr_50
fireint[4, 0, : , :] =csvDataAr_60
fireint[5, 0, : , :] =csvDataAr_70
fireint[5, 0, : , :] =csvDataAr_80
fireint[5, 0, : , :] =csvDataAr_90



#### RCP 26
# reading in csvs for compilation
data_cs_20 = pd.read_csv('J:/Cai_data/TCFD/BrdEvrgrnTrTropical/nppchng_CARAIB_rcp60_2029.csv')
data_dr_20 = data_cs_20.iloc[:,1:]
data_cs_30 = pd.read_csv('J:/Cai_data/TCFD/BrdEvrgrnTrTropical/nppchng_CARAIB_rcp60_3039.csv')
data_dr_30 = data_cs_30.iloc[:,1:]
data_cs_40 = pd.read_csv('J:/Cai_data/TCFD/BrdEvrgrnTrTropical/nppchng_CARAIB_rcp60_4049.csv')
data_dr_40 = data_cs_40.iloc[:,1:]
data_cs_50 = pd.read_csv('J:/Cai_data/TCFD/BrdEvrgrnTrTropical/nppchng_CARAIB_rcp60_5059.csv')
data_dr_50 = data_cs_50.iloc[:,1:]
data_cs_60 = pd.read_csv('J:/Cai_data/TCFD/BrdEvrgrnTrTropical/nppchng_CARAIB_rcp60_6069.csv')
data_dr_60 = data_cs_60.iloc[:,1:]
data_cs_70 = pd.read_csv('J:/Cai_data/TCFD/BrdEvrgrnTrTropical/nppchng_CARAIB_rcp60_7079.csv')
data_dr_70 = data_cs_70.iloc[:,1:]
data_cs_80 = pd.read_csv('J:/Cai_data/TCFD/BrdEvrgrnTrTropical/nppchng_CARAIB_rcp60_8089.csv')
data_dr_80 = data_cs_80.iloc[:,1:]
data_cs_90 = pd.read_csv('J:/Cai_data/TCFD/BrdEvrgrnTrTropical/nppchng_CARAIB_rcp60_9099.csv')
data_dr_90 = data_cs_90.iloc[:,1:]

csvDataAr_20 = pd.DataFrame(data_dr_20).to_numpy()
csvDataAr_30 = pd.DataFrame(data_dr_30).to_numpy()
csvDataAr_40 = pd.DataFrame(data_dr_40).to_numpy()
csvDataAr_50 = pd.DataFrame(data_dr_50).to_numpy()
csvDataAr_60 = pd.DataFrame(data_dr_60).to_numpy()
csvDataAr_70 = pd.DataFrame(data_dr_70).to_numpy()
csvDataAr_80 = pd.DataFrame(data_dr_80).to_numpy()
csvDataAr_90 = pd.DataFrame(data_dr_90).to_numpy()
#csvDataArRshp = np.reshape(csvDataAr, (len(origLat),len(origLon)))

fireint[0, 1, : , :] =csvDataAr_20     # ('stYr','RCP', 'lat','lon',)
fireint[1, 1, : , :] =csvDataAr_30
fireint[2, 1, : , :] =csvDataAr_40
fireint[3, 1, : , :] =csvDataAr_50
fireint[4, 1, : , :] =csvDataAr_60
fireint[5, 1, : , :] =csvDataAr_70
fireint[5, 1, : , :] =csvDataAr_80
fireint[5, 1, : , :] =csvDataAr_90





### finish reading in the data
##########################################



ds.close()














