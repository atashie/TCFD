# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 14:27:11 2022

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
data_nc = Dataset('J:/Cai_data/TCFD/Flash Floods/clm45_gfdl-esm2m_ewembi_rcp26_2005soc_co2_maxdis_global_monthly_2006_2099.nc4')
    # accessing the variables
origLon = data_nc.variables['lon']
origLat = data_nc.variables['lat']
maxdis = data_nc.variables['maxdis']       # note time is in 3 dimensions, lon lat and time



fn = 'j:/Cai_data/TCFD/Flash Floods/futureFloodFreq_RCP26et60_2020sto2060s_CLM45_18apr2022.nc'

ds = nc.Dataset(fn,'w',format='NETCDF4')

stYr = ds.createDimension('stYr', 5)
HSP = ds.createDimension('HSP', 6)
RCP = ds.createDimension('RCP', 2)
lon = ds.createDimension('lon', len(origLon[:]))
lat = ds.createDimension('lat', len(origLat[:]))

stYrs = ds.createVariable('stYr', 'f4', ('stYr',))
HSPs = ds.createVariable('HSP', 'f4', ('HSP',))
RCPs = ds.createVariable('RCP', 'f4', ('RCP',))
lons = ds.createVariable('lon', 'f4', ('lon',))
lats = ds.createVariable('lat', 'f4', ('lat',))

returnPeriod = ds.createVariable('returnPeriod', 'f4', ('stYr','HSP', 'RCP', 'lat','lon',))
returnPeriod.units = 'Years'

stYrs[:] = np.array([2020,2030,2040,2050,2060])
HSPs[:] = np.array([5,10,20,30,50,100])
RCPs[:] = np.array([26,60])
lons[:] = np.array(origLon[:])    # np.arange(-179.75,179.76,0.5)#
lats[:] = np.array(origLat[:]) # this reverses the array, which appears high to low np.arange(-89.75,89.76,0.5)# 


##########################################
### now reading in the data

  ### RCP 26, 2029
# reading in csvs for compilation
data_cs_5 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_5yrflod_CLM_2029.csv')
data_dr_5 = data_cs_5.iloc[:,1:]
data_cs_10 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_10yrflod_CLM_2029.csv')
data_dr_10 = data_cs_10.iloc[:,1:]
data_cs_20 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_20yrflod_CLM_2029.csv')
data_dr_20 = data_cs_20.iloc[:,1:]
data_cs_30 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_30yrflod_CLM_2029.csv')
data_dr_30 = data_cs_30.iloc[:,1:]
data_cs_50 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_50yrflod_CLM_2029.csv')
data_dr_50 = data_cs_50.iloc[:,1:]
data_cs_100 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_100yrflod_CLM_2029.csv')
data_dr_100 = data_cs_100.iloc[:,1:]

csvDataAr_5 = pd.DataFrame(data_dr_5).to_numpy()
csvDataAr_10 = pd.DataFrame(data_dr_10).to_numpy()
csvDataAr_20 = pd.DataFrame(data_dr_20).to_numpy()
csvDataAr_30 = pd.DataFrame(data_dr_30).to_numpy()
csvDataAr_50 = pd.DataFrame(data_dr_50).to_numpy()
csvDataAr_100 = pd.DataFrame(data_dr_100).to_numpy()
#csvDataArRshp = np.reshape(csvDataAr, (len(origLat),len(origLon)))

returnPeriod[0, 0, 0, : , :] =csvDataAr_5     # ('stYr','HSP', 'RCP', 'lat','lon',)
returnPeriod[0, 1, 0, : , :] =csvDataAr_10
returnPeriod[0, 2, 0, : , :] =csvDataAr_20
returnPeriod[0, 3, 0, : , :] =csvDataAr_30
returnPeriod[0, 4, 0, : , :] =csvDataAr_50
returnPeriod[0, 5, 0, : , :] =csvDataAr_100


  ###  RCP 26, 3039
# reading in csvs for compilation
data_cs_5 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_5yrflod_CLM_3039.csv')
data_dr_5 = data_cs_5.iloc[:,1:]
data_cs_10 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_10yrflod_CLM_3039.csv')
data_dr_10 = data_cs_10.iloc[:,1:]
data_cs_20 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_20yrflod_CLM_3039.csv')
data_dr_20 = data_cs_20.iloc[:,1:]
data_cs_30 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_30yrflod_CLM_3039.csv')
data_dr_30 = data_cs_30.iloc[:,1:]
data_cs_50 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_50yrflod_CLM_3039.csv')
data_dr_50 = data_cs_50.iloc[:,1:]
data_cs_100 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_100yrflod_CLM_3039.csv')
data_dr_100 = data_cs_100.iloc[:,1:]

csvDataAr_5 = pd.DataFrame(data_dr_5).to_numpy()
csvDataAr_10 = pd.DataFrame(data_dr_10).to_numpy()
csvDataAr_20 = pd.DataFrame(data_dr_20).to_numpy()
csvDataAr_30 = pd.DataFrame(data_dr_30).to_numpy()
csvDataAr_50 = pd.DataFrame(data_dr_50).to_numpy()
csvDataAr_100 = pd.DataFrame(data_dr_100).to_numpy()
#csvDataArRshp = np.reshape(csvDataAr, (len(origLat),len(origLon)))

returnPeriod[1, 0, 0, : , :] =csvDataAr_5     # ('stYr','HSP', 'RCP', 'lat','lon',)
returnPeriod[1, 1, 0, : , :] =csvDataAr_10
returnPeriod[1, 2, 0, : , :] =csvDataAr_20
returnPeriod[1, 3, 0, : , :] =csvDataAr_30
returnPeriod[1, 4, 0, : , :] =csvDataAr_50
returnPeriod[1, 5, 0, : , :] =csvDataAr_100


  ###  RCP 26, 4049
# reading in csvs for compilation
data_cs_5 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_5yrflod_CLM_4049.csv')
data_dr_5 = data_cs_5.iloc[:,1:]
data_cs_10 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_10yrflod_CLM_4049.csv')
data_dr_10 = data_cs_10.iloc[:,1:]
data_cs_20 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_20yrflod_CLM_4049.csv')
data_dr_20 = data_cs_20.iloc[:,1:]
data_cs_30 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_30yrflod_CLM_4049.csv')
data_dr_30 = data_cs_30.iloc[:,1:]
data_cs_50 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_50yrflod_CLM_4049.csv')
data_dr_50 = data_cs_50.iloc[:,1:]
data_cs_100 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_100yrflod_CLM_4049.csv')
data_dr_100 = data_cs_100.iloc[:,1:]

csvDataAr_5 = pd.DataFrame(data_dr_5).to_numpy()
csvDataAr_10 = pd.DataFrame(data_dr_10).to_numpy()
csvDataAr_20 = pd.DataFrame(data_dr_20).to_numpy()
csvDataAr_30 = pd.DataFrame(data_dr_30).to_numpy()
csvDataAr_50 = pd.DataFrame(data_dr_50).to_numpy()
csvDataAr_100 = pd.DataFrame(data_dr_100).to_numpy()
#csvDataArRshp = np.reshape(csvDataAr, (len(origLat),len(origLon)))

returnPeriod[2, 0, 0, : , :] =csvDataAr_5     # ('stYr','HSP', 'RCP', 'lat','lon',)
returnPeriod[2, 1, 0, : , :] =csvDataAr_10
returnPeriod[2, 2, 0, : , :] =csvDataAr_20
returnPeriod[2, 3, 0, : , :] =csvDataAr_30
returnPeriod[2, 4, 0, : , :] =csvDataAr_50
returnPeriod[2, 5, 0, : , :] =csvDataAr_100


  ###  RCP 26, 5059
# reading in csvs for compilation
data_cs_5 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_5yrflod_CLM_5059.csv')
data_dr_5 = data_cs_5.iloc[:,1:]
data_cs_10 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_10yrflod_CLM_5059.csv')
data_dr_10 = data_cs_10.iloc[:,1:]
data_cs_20 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_20yrflod_CLM_5059.csv')
data_dr_20 = data_cs_20.iloc[:,1:]
data_cs_30 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_30yrflod_CLM_5059.csv')
data_dr_30 = data_cs_30.iloc[:,1:]
data_cs_50 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_50yrflod_CLM_5059.csv')
data_dr_50 = data_cs_50.iloc[:,1:]
data_cs_100 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_100yrflod_CLM_5059.csv')
data_dr_100 = data_cs_100.iloc[:,1:]

csvDataAr_5 = pd.DataFrame(data_dr_5).to_numpy()
csvDataAr_10 = pd.DataFrame(data_dr_10).to_numpy()
csvDataAr_20 = pd.DataFrame(data_dr_20).to_numpy()
csvDataAr_30 = pd.DataFrame(data_dr_30).to_numpy()
csvDataAr_50 = pd.DataFrame(data_dr_50).to_numpy()
csvDataAr_100 = pd.DataFrame(data_dr_100).to_numpy()
#csvDataArRshp = np.reshape(csvDataAr, (len(origLat),len(origLon)))

returnPeriod[3, 0, 0, : , :] =csvDataAr_5     # ('stYr','HSP', 'RCP', 'lat','lon',)
returnPeriod[3, 1, 0, : , :] =csvDataAr_10
returnPeriod[3, 2, 0, : , :] =csvDataAr_20
returnPeriod[3, 3, 0, : , :] =csvDataAr_30
returnPeriod[3, 4, 0, : , :] =csvDataAr_50
returnPeriod[3, 5, 0, : , :] =csvDataAr_100


  ###  RCP 26, 6069
# reading in csvs for compilation
data_cs_5 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_5yrflod_CLM_6069.csv')
data_dr_5 = data_cs_5.iloc[:,1:]
data_cs_10 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_10yrflod_CLM_6069.csv')
data_dr_10 = data_cs_10.iloc[:,1:]
data_cs_20 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_20yrflod_CLM_6069.csv')
data_dr_20 = data_cs_20.iloc[:,1:]
data_cs_30 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_30yrflod_CLM_6069.csv')
data_dr_30 = data_cs_30.iloc[:,1:]
data_cs_50 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_50yrflod_CLM_6069.csv')
data_dr_50 = data_cs_50.iloc[:,1:]
data_cs_100 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp26_100yrflod_CLM_6069.csv')
data_dr_100 = data_cs_100.iloc[:,1:]

csvDataAr_5 = pd.DataFrame(data_dr_5).to_numpy()
csvDataAr_10 = pd.DataFrame(data_dr_10).to_numpy()
csvDataAr_20 = pd.DataFrame(data_dr_20).to_numpy()
csvDataAr_30 = pd.DataFrame(data_dr_30).to_numpy()
csvDataAr_50 = pd.DataFrame(data_dr_50).to_numpy()
csvDataAr_100 = pd.DataFrame(data_dr_100).to_numpy()
#csvDataArRshp = np.reshape(csvDataAr, (len(origLat),len(origLon)))

returnPeriod[4, 0, 0, : , :] =csvDataAr_5     # ('stYr','HSP', 'RCP', 'lat','lon',)
returnPeriod[4, 1, 0, : , :] =csvDataAr_10
returnPeriod[4, 2, 0, : , :] =csvDataAr_20
returnPeriod[4, 3, 0, : , :] =csvDataAr_30
returnPeriod[4, 4, 0, : , :] =csvDataAr_50
returnPeriod[4, 5, 0, : , :] =csvDataAr_100




  ### RCP 60, 2029
# reading in csvs for compilation
data_cs_5 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_5yrflod_CLM_2029.csv')
data_dr_5 = data_cs_5.iloc[:,1:]
data_cs_10 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_10yrflod_CLM_2029.csv')
data_dr_10 = data_cs_10.iloc[:,1:]
data_cs_20 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_20yrflod_CLM_2029.csv')
data_dr_20 = data_cs_20.iloc[:,1:]
data_cs_30 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_30yrflod_CLM_2029.csv')
data_dr_30 = data_cs_30.iloc[:,1:]
data_cs_50 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_50yrflod_CLM_2029.csv')
data_dr_50 = data_cs_50.iloc[:,1:]
data_cs_100 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_100yrflod_CLM_2029.csv')
data_dr_100 = data_cs_100.iloc[:,1:]

csvDataAr_5 = pd.DataFrame(data_dr_5).to_numpy()
csvDataAr_10 = pd.DataFrame(data_dr_10).to_numpy()
csvDataAr_20 = pd.DataFrame(data_dr_20).to_numpy()
csvDataAr_30 = pd.DataFrame(data_dr_30).to_numpy()
csvDataAr_50 = pd.DataFrame(data_dr_50).to_numpy()
csvDataAr_100 = pd.DataFrame(data_dr_100).to_numpy()
#csvDataArRshp = np.reshape(csvDataAr, (len(origLat),len(origLon)))

returnPeriod[0, 0, 1, : , :] =csvDataAr_5     # ('stYr','HSP', 'RCP', 'lat','lon',)
returnPeriod[0, 1, 1, : , :] =csvDataAr_10
returnPeriod[0, 2, 1, : , :] =csvDataAr_20
returnPeriod[0, 3, 1, : , :] =csvDataAr_30
returnPeriod[0, 4, 1, : , :] =csvDataAr_50
returnPeriod[0, 5, 1, : , :] =csvDataAr_100


  ###  RCP 60, 3039
# reading in csvs for compilation
data_cs_5 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_5yrflod_CLM_3039.csv')
data_dr_5 = data_cs_5.iloc[:,1:]
data_cs_10 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_10yrflod_CLM_3039.csv')
data_dr_10 = data_cs_10.iloc[:,1:]
data_cs_20 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_20yrflod_CLM_3039.csv')
data_dr_20 = data_cs_20.iloc[:,1:]
data_cs_30 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_30yrflod_CLM_3039.csv')
data_dr_30 = data_cs_30.iloc[:,1:]
data_cs_50 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_50yrflod_CLM_3039.csv')
data_dr_50 = data_cs_50.iloc[:,1:]
data_cs_100 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_100yrflod_CLM_3039.csv')
data_dr_100 = data_cs_100.iloc[:,1:]

csvDataAr_5 = pd.DataFrame(data_dr_5).to_numpy()
csvDataAr_10 = pd.DataFrame(data_dr_10).to_numpy()
csvDataAr_20 = pd.DataFrame(data_dr_20).to_numpy()
csvDataAr_30 = pd.DataFrame(data_dr_30).to_numpy()
csvDataAr_50 = pd.DataFrame(data_dr_50).to_numpy()
csvDataAr_100 = pd.DataFrame(data_dr_100).to_numpy()
#csvDataArRshp = np.reshape(csvDataAr, (len(origLat),len(origLon)))

returnPeriod[1, 0, 1, : , :] =csvDataAr_5     # ('stYr','HSP', 'RCP', 'lat','lon',)
returnPeriod[1, 1, 1, : , :] =csvDataAr_10
returnPeriod[1, 2, 1, : , :] =csvDataAr_20
returnPeriod[1, 3, 1, : , :] =csvDataAr_30
returnPeriod[1, 4, 1, : , :] =csvDataAr_50
returnPeriod[1, 5, 1, : , :] =csvDataAr_100


  ###  RCP 60, 4049
# reading in csvs for compilation
data_cs_5 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_5yrflod_CLM_4049.csv')
data_dr_5 = data_cs_5.iloc[:,1:]
data_cs_10 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_10yrflod_CLM_4049.csv')
data_dr_10 = data_cs_10.iloc[:,1:]
data_cs_20 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_20yrflod_CLM_4049.csv')
data_dr_20 = data_cs_20.iloc[:,1:]
data_cs_30 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_30yrflod_CLM_4049.csv')
data_dr_30 = data_cs_30.iloc[:,1:]
data_cs_50 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_50yrflod_CLM_4049.csv')
data_dr_50 = data_cs_50.iloc[:,1:]
data_cs_100 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_100yrflod_CLM_4049.csv')
data_dr_100 = data_cs_100.iloc[:,1:]

csvDataAr_5 = pd.DataFrame(data_dr_5).to_numpy()
csvDataAr_10 = pd.DataFrame(data_dr_10).to_numpy()
csvDataAr_20 = pd.DataFrame(data_dr_20).to_numpy()
csvDataAr_30 = pd.DataFrame(data_dr_30).to_numpy()
csvDataAr_50 = pd.DataFrame(data_dr_50).to_numpy()
csvDataAr_100 = pd.DataFrame(data_dr_100).to_numpy()
#csvDataArRshp = np.reshape(csvDataAr, (len(origLat),len(origLon)))

returnPeriod[2, 0, 1, : , :] =csvDataAr_5     # ('stYr','HSP', 'RCP', 'lat','lon',)
returnPeriod[2, 1, 1, : , :] =csvDataAr_10
returnPeriod[2, 2, 1, : , :] =csvDataAr_20
returnPeriod[2, 3, 1, : , :] =csvDataAr_30
returnPeriod[2, 4, 1, : , :] =csvDataAr_50
returnPeriod[2, 5, 1, : , :] =csvDataAr_100


  ###  RCP 60, 5059
# reading in csvs for compilation
data_cs_5 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_5yrflod_CLM_5059.csv')
data_dr_5 = data_cs_5.iloc[:,1:]
data_cs_10 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_10yrflod_CLM_5059.csv')
data_dr_10 = data_cs_10.iloc[:,1:]
data_cs_20 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_20yrflod_CLM_5059.csv')
data_dr_20 = data_cs_20.iloc[:,1:]
data_cs_30 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_30yrflod_CLM_5059.csv')
data_dr_30 = data_cs_30.iloc[:,1:]
data_cs_50 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_50yrflod_CLM_5059.csv')
data_dr_50 = data_cs_50.iloc[:,1:]
data_cs_100 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_100yrflod_CLM_5059.csv')
data_dr_100 = data_cs_100.iloc[:,1:]

csvDataAr_5 = pd.DataFrame(data_dr_5).to_numpy()
csvDataAr_10 = pd.DataFrame(data_dr_10).to_numpy()
csvDataAr_20 = pd.DataFrame(data_dr_20).to_numpy()
csvDataAr_30 = pd.DataFrame(data_dr_30).to_numpy()
csvDataAr_50 = pd.DataFrame(data_dr_50).to_numpy()
csvDataAr_100 = pd.DataFrame(data_dr_100).to_numpy()
#csvDataArRshp = np.reshape(csvDataAr, (len(origLat),len(origLon)))

returnPeriod[3, 0, 1, : , :] =csvDataAr_5     # ('stYr','HSP', 'RCP', 'lat','lon',)
returnPeriod[3, 1, 1, : , :] =csvDataAr_10
returnPeriod[3, 2, 1, : , :] =csvDataAr_20
returnPeriod[3, 3, 1, : , :] =csvDataAr_30
returnPeriod[3, 4, 1, : , :] =csvDataAr_50
returnPeriod[3, 5, 1, : , :] =csvDataAr_100


  ###  RCP 60, 6069
# reading in csvs for compilation
data_cs_5 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_5yrflod_CLM_6069.csv')
data_dr_5 = data_cs_5.iloc[:,1:]
data_cs_10 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_10yrflod_CLM_6069.csv')
data_dr_10 = data_cs_10.iloc[:,1:]
data_cs_20 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_20yrflod_CLM_6069.csv')
data_dr_20 = data_cs_20.iloc[:,1:]
data_cs_30 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_30yrflod_CLM_6069.csv')
data_dr_30 = data_cs_30.iloc[:,1:]
data_cs_50 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_50yrflod_CLM_6069.csv')
data_dr_50 = data_cs_50.iloc[:,1:]
data_cs_100 = pd.read_csv('J:/Cai_data/TCFD/Flash Floods/rcp60_100yrflod_CLM_6069.csv')
data_dr_100 = data_cs_100.iloc[:,1:]

csvDataAr_5 = pd.DataFrame(data_dr_5).to_numpy()
csvDataAr_10 = pd.DataFrame(data_dr_10).to_numpy()
csvDataAr_20 = pd.DataFrame(data_dr_20).to_numpy()
csvDataAr_30 = pd.DataFrame(data_dr_30).to_numpy()
csvDataAr_50 = pd.DataFrame(data_dr_50).to_numpy()
csvDataAr_100 = pd.DataFrame(data_dr_100).to_numpy()
#csvDataArRshp = np.reshape(csvDataAr, (len(origLat),len(origLon)))

returnPeriod[4, 0, 1, : , :] =csvDataAr_5     # ('stYr','HSP', 'RCP', 'lat','lon',)
returnPeriod[4, 1, 1, : , :] =csvDataAr_10
returnPeriod[4, 2, 1, : , :] =csvDataAr_20
returnPeriod[4, 3, 1, : , :] =csvDataAr_30
returnPeriod[4, 4, 1, : , :] =csvDataAr_50
returnPeriod[4, 5, 1, : , :] =csvDataAr_100


### finish reading in the data
##########################################



ds.close()


































#######################################################
## netCDF4 data ex from APHRODITE's Water Resources




#########################################################
## reading in data using netCDf4
from netCDF4 import Dataset
import pandas as pd

    # reading in the dataset
data = Dataset(fn)
    # checking names of variables
print(data.variables.keys())
    # accessing the variables
lon = data.variables['lon']
lat = data.variables['lat']
stYr = data.variables['stYr']       
HSP = data.variables['HSP']       
RCP = data.variables['RCP']       
returnPeriod = data.variables['returnPeriod']

    # accessing data from the variables
stYr_data = data.variables['stYr'][:]
lon_data = data.variables['lon'][:]
lat_data = data.variables['lat'][:]
HSP_data = data.variables['HSP'][:]
RCP_data = data.variables['RCP'][:]
returnPeriod_data = data.variables['returnPeriod'][:]

    # target coordinates
lat_katmandu = 27.697817
lon_katmandu = 85.329806

  
##########################################################
## plotting netCDF data
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

mp= Basemap(projection ='merc',
            llcrnrlon = -180,
            llcrnrlat = -70,
            urcrnrlon = 180,
            urcrnrlat = 85,
            resolution = 'i')

lon,lat = np.meshgrid(lon_data, lat_data)
x,y = mp(lon,lat)

c_scheme = mp.pcolor(x,y, np.squeeze(returnPeriod[4,0,1,:,:]), cmap='jet')
mp.drawcoastlines()
mp.drawstates()
mp.drawcountries()
cbar = mp.colorbar(c_scheme, location = 'right', pad='3%')
plt.title('Average Temperature on 01JAN1980')
plt.show()

##########################################################
## plotting netCDF data as a timelapse
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

mp= Basemap(projection ='merc',
            llcrnrlon = 42,
            llcrnrlat = -2,
            urcrnrlon = 106,
            urcrnrlat = 38,
            resolution = 'i')

lon,lat = np.meshgrid(lon_data, lat_data)
x,y = mp(lon,lat)

days = np.arange(0,365)

for i in days:
    c_scheme = mp.pcolor(x,y, np.squeeze(tave[i,:,:]), cmap='jet')
    mp.drawcoastlines()
    mp.drawstates()
    mp.drawcountries()

    cbar = mp.colorbar(c_scheme, location = 'right', pad='3%')
    this_day = i+1

    plt.title('Average Temperature: Day ' + str(this_day) + ' of 1980')
    plt.clim(-40,40)
    plt.savefig(save_path + 'jpgs\\' +  str(this_day) + '.jpg')
    plt.show()
    plt.clf()   # clearing the figure


import PIL

image_frames = []

for k in (days+1):
    new_frame = PIL.Image.open(save_path + 'jpgs\\' +  str(k) + '.jpg')
    image_frames.append(new_frame)  # storing as a list

image_frames[0].save(save_path + 'jpgs\\' + "temperature_timelapse.gif", format= 'GIF',
                     append_images = image_frames[1:],
                     save_all = True, duration = 300, loop = 0)



##########################################################
## extracting data from multiple netCDFs into a csv
import glob
from netCDF4 import Dataset

directory_path = 'C:\\Users\\arik\\Documents\\Other Research\\climai\\aphro\\'
all_years = []


for this_file in glob.glob(directory_path + '*.nc'):  # wow rad function wiwh this existed (or Iknew about) in R
    print(this_file)
    the_data = Dataset(this_file, 'r')
    the_time = the_data.variables['time']
    the_year = the_time.units[14:18]
    all_years.append(the_year)

start_year = min(all_years)
end_year = max(all_years)
date_range = pd.date_range(start = str(start_year) + '-01-01',
                           end = str(end_year) + '-12-31',
                           freq = 'D')
new_df = pd.DataFrame(None, columns = ['Temp'], index = date_range)

# defining the lat lon for location of interest
lat_ktmndu = 27.698
lon_ktmndu = 85.330

all_years.sort()    # ensuring that the years are sorted

for yr in all_years:
    # reading in the data
    data = Dataset(directory_path + 'APHRO_MA_TAVE_025deg_V1808.' + str(yr) + '.nc', 'r')
    
    # storing lat lon data of netCDF file into variables
    lat = data.variables['lat'][:]
    lon = data.variables['lon'][:]
    
    # nearest neighbor
    cls_lat = abs(lat - lat_ktmndu)
    cls_lon = abs(lon - lon_ktmndu)
    
    #  identify the index of the min value for lat and lon
    min_index_lat = cls_lat.argmin()
    max_index_lon = cls_lon.argmin()
    
    temp = data.variables['tave']














##############################################
## xarray as netCDF reader
import xarray as xr
this_path = r'C:\Users\arik\Documents\Other Research\climai'
NOA_dat = r"\NOAA_GFS_2m_temperature_20201201Z0000.nc"
ERA_dat = r"\ECMWF_ERA5_2m_temperature_20201201Z0000.nc"

ds_ERA = xr.open_dataset(filename_or_obj = this_path + ERA_dat) # temp in K
ds_NOA = xr.open_dataset(filename_or_obj = this_path + NOA_dat) # temp in K

print(ds_ERA)   # data are global in lat lon; crop for analysis
min_lon = 360-95
max_lon = 360-60
min_lat = 23.5
max_lat = 47
cropped = (ds_ERA.latitude > min_lat) & (ds_ERA.latitude < max_lat) & (ds_ERA.longitude > min_lon) & (ds_ERA.longitude < max_lon)

ds_ERA_crp = ds_ERA.where(cropped, drop=True)     #cropping to target region
lat_ERA = ds_ERA_crp['latitude'].data
lon_ERA = ds_ERA_crp['longitude'].data
tmp_ERA = ds_ERA_crp['temperature'].data - 273.15 # convert temp to C

print(ds_NOA)   # lower (1/2) resolution than ERA
ds_NOA_intrp = ds_NOA.interp_like(ds_ERA, method='linear')   # linear interp to coordinates of ERA
ds_NOA_intrp_crp = ds_NOA_intrp.where(cropped, drop=True)
tmp_NOA = ds_NOA_intrp_crp['temperature'].data  - 273.15 # convert temp to C
#lat_NOA = ds_NOA['latitude'].data # saving for comparison with interpolated NOA array
#lon_NOA = ds_NOA['longitude'].data

    # calculating error as simple difference
ds_error = tmp_NOA - tmp_ERA
print(ds_error)



##############################################
## data plotting
    # import basemap
import matplotlib.colors as colors
import matplotlib.pyplot as plt # note check dependency errors, potentially reinstall via anaconda console admin
from mpl_toolkits.basemap import Basemap

    # initializing basemaps and coordinates  
bsmp = Basemap(projection = 'merc',
               llcrnrlon = min_lon +.25,
               llcrnrlat = min_lat +.25,
               urcrnrlon = max_lon -.25,
               urcrnrlat = max_lat -.25,
               resolution = 'i')

lon,lat = np.meshgrid(lon_ERA, lat_ERA)
x,y = bsmp(lon, lat)


    # plotting ERA data
ds_ERA_crp.temperature.mean() - 273.15
ds_NOA_intrp_crp.temperature.mean() - 273.15
ds_error.mean()
print(round(np.mean(ds_NOA_intrp_crp.temperature).values + 1, 1))


col_schm = bsmp.pcolor(x, y,
                       np.squeeze(tmp_ERA),
                       cmap = 'jet')
bsmp.drawcoastlines(color='black', linewidth=1.2)
bsmp.drawcountries(color='black', linewidth=1.2)
bsmp.drawrivers(color='blue',linewidth=0.5)
bsmp.drawstates(color='black', linewidth=1)
cbar = bsmp.colorbar(col_schm, location = 'right', pad = '5%')

plt.title('ERA5 Average Temperature (C) on 01-DEC-2020')
plt.figtext(0.4, .07,#.25, 
            "mean="   + str(round(np.mean(tmp_ERA), 1)) + 'C, ' +
            'std=' + str(round(np.std(tmp_ERA), 1)) + 'C')
plt.savefig(r"C:\Users\arik\Documents\Other Research\climai\ERA5_avg_T.png", dpi=1000, format='png')
plt.show()


    # plotting NOAA data
col_schm = bsmp.pcolor(x, y,
                       np.squeeze(tmp_NOA),
                       cmap = 'jet')
bsmp.drawcoastlines(color='black', linewidth=1.2)
bsmp.drawcountries(color='black', linewidth=1.2)
bsmp.drawrivers(color='blue',linewidth=0.5)
bsmp.drawstates(color='black', linewidth=1)
cbar = bsmp.colorbar(col_schm, location = 'right', pad = '5%')

plt.title('NOAA GFS Average Temperature (C) on 01-DEC-2020')
plt.figtext(0.4, .07,#.25, 
            "mean="   + str(round(np.mean(tmp_NOA), 1)) + 'C, ' +
            'std=' + str(round(np.std(tmp_NOA), 1)) + 'C')
plt.savefig(r"C:\Users\arik\Documents\Other Research\climai\NOAA_avg_T.png", dpi=1000, format='png')
plt.show()


    # plotting error
col_schm = bsmp.pcolor(x, y,
                       np.squeeze(ds_error),
                       cmap = 'seismic',
                       norm=colors.CenteredNorm())
bsmp.drawcoastlines(color='black', linewidth=1.2)
bsmp.drawcountries(color='black', linewidth=1.2)
bsmp.drawrivers(color='blue',linewidth=0.5)
bsmp.drawstates(color='black', linewidth=1)
cbar = bsmp.colorbar(col_schm, location = 'right', pad = '5%')

plt.title('Prediction Temperature (C) Error on 01-DEC-2020')
plt.figtext(0.4, .07,#.25, 
            "mean=" + str(round(np.mean(ds_error), 1)) + 'C, ' +
            'std=' + str(round(np.std(ds_error), 1)) + 'C')
plt.savefig(r"C:\Users\arik\Documents\Other Research\climai\error_T.png", dpi=1000, format='png')
plt.show()

