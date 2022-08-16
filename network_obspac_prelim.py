# -*- coding: utf-8 -*-
"""
Created on Wed May 18 19:20:16 2022

@author: jhask
"""
import xarray as xr 
import numpy as np 
import pandas as pd 
import sys 
import glob
import os 
sys.path.insert(0, 'C://Users//jhask//OneDrive//Documents//Research//Projects//MIT//pOrgNO3//code//analysis//gcpy_campaigns//')
import obspack_io as obs
import utils as ut

sys.path.insert(0, 'C:\\Users\\jhask\\OneDrive\\Documents\\Research\\Projects\\MIT\\pOrgNO3\\code\\analysis\\gcpy_campaigns\\')
import obspack_io as obs
pth= 'C://Users//jhask//OneDrive//Desktop//test//'

#Restore the data we loaded for fastness.
ds=xr.open_dataset('C://Users//jhask//OneDrive//Documents//Research//Projects//MIT//UROP//code//NADP_NTN_native_All.nc')
ds = ds.rename_dims({'index':'time'}) # because obs contains repeated values if both sampled at same time. 
ds.coords['time']=('time', ds.index) 
ds=ds.reset_coords(names='index', drop=True)
bb=ds.sel(time=ds.time[-1])

df=bb.to_dataframe()
df2= df[df[['longitude','latitude','elevation']].notnull().all(1)]
df3=df2[['time','longitude','latitude','elevation']].copy()
tms= pd.to_datetime(ds.time.values).to_series().reset_index(drop=True)


obs.write_obspack_inputs_flights(df3.index, df3.latitude , df3.longitude , df3.elevation,
                               tms, sample_stragety=2, 
                               outpath=pth) 



