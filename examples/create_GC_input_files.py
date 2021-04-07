"""
Script to make ObsPack & Planelight.dat input files for individual campaigns.

Created on Sun Feb 14 19:25:16 2021

@author: Dr. Jessica D. Haskins
"""

import xarray as xr
import pandas as pd
import numpy as np 
import os 

import write_obspack_inputs as obs
import write_planeflight_inputs as pln

# Get the current working directory
cwd = os.getcwd()

# =============================================================================
# #########################  SENEX planeflight.date files  ####################
# =============================================================================

# Open SENEX merged data. Convert times to a pandas datetime array. 
senex=cwd+"/examples/datafiles_for_examples/SENEX.nc"
sn=xr.open_dataset(senex) 
dt=pd.to_datetime(sn.time.values).to_series().reset_index(drop=True)

#  OPTION 1:  Make a planeflight.dat file for all SENEX flights. Sample just NO, O3, CO. 
os.mkdir(cwd+'\\example1\\') # Make directory to hold output. 
pln.make_planeflightdat_files(outpath=cwd+'\\example1\\',
                            datetimes=dt,
                            lat_arr=sn.GpsLat.values, 
                            lon_arr=sn.GpsLon.values, 
                            pres_arr=sn.StaticPrs.values,
                            typestr='aSENEX', # Don't make typestrs begin with "S" or errors! 
                            tracers=['NO', 'O3', 'CO'],
                            username= 'jdh',
                            drop_dupes = False,
                            diags = ['all'])

#  OPTION 2:  Make a planeflight.dat file for all SENEX flights. 
# Sample all tracers in your input file and all optional diagnostics
os.mkdir(cwd+'\\example2\\') # Make directory to hold output. 
pln.make_planeflightdat_files(outpath=cwd+'\\example2\\',
                            input_file=cwd+'\\examples\\datafiles_for_examples\\input.geos', 
                            typestr='aSENEX', # Don't make typestrs begin with "S" or errors! 
                            datetimes=dt,
                            lat_arr=sn.GpsLat.values, 
                            lon_arr=sn.GpsLon.values, 
                            pres_arr=sn.StaticPrs.values,
                            username= 'jdh',
                            drop_dupes = False,
                            diags = ['all'],
                            print_diag_options= False)
 
#  OPTION 3:  Make a planeflight.dat file for all SENEX flights. 
# Sample all tracers in your input file except those in "minus". 
os.mkdir(cwd+'\\example3\\') # Make directory to hold output. 

minus=["AODC_SULF","AODC_BLKC","AODC_ORGC","AODC_SALA","AODC_SALC",
"AODC_DUST","AODB_SULF","AODB_BLKC","AODB_ORGC","AODB_SALA","AODB_SALC",
"AODB_DUST","HG2_FRACG" ,"HG2_FRACP","GMAO_ICE00","GMAO_ICE10","GMAO_ICE20", 
"GMAO_ICE30", "GMAO_ICE40","GMAO_ICE50", "GMAO_ICE60","GMAO_ICE70", "GMAO_ICE80","GMAO_ICE90"]

pln.make_planeflightdat_files(outpath=cwd+'\\example3\\',
                            input_file=cwd+'\\examples\\datafiles_for_examples\\input.geos', 
                            campaign='aSENEX',
                            datetimes=dt,
                            lat_arr=sn.GpsLat.values, 
                            lon_arr=sn.GpsLon.values, 
                            pres_arr=sn.StaticPrs.values,
                            username= 'jdh',
                            drop_dupes = False,
                            diags = ['all'],
                            diags_minus=minus,
                            print_diag_options= False)


# =============================================================================
#    ################  SOAS Centreville  ObsPack  Ground Sites ################
# =============================================================================
# Create ObsPack files for the Centreville, Alabama site during SOAS 
# sampling the model there every hour of the campaign between 6/1/2013 and 7/15/2013.

os.mkdir(cwd+'\\SOAS\\') # Make directory to hold output. 

filenames = obs.write_obsPack_inputs_ground('SOAS-Ground', lat=32.903281, 
                                      lon=-87.249942, alt=2,
                                      datestart='20130601 00:00:00',
                                      dateend='20130715 00:00:00', 
                                      samplefreq=3600,
                                      sample_stragety=2, 
                                      outpath=cwd+'\\SOAS\\')



# =============================================================================
#   ####################    SENEX Flights ObsPack     #####################
# =============================================================================

# Open SENEX data to get lat, long, alt info. 
# Open SENEX merged data. Convert times to a pandas datetime array. 
senex=cwd+"\\examples\\datafiles_for_examples\\SENEX.nc"
sn=xr.open_dataset(senex) 

# Make sure my alt, lat, lon, time data is formatted correctly: 
alt = sn.GpsAlt.values # Pull alt from ds
alt[alt < 0] = np.NaN  # Don't allow negatives. 
# Create smaller xarray with ONLY the data we need, drop any NaNs from it! 
sen =  xr.Dataset({'alt':xr.DataArray( data= alt , dims=['obs']), 
                      'lat':xr.DataArray( data= sn.GpsLat.values , dims=['obs']),
                      'lon':xr.DataArray( data= sn.GpsLon.values , dims=['obs']),
                      'time':xr.DataArray( data= sn.time.values , dims=['obs'])}).dropna(dim='obs')
# Get time to not be an index and not be weird. 
time= pd.to_datetime(sen.time.values).to_series().reset_index(drop=True)

os.mkdir(cwd+'\\SENEX\\') # Make directory to hold output. 
# Now make ObsPack files for all flights: 
filenames = obs.write_obspack_inputs_flights('SENEX',lat=sen.lat.values, 
                                                  lon=sen.lon.values,
                                                  alt=sen.alt.values,
                                                  time=time, 
                                                  sample_stragety=4,
                                                  outpath=cwd+'\\SENEX\\')

# =============================================================================
#############    Combine SENEX and SOAS files on common dates   ###############
# =============================================================================
# I want to sample GEOS-Chem at the SOAS Centerville site and 
# along the SENEX flight path during June- July of 2013. Some flights took
# place on the same days as continous ground monitoring, so we have mutliple files for 
# the same day which isn't allowed as input, so we need to combine them.
# (( You could seperate the output by either using the lat/lon of the ground site 
# or probably, easier, by filtering by the obspack_ID string (which is unique
# for those from SOAS vs those from SENEX.)))

os.mkdir(cwd+'\\both\\') # Make directory to hold output. 

#Folders to each of their indvidual ObsPack Files 
senex_folder = cwd+'\\SENEX\\'
soas_folder = cwd+'\\SOAS\\'
pathy= cwd+'\\both\\' # output folder for all files.

commons= obs.combine_common_ObsPacks(soas_folder, senex_folder, outpath=pathy, 
                            copy_not_common= True)

# Open a file and take a peek for sanity! 
t=xr.open_dataset(pathy+commons[0])

lat= t.latitude.values
lon= t.longitude.values
alt= t.altitude.values
tm= t.time_components.values
obspack_id= str(t.obspack_id.values)
sample= t.CT_sampling_strategy.values
