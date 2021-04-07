"""
Scripts to make readable Obspack input files for GEOS-Chem for planes or ground sites
and read output from GEOS-Chem, concatonate files into a single .nc file for easy use.

Created on Sun Feb 14 12:05:25 2021.

@author: Dr. Jessica D.Haskins
"""
import pandas as pd
import numpy as np
import xarray as xr
import os 
import datetime
import glob
import sys
import utils as ut
from shutil import copyfile

def write_obspack_inputs_flights(sitename: str, lat , lon , alt,
                               time, sample_stragety: int, outpath: str):
    """Create ObsPack input files for GEOS Chem v13.0.0. for a Flight.
    
    # =========================================================================
    #                             INPUTS
    # =========================================================================
    # sitename = string with sitename used in obspack ID, (e.g.'SOAS-Tower')
    # lat = np array with latitudes of plane
    # lon = np array with longitudes of plane
    # alt = np array with plane height in meters above sea level
    # time  =  datetime series (usually pandas) witht times you want samples at.
    #
    # Sample Stragety = Integer specifiing how the model will average to
    #                   the timebase given (calc'd by  datestart,dateend,
    #                   and freq). Only Valid option are:
    #                             1  =   4-hour avg;
    #                             2  =   1-hour avg;
    #                             3  =   90-min avg;
    #                             4  =   instantaneous
    #
    # outfile_path  = string containing absolute path where the netcdf files
    #                 will be written.
    #
    # =========================================================================
    #                             OUTPUTS
    # =========================================================================
    #
    #  out= Returns list of netcdf files written to outpath for each day in
    #       the time period that  can be used as GEOS Chem inputs.
    #
    #  NOTE: The output files MUST be named as they are here in order to be
    #        read in by GEOS-Chem properly. So do NOT change the names.
    #        The convention is as follows:
    #        filename = 'obspack_input.YYYYMMDD.nc'
    """
    # Create dataframe with Y M D H M S as columns with integer values, 
    # convert to 2D numpy. Contains entire time period. 
    df=pd.DataFrame() 
    df['year'] = time.dt.strftime("%Y").astype(int)
    df['mon'] =time.dt.strftime("%m").astype(int)
    df['day'] = time.dt.strftime("%d").astype(int)
    df['hr'] =time.dt.strftime("%H").astype(int)
    df['min'] = time.dt.strftime("%M").astype(int)
    df['s'] = time.dt.strftime("%S").astype(int)
    time_components = df.to_numpy()
    
    # Get the Date only, no HMS so we can select points from a specifc day
    date_only=time.dt.strftime("%Y%m%d") # List of all days, for full time.
    indv_days= np.unique(date_only) # List of Unique days for whole period. 
    
    # Check to make sure passed values are valid. 
    if (any(abs(lat)> 90)) or  (np.isfinite(lat).any()  ==False or (any(np.isnan(lat))==True)):
        sys.exit('Latitudes for '+sitename+ 'are either > 90 or contain non -Finite values')
    if (any(abs(lon)> 180)) or  (any(np.isfinite(lon))  ==False) or (any(np.isnan(lon))==True):
        sys.exit('Longitudes for '+sitename+' are either > 180 or contain non -Finite values')
    if (any(alt < 0)) or  (any(np.isfinite(alt))  ==False) or (any(np.isnan(alt))==True):
        sys.exit('Altitude for '+sitename+' are < 0 or contain non -Finite values')
    if (len(time_components[np.isnan(df)])> 0):
        print(time_components[np.isnan(time_components)])
        sys.exit('Times for '+sitename+' contain non -Finite values')
        
    # Create unique obspack array of ID strings: for entire time period:
    n = np.arange(0, len(date_only)).astype(str)  # Create a uniq obs #
    prefix = sitename + '_from_' + \
        date_only.min()+ '_to_' + date_only.max()+ '_n'
    ids = np.array([(prefix + stri) for stri in n]).astype('|S200') # Make ID. 
    
    # ==================   Begin making Indv Files:  ==========================
    all_files = list()  # empty list to contain netcdf filenames generated.
    
    # Loop over each day you need an ObsPack File for. Grab time, lat lon, alt, ids, 
    # only on this date. Create the xarray and then save the .nc file. 
    for i in range(0, len(indv_days)): 
        # Get index of obs points on this day from larger arrays. 
        inds= np.where(date_only == indv_days[i])[0] 
        
        these_times= time_components[inds, :] 
        these_lats= lat[inds].astype(np.float32)
        these_lons= lon[inds].astype(np.float32)
        these_alts= alt[inds].astype(np.float32)
        these_ids= ids[inds].astype('|S200')
    
        # ====================================================================
        # =========    Make Big XArray with everything required.      ========
        # ====================================================================
        # Format of vars & attributes comes from GC example shown here:
        #   http://wiki.seas.harvard.edu/geos-chem/index.php/ObsPack_diagnostic
        
        ds = xr.Dataset({
            'latitude': xr.DataArray(
                        data=these_lats,
                        dims=['obs'],
                        attrs={
                            "units": "degrees_north",
                            "_FillValue": -1.0e+34,
                            "long_name": "Sample latitude"
                        }),
            'longitude': xr.DataArray(
                data=these_lons,
                dims=['obs'],
                attrs={
                    "units": "degrees_east",
                    "_FillValue": -1.0e+34,
                    "long_name": "Sample longitude"
                }),
            'altitude': xr.DataArray(
                data=these_alts,
                dims=['obs'],
                attrs={
                    "units": "meters",
                    "_FillValue": -1.0e+34,
                    "long_name": "sample altitude in meters above sea level",
                    "comment": "Altitude is elevation plus sample intake height in meters above sea level."
                }),
            # Array of Y M D H M S to sample Must be UTC.
            'time_components': xr.DataArray(
                data=these_times,
                dims=['obs', 'calendar_components'],
                attrs={
                    "_FillValue": -9,
                    "long_name": " Calendar time components as integers. Times and dates are UTC.",
                    "order": "year, month, day, hour, minute, second",
                    "comment": "Calendar time components as integers.  Times and dates are UTC."
                }),
            # And now also pass our unique IDs (for each sample point).
            # Dimension is observation #, 200 len string.
            'obspack_id': xr.DataArray(
                data=these_ids,
                dims=['obs'],
                attrs={
                    "long_name": "Unique ObsPack observation id",
                    "comment": "Unique observation id string that includes obs_id, dataset_id and obspack_num."
                }),
            # Duplicate the int that indicates sampling stragety as many times
            # as you want to sample the model.
            'CT_sampling_strategy': xr.DataArray(
                data=np.full(len(inds), sample_stragety, dtype=int),
                dims=['obs'],
                attrs={
                    "_FillValue": -9,
                    "long_name": "Unique ObsPack observation id",
                    "values": "How to sample model. 1=4-hour avg; 2=1-hour avg; 3=90-min avg; 4=instantaneous"
                })
        })
        
        # ====================================================================
        # NOTE: FILES MUST BE NAMED LIKE THIS IN ORDER TO BE READ IN GEOS CHEM
        #       OBSPACK_MOD.F90
        filename = 'obspack_input.'+ str(indv_days[i])+ '.nc'
        # ====================================================================
    
        # Convert our Xarray Data set to a netcdf file w/ this name at outpath
        ds.to_netcdf(outpath + filename)
        
        # Append name of this file to our list of file names that we return
        # to the user.
        all_files.append(outpath + filename)
        
    return all_files
    

def write_obspack_inputs_ground(sitename: str, lat: int, lon: int, alt: int,
                               datestart: str, dateend: str, samplefreq: int,
                               sample_stragety: int, outpath: str):
    """Create ObsPack input files for GEOS Chem v13.0.0. for a single gound site.
    
    # =========================================================================
    #                             INPUTS
    # =========================================================================
    # sitename = string with sitename used in obspack ID, (e.g.'SOAS-Tower')
    # lat = integer with latitude of obs site in degrees north (e.g. 35.45)
    # lon = integer with longitude of obs site in degrees east
    # alt = integer with elevation plus sample intake height in meters above
    #       sea level
    # datestart = string format YYYYMMDD HH:MM:SS' indicating when you need
    #              obsPack netcdf  file to begin
    # dateend = string format YYYYMMDD HH:MM:SS' indicating when you need
    #              obsPack netcdf  file to stop
    # freq = integer for # of seconds you'd like to step betweeen datestart,
    #         dateend (e.g. 3600 for hourly steps)
    # Sample Stragety = Integer specifiing how the model will average to
    #                   the timebase given (calc'd by  datestart,dateend,
    #                   and freq). Only Valid option are:
    #                             1  =   4-hour avg;
    #                             2  =   1-hour avg;
    #                             3  =   90-min avg;
    #                             4  =   instantaneous
    #
    # outfile_path  = string containing absolute path where the netcdf files
    #                 will be written.
    #
    # =========================================================================
    #                             OUTPUTS
    # =========================================================================
    #
    #  out= Returns list of netcdf files written to outpath for each day in
    #       the time period that  can be used as GEOS Chem inputs.
    #
    #  NOTE: The output files MUST be named as they are here in order to be
    #        read in by GEOS-Chem properly. So do NOT change the names.
    #        The convention is as follows:
    #        filename = 'obspack_input.YYYYMMDD.nc'
    #
    # =========================================================================
    #                              Example
    # =========================================================================
    #
    # SOAS Centerville Ground collection site was at lat=32.903281,
    # lon=-87.249942, alt=125. Let's sample the model there ever hour
    # (outputting the hourly average).This code snippiet will make ALL the
    # .nc files we need for ObsPack input for the length of our run (6/1/2013-
    # 7/15/2013) and save them to my desktop.
    #
    # import obspack_fill as obs # place function in same path as call.
    # import xarray
    #
    # filename = obs.make_ObsPack_Input_netcdfs('SOAS-Ground', lat=32.903281,
    #                   lon=-87.249942, alt=125,datestart='20130601 00:00:00',
    #                   dateend='20130716 00:00:00', samplefreq=3600,
    #                   sample_stragety=2, outpath='C:/Users/jhask/Desktop/')
    #
    #  # Open the first file with xarray and print info about it.
    #  dat = xr.open_dataset(filename[0])
    #  print(dat)
    #
    """
    # Check to make sure passed values are valid. 
    if (abs(lat)> 90) or  (np.isfinite(lat) is False):
        sys.exit('Latitudes for '+sitename+'are either > 90 or contain non -Finite values')
    if (abs(lon)> 180) or  (np.isfinite(lon)  is False):
        sys.exit('Longitudes for '+sitename+' are either > 180 or contain non -Finite values')
    if (alt < 0) or  (np.isfinite(alt)  is False):
        sys.exit('Altitudes for '+sitename+' are < 0 or contain non -Finite values')
        
    # Figure out # of days between start & Enddate (e.g. # of NC files we need)
    start = datetime.datetime.strptime(datestart, '%Y%m%d %H:%M:%S')
    endd = datetime.datetime.strptime(dateend, '%Y%m%d %H:%M:%S')
    delta = endd - start

    if (delta.seconds != 0):
        days = delta.days + 1  # add a day if there's some leftover time
    else:
        days = delta.days  # otherwise stick with the round number.
    starts = pd.date_range(datestart, dateend, freq='1D').to_series()
    starts2Use = starts[:days]  # only use the # of days we decided on.
    ends2Use = starts2Use + \
        datetime.timedelta(days=1) - datetime.timedelta(seconds=1)
    # don't do the whole day if the stop time is shorter.
    if ends2Use[-1] > endd:
        ends2Use[-1] = endd
    
    # ===================================================================
    all_files = list()  # empty list to contain netcdf filenames generated.
    
    # Loop over # of days we need a file for and make a netcdf!
    for t in range(0, len(starts2Use)):
        # Create the time components variable.
        dfi = pd.DataFrame()  # make empty data frame.
        dts = pd.date_range(  # Get datetimes in our day at sample freq
            starts2Use[t], ends2Use[t], freq=str(samplefreq) + 's')
        
        # Fill columns of dataframe with datetime objs
        dfi['datetime'] = dts
        
        # Strip YMD HMS out of the datetime as integers.
        dfi['year'] = dfi['datetime'].dt.strftime("%Y").astype(int)
        dfi['mon'] = dfi['datetime'].dt.strftime("%m").astype(int)
        dfi['day'] = dfi['datetime'].dt.strftime("%d").astype(int)
        dfi['hr'] = dfi['datetime'].dt.strftime("%H").astype(int)
        dfi['min'] = dfi['datetime'].dt.strftime("%M").astype(int)
        dfi['s'] = dfi['datetime'].dt.strftime("%S").astype(int)
        dfi = dfi.drop(columns='datetime')  # drop datetime column
        
        time_components = dfi.to_numpy()  # convert to a numpy array b4 xarray
        
        # =====================================================================
        # Create unique obspack array of ID strings:
        n = np.arange(0, len(dts)).astype(str)  # Create a uniq obs #
        # Create a prefix containing the sitename, start/stop dates
        prefix = sitename + '_from_' + \
            datestart.split(' ')[0] + '_to_' + dateend.split(' ')[0] + '_n'
            
        # Make obspack IDs. 
        ids = np.array([(prefix + stri) for stri in n]).astype('|S200')
        
        # ====================================================================
        # =========    Make Big XArray with everything required.      ========
        # ====================================================================
        
        # Format of vars & attributes comes from GC example shown here:
        #   http://wiki.seas.harvard.edu/geos-chem/index.php/ObsPack_diagnostic
        
        ds = xr.Dataset({
            # Duplicate lat as many times as you collect obs points. Station
            # isnt' moving.
            'latitude': xr.DataArray(
                        data=np.full(len(dts), lat, dtype=np.float32),
                        dims=['obs'],
                        attrs={
                            "units": "degrees_north",
                            "_FillValue": -1.0e+34,
                            "long_name": "Sample latitude"
                        }),
            # Duplicate lon as many times as you collect obs points.
            'longitude': xr.DataArray(
                data=np.full(len(dts), lon, dtype=np.float32),
                dims=['obs'],
                attrs={
                    "units": "degrees_east",
                    "_FillValue": -1.0e+34,
                    "long_name": "Sample longitude"
                }),
            # Duplicate alt as many times as you collect obs points.
            'altitude': xr.DataArray(
                data=np.full(len(dts), alt, dtype=np.float32),
                dims=['obs'],
                attrs={
                    "units": "meters",
                    "_FillValue": -1.0e+34,
                    "long_name": "sample altitude in meters above sea level",
                    "comment": "Altitude is elevation plus sample intake height in meters above sea level."
                }),
            # This is the only thing changing- array of Y M D H M S to sample
            # this lat,on, alt point. Must be UTC.
            'time_components': xr.DataArray(
                data=time_components,
                dims=['obs', 'calendar_components'],
                attrs={
                    "_FillValue": -9,
                    "long_name": " Calendar time components as integers. Times and dates are UTC.",
                    "order": "year, month, day, hour, minute, second",
                    "comment": "Calendar time components as integers.  Times and dates are UTC."
                }),
            # And now also pass our unique IDs (for each sample point).
            # Dimension is obs #, 200 len string.
            'obspack_id': xr.DataArray(
                data=ids,
                dims=['obs'],
                attrs={
                    "long_name": "Unique ObsPack observation id",
                    "comment": "Unique observation id string that includes obs_id, dataset_id and obspack_num."
                }),
            # Duplicate the int that indicates sampling stragety as many times
            # as you want to sample the model.
            'CT_sampling_strategy': xr.DataArray(
                data=np.full(len(dts), sample_stragety, dtype=int),
                dims=['obs'],
                attrs={
                    "_FillValue": -9,
                    "long_name": "Unique ObsPack observation id",
                    "values": "How to sample model. 1=4-hour avg; 2=1-hour avg; 3=90-min avg; 4=instantaneous"
                })
        })
        # ====================================================================
        # NOTE: FILES MUST BE NAMED LIKE THIS IN ORDER TO BE READ IN GEOS CHEM!
        filename = 'obspack_input.'+ str(starts2Use[t]).split(' ')[0].replace('-', '') + '.nc'
        # ====================================================================
        
        # Convert our Xarray Data set to a netcdf file w/ this name at outpath
        ds.to_netcdf(outpath + filename)
        
        # Append name of this file to our list of file names that we return
        # to the user.
        all_files.append(outpath + filename)
    
    return all_files


def combine_common_obspack_inputs(folder_1: str , folder_2:str , outpath: str, 
                            copy_not_common : bool = False):
    """Look in two folders for ObsPack Files and combine common date files.
    
    This function is usefule if you want to sample the model on the same date 
    at both a ground site and flight based path. For example, maybe  I want to 
    sample the GEOS-Chem at the SOAS Centerville site and along the SENEX 
    flight path during June- July of 2013. Some flights took place on the 
    same days as continous ground monitoring. I used the  
    write_obsPack_inputs_ground() function to make ObsPack files for SOAS 
    and I used the write_obspack_inputs_flights() function for SENEX, and then
    this function to combine the files for the same dates (in time order). 
    
    You could seperate the output by either using the lat/lon of the ground site 
    or probably, easier, by filtering by the obspack_ID string (which is unique
    for those from SOAS vs those from SENEX.)
    
    Args: 
    ----
        folder_1: string path to a folder containing ObsPack files. 
                    (outpath for make Obspack Inputs)
        folder_2: string path to a folder containing ObsPack files. 
                    (outpath for make Obspack Inputs)
        outpath: string to path you want to store the combined ObsPack files in.
        
        copy_not_common: Boolean if you want to copy files that the two folders
            don't have in common to this outpath. So have a folder of all 
            ObsPack files from Folder 1, folder 2, withe the commons combined 
            and the not commons jsut copied over. Default is False. 
            
    Output: 
    ------
        files_combo= list of files that were combined and saved at outpath.
    """
    # Get all file names in these source folders: 
    filenames_1=glob.glob(folder_1+"*.nc")
    filenames_2=glob.glob(folder_2+"*.nc")
    
    # Split the path of the filename, and compare just the dates of the files: 
    dates1= [filenames_1[i].split('obspack_input.')[1] for i in range(0,len(filenames_1))]
    dates2= [filenames_2[i].split('obspack_input.')[1] for i in range(0,len(filenames_2))]
    
    # Get list of Common dates we have 2 ObsPack files for: 
    common = list(set(dates1) & set(dates2))
    common = ['obspack_input.' + sub for sub in common] 
    
    # Loop over each file we have 2 ObsPack files for: 
    for i in range(0, len(common)):  
        # Full path of each file for the same date: 
        indv_filename1= folder_1+ common[i]
        indv_filename2= folder_2+ common[i]
        
        # Open each data set: 
        one = xr.open_dataset(indv_filename1)
        two = xr.open_dataset(indv_filename2)
        
        # Combine the two data sets: 
        both= xr.concat([one, two], dim='obs')
        
        # Take time componenets, make datetimes so we can sort "both" in time.
        for j in range(0, len(both.latitude)):
            tms= both.time_components.values.astype(int)
            out =datetime.datetime(tms[j,0],tms[j,1],
                                    tms[j,2],tms[j,3],tms[j,4],tms[j,5])
            if j ==0 : 
                timeline=np.array(out)
            else: 
                timeline=np.append(timeline, out)
        ind= np.argsort(timeline) # indexes of "both" sorted in time.
        
        # Re-Index "Both" so that everything is in time order. 
        both.latitude.values = both.latitude.values[ind]
        both.longitude.values = both.longitude.values[ind]
        both.altitude.values = both.altitude.values[ind]
        both.obspack_id.values = both.obspack_id.values[ind]
        both.CT_sampling_strategy.values = both.CT_sampling_strategy.values[ind]
        both.time_components.values = both.time_components.values[ind, :]

        # Convert combined DataSet to a netcdf file w/ this name at outpath
        both.to_netcdf(outpath+ common[i])
        
    print('Common files saved at: ', outpath) # Let 'em know.
    
    if copy_not_common is True: 
        # Get list of dates we don't have 2 ObsPack files for: 
        not_common = list(set(dates1) ^ set(dates2))
        not_common = ['obspack_input.' + sub for sub in not_common] 
        for i in range(0, len(not_common)):  
            # Two possible paths to this file in each folder: 
            indv_filename1= folder_1+ not_common[i]
            indv_filename2= folder_2+ not_common[i]
            
            # Copy it over to the new outpath from its source when you find it: 
            if indv_filename1 in filenames_1: 
                copyfile(indv_filename1, outpath+not_common[i])
                
            if indv_filename2 in filenames_2: 
                copyfile(indv_filename2, outpath+not_common[i])
                
        print('Uncommon files copied over to: ', outpath)
                
    return common


def _switch_obs_2_time_dim(ds):
    """Function to create a single time variable that is the midpoint of the 
    ObsPack averaging interval, and make it the xarray coordinate. """ 
    # Get the midpoint of the average pulled from the model:
    midpoint = pd.to_datetime(ds.averaging_interval_start.data) + \
        np.asarray(ds.averaging_interval.data) / 2
    
    # Make it the time midpoint a new variable in the dataset.
    t = midpoint.to_series().reset_index(drop=True)
    ds['time'] = ("obs", t)
    
    # Tell xarray that we want time to be a coordinate.
    ds = ds.set_coords('time')
    
    # And tell it to replace Obs # with time as the preferred dimension.
    ds = ds.swap_dims({"obs": "time"})
    
    return ds

    
def read_and_concat_output_files(path_to_dir: str, outdir: str = None,
                   outfile: str = None):
    """
    Concatonate output ObsPack files into a single file from a directory.
    # If you only want to open a single file, try xarray.open_dataset().
    
    Args
    ----
        path_to_dir = String with filepath to folder containing GEOS Chem output 
                      Obspack netcdf files
    
        outdir(optional) = # String path to where concatonated file will be saved 
    
        outfilename (optional) = string name of output file, no extension needed.
    """
    # If outdir not set, then set it to the same as the input file path.
    outdir = path_to_dir if outdir is None else outdir
    
    # If outfilename not set, then set it to be concat_ObsPack.nc
    outfile = 'concat_ObsPack' if outfile is None else outfile
    
    # Get a list of variables that GCPy should not read.
    # These are mostly variables introduced into GCHP with the MAPL v1.0.0
    # update.  These variables contain either repeated or non-standard
    # dimensions that can cause problems in xarray when combining datasets.
    skip_vars = ["anchor","ncontact","orientation","contacts","cubed_sphere"]
    
    # Look for all the ObsPack netCDF files in the path
    file_list = ut._find_files_in_dir(path_to_dir, ['ObsPack', '.nc4'])
    
    # Return a single xarray Dataset containing data from all files
    # Specify a pre-processing function that makes time midpoint the coordinate.
    ds = xr.open_mfdataset(file_list, drop_variables=skip_vars,
                                 preprocess=_switch_obs_2_time_dim)

    # Specify the path and filename for the concatenated data
    out = os.path.join(outdir, outfile+'.nc')
    
    # Write concatenated data to a netCDF file
    ds.to_netcdf(out)
    
    print('File saved to:' + out)
    
    return out 
