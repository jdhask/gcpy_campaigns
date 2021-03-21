"""
Concatenates several individual ObsPack ouput netCDF filesinto a single netCDF file.

Created on Sun Feb 14 18:16:25 2021.

@author: Dr. Jessica D. Haskins
"""

import os
import numpy as np
import xarray as xr
import pandas as pd

def _find_files_in_dir(path, substrs):
    """
    Return a list of all files in a directory that match substrings.

    Args
    ----
        path : str
            Path to the directory in which to search for files.

        substrs : list of str
            List of substrings used in the search for files.

    Returns
    -------
        file_list : list of str
            List of files in the directory (specified by path)
            that match all substrings (specified in substrs).
    """
    # Initialize
    file_list = []

    # Walk through the given data directory.  Then for each file found,
    # add it to file_list if it matches text in search_list.
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))

    # Return an alphabetically sorted list of files
    file_list.sort()
    return file_list


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


def obspack_concat_netcdfs(path_to_dir: str, outdir: str = None,
                   outfile: str = None):
    """
    Concatonate output ObsPack NetCDF files into a single file with 
    averaging interval midpoint as a coordinate.
    
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
    
    # Look for all the netCDF files in the path
    file_list = _find_files_in_dir(path_to_dir, 'ObsPack')
    
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
