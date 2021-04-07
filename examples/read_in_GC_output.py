# -*- coding: utf-8 -*-
"""
Script with examples of reading in ObsPack & Planelight.dat output files
from GEOS-Chem. 

@author: Dr. Jessica D. Haskins
"""

import xarray as xr
import pandas as pd
import matplotlib as plt 
import os 

import write_obspack_inputs as obs
import write_planeflight_inputs as pln

# Get the current working directory
cwd = os.getcwd()

# =============================================================================
# ############### Read in & Concatenate output ObsPack files  #################
# =============================================================================  

# Path to folder containing example GEOS-Chem ObsPack Output 
filepath = cwd+'/datafiles_for_examples/'

# Concatonate our ObsPack netcdf files. (only need to do once after you get output...)
fileconcat = obs.read_and_concat_output_files(filepath, outfile='all_ObsPacks')

# fileconcat=filepath+'/all_ObsPacks').nc' Use this to restore after concatenating! 
obsdat=xr.open_dataset(fileconcat) # open obspack file 

# Plot the outputted O3 from ObsPack
tracer='O3' 
fig, ax = plt.subplots(figsize=(15, 9))
ax.plot(obsdat.time, obsdat[tracer], label='ObsPack',linewidth=2) 
ax.legend()
plt.tight_layout()
        
plt.show()
            

# =============================================================================
# ############### Read in & Concatenate output plane.log files  #################
# =============================================================================  

#  OPTION 1: Read in a single plane.log files as a pandas dataframe, 
pdat= pln.read_planelog(filepath+'/plane.log.20130603') # open planelog file

############################################################################
#  OPTION 2: Read in several plane.log files as a pandas dataframe, 
# concatenate data from all flights into single dataframe, and save 
# all as a pickle! 
############################################################################


# Path to folder containing GEOS-Chem ObsPack Output 
filepath = cwd+'/datafiles_for_examples/'

# Concatonate our Plane.log files. (only need to do once after you get output...)
p= pln.planelog_read_and_concat(filepath) 

# Use this to restore saved file of concatenated planelogs! 
pdat = pd.read_pickle(filepath+'planelog_concat') 

# Plot the outputted O3 from plane.log
tracer='O3' 
fig, ax = plt.subplots(figsize=(15, 9))
ax.plot(pdat.time, pdat[tracer], label='Plane.log',linewidth=2) 
ax.legend()
plt.tight_layout()
        
plt.show()
