"""
Script to make readable Planelight.dat.YYYYMMDD input files for GEOS-Chem.

Created on Sun Mar 21 14:21:14 2021

@author: Dr. Jessica D. Haskins
"""

import pandas as pd 
from datetime import datetime
import numpy as np
import sys
import os

def _build_species_list_from_input(input_filename: str):
    """Parse an input.geos file and build a list of advected species."""
    inputf = open(input_filename, 'r')  # Open the input file.
    
    tracers = []  # Make an empty list that will contain advected species.
    adv_species = False  # Line 1 doesn't contain Adv. Sepcies.
    count = 0  # Initialize line counting variable.

    while True:  # while not at the end of the file...
        count += 1  # update line counter
        line = inputf.readline()  # Read next line from file
        
        if 'END OF FILE' in line: break
        
        # Next lines will NOT contain species names if in Transport Menu.
        adv_species = False if 'TRANSPORT MENU' in line else adv_species
            
        if adv_species is True:  # Append speices name to list of tracers!
            if ":" in line:  # don't parse lines that aren't "assignments"
                tracers=tracers+[(line.strip().split(':')[1]).strip()]
        
        # Next lines will contain species names, tell it to expect them!
        adv_species = True if 'ADVECTED SPECIES MENU' in line else adv_species
        
    return tracers

def _get_optional_diags(print_diag_options: bool = False,) :
    # List of all optional Diagnostic Quantites you can get with Planeflight 
    # grabbed from: http://wiki.seas.harvard.edu/geos-chem/index.php/Planeflight_diagnostic
    # and from inside planeflight_mod.f90
    diag_dict= dict({"RO2": "Concentration of RO2 family", 
                     "AN": "Concentration of AN family", 
                     "NOy": "Concentration of NOy family",
                     "GMAO_TEMP":"Temperature",
                     "GMAO_ABSH":"Absolute humidity",
                     "GMAO_SURF":"Aerosol surface area",
                     "GMAO_PSFC":"Surface pressure", 
                     "GMAO_UWND":"Zonal winds",
                     "GMAO_VWND":"Meridional winds", 
                     "GMAO_IIEV":"GEOS-Chem grid box index, longitude", 
                     "GMAO_JJEV":"GEOS-Chem grid box index, latitude",
                     "GMAO_LLEV":"GEOS-Chem grid box index, altitude", 
                     "GMAO_RELH":"Relative humidity", 
                     #"GMAO_PVRT":"Ertel's potential vorticity", # Currently disabled, State_Met%PV is not defined.
                     "GMAO_PSLV":"Sea level pressure", 
                     "GMAO_AVGW":"Water vapor mixing ratio", 
                     "GMAO_THTA":"Potential temperature",
                     "GMAO_PRES":"Pressure at center of grid box", 
                     "AODC_SULF":"Column aerosol optical depth for sulfate", 
                     "AODC_BLKC":"Column aerosol optical depth for black carbon",
                     "AODC_ORGC":"Column aerosol optical depth for organic carbon",
                     "AODC_SALA":"Column aerosol optical depth for accumulation mode sea salt",
                     "AODC_SALC":"Column aerosol optical depth for coarse mode sea salt", 
                     "AODC_DUST":"Column aerosol optical depth for dust", 
                     "AODB_SULF":"Column aerosol optical depth for sulfate below the aircraft",
                     "AODB_BLKC":"Column aerosol optical depth for black carbon below the aircraft",
                     "AODB_ORGC":"Column aerosol optical depth for organic carbon below the aircraft",
                     "AODB_SALA":"Column aerosol optical depth for accumulation mode sea salt below the aircraft",
                     "AODB_SALC":"Column aerosol optical depth for coarse mode sea salt below the aircraft",
                     "AODB_DUST":"Column aerosol optical depth for dust below the aircraft",
                     "HG2_FRACG":"Fraction of Hg(II) in the gas phase", "HG2_FRACP":"Fraction of Hg(II) in the particle phase",
                     "ISOR_HPLUS":"ISORROPIA H+", 
                     "ISOR_PH":"ISORROPIA pH (non-ideal system, so pH can be negative)",
                     "ISOR_AH2O":"ISORROPIA aerosol water", 
                     "ISOR_HSO4":"ISORROPIA bifulfate",
                     "TIME_LT":"Local time", 
                     "GMAO_ICE00": "the fraction of each grid box that has 0% to +10% of sea ice coverage",
                     "GMAO_ICE10" :"the fraction of each grid box that has 10% to +20% of sea ice coverage",
                     "GMAO_ICE20" : "the fraction of each grid box that has 20% to +30% of sea ice coverage",
                     "GMAO_ICE30" : "the fraction of each grid box that has 30% to +40% of sea ice coverage",
                     "GMAO_ICE40" : "the fraction of each grid box that has 40% to +50% of sea ice coverage",
                     "GMAO_ICE50" : "the fraction of each grid box that has 50% to +60% of sea ice coverage",
                     "GMAO_ICE60": "the fraction of each grid box that has 60% to +70% of sea ice coverage",
                     "GMAO_ICE70" : "the fraction of each grid box that has 70% to +80% of sea ice coverage",
                     "GMAO_ICE80": "the fraction of each grid box that has 80% to +90% of sea ice coverage",
                     "GMAO_ICE90": "the fraction of each grid box that has 90% to +100% of sea ice coverage",
                     #"GAMM_EPOX":"Uptake coefficient for EPOX", # don't currently work with flexchem
                     #"GAMM_IMAE":"Uptake coefficient for IMAE", 
                     #"GAMM_ISOPN":"Uptake coefficient for ISOPN", 
                     #"GAMM_DHDN":"Uptake coefficient for DHDN", 
                     #"GAMM_GLYX":"Uptake coefficient for GLYX", 
                     "AQAER_RAD":"Aqueous aerosol radius", 
                     "AQAER_SURF":"Aqueous aerosol surface area"})
    
    # Print off the optional diagnostics... 
    if print_diag_options is True: 
        k=[*diag_dict]
        print('Diagnositic options for planeflight.dat are as follows:')
        for i in range(0, len(k)):
            v=diag_dict.get(k[i])
            print("{:<15} {:<100}".format(k[i], v))
            
    return diag_dict
            
def make_planeflightdat_files(outpath: str, 
                           datetimes,
                           lat_arr, 
                           lon_arr, 
                           pres_arr= [],
                           alt_arr = [], 
                           tracers: list = [],
                           input_file: str ='',
                           campaign: str= '',
                           username: str = 'user',
                           overwrite: bool = True,
                           drop_dupes: bool = False,
                           diags: list = ['all'],
                           diags_minus: list =[],
                           print_diag_options: bool = False):
    """Function to create planeflight.dat files in correct format for GEOS-Chem.
    
    Args:
    ----
    outpath: string of path to save the output planeflight.dat files at. 
    datetimes: Pandas series of datetimes where obs to be collected at. 
    lat_arr: Array of latitudes where observations take at (degrees N)
    lon_arr: Array of longitudes where observations take at  (degrees N)
    
# Planeflight can take either alt or pressure. Must specify one. 
    pres_arr: Array of pressures in hPa where observations take at 
    alt_arr: Array of altitudes where observations take at  (meters)
    
# You can tell this function what tracers to sample with planeflight
# by either passing them as a list or by passing your input file. 
# That option will use all advected species.Not passing either arg will
# make a file without any tracers.
    input_file: String of path to your GEOS-Chem input file (to read in tracer names).
    tracers: List of tracers you want to sample using planeflight.dat
    
    Optional Arguements:
    -------------------
    campaign:   String of the "campaign" or the "type" of aircraft obs collected on. This is printed in the file.
    username:   String of user who created files. Optional. Arg in header of file.
    overwrite:  Boolean of whether to overwrite existing files at outpath with this name or not. 
    drop_dupes: Boolean of whether to drop duplicate rows 
    diags:      List of additional diagnostics to include (beyond tracers). Default option is 
               to include ALL available additional diagnositcs. 
    diags_minus: List of diagnostics you don't want to include (e.g. if you use "all").
    print_diag_options: Boolean of whether to print other available diagnostic options.
    
    Output:
    ------
    Files named planeflight.dat.YYYYMMDD written to outpath. 
    """
    if (len(pres_arr)==0) and (len(alt_arr)==0): 
        sys.exit('Please specify either an altitude or pressure array.')
    if (len(pres_arr)>0) and (len(alt_arr)>0): 
        sys.exit('Please specify either an altitude or pressure array, not both.'+\
                 'Planeflight.dat converts altitutdes to pressures. ')
    if len(campaign) > 7: 
        sys.exit('Please change the campaign string used as a type to be'+\
              'less than 7 chars, which is the max allowed by GEOS-Chem.')
            
    if len(diags)>0: # If the has user asked to include specific diagnostics... 
        # Build the dictionary of optional diagnostics. Print if user asks. 
        diag_dict= _get_optional_diags(print_diag_options)
        
        #If the user says to use all diagnostic quantities, then grab them from the dict. 
        if diags[0].lower() =='all': 
             optionals=[*diag_dict]
        else: # otherwise just use their inputted list! 
            optionals=diags
        
        if len(diags_minus): 
            optionals = [j for j in optionals if j not in diags_minus]
        
    # Make list of all tracers from input file plus the optional diagnostic quantites you want
    if input_file != '': #either by reading the input file 
        tracer_list= optionals+ _build_species_list_from_input(input_file)
    else: # or by using the tracesr they gave you. 
        tracer_list= optionals+ tracers
    ntracers=str(int(len(tracer_list))) # count the number of quantities. 
    
    # Designate a  few variables that we'll use to make the header lines of the files
    today= str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) # Time stamp when you made this file.
    spacer='-------------------------------------------------------------------------------'
    title = '  Now give the times and locations of the flight' 

    # Parse the passed dates, get list of indivudal dates we have obs for: 
    all_dates= datetimes.dt.date
    unq_dates=np.unique( datetimes.dt.date)

    # Loop over each unique day during campaign to make a planeflight.dat file for: 
    for i in range(0, len(unq_dates)): 
        # Create a filesname based on the date: 
        YYMMDD=unq_dates[i].strftime('%Y%m%d')
        filename='Planeflight.dat.'+ YYMMDD 
        
        #==========================================================================
        # Build arrays of all the data we need in Planeflight.dat for this date.
        # Rount to decmial places allowed in GEOS-Chem input for these vars. 
        #==========================================================================
        # Get indexes in larger array of flight obs that took place on this date.
        inds = np.where(all_dates== unq_dates[i])[0] 
        
        points= np.arange(1, len(inds)+1).astype(str)# Values for column "Points"
        obs= np.full(len(inds),9999.000) ## Observation value from the flight campaign
        typez= np.full(len(inds),'{typ: >6}'.format(typ=campaign)) #type  is allowed 7 chars. 
        day=   datetimes.dt.strftime('%d-%m-%Y')[inds]  # Day, month, and year (GMT date) for each flight track point
        tms=   datetimes.dt.strftime('%H:%M')[inds]  # Hour and minute (GMT time) for each flight track point
        lats=np.around( lat_arr[inds], decimals=2) # Latitude (-90 to 90 degrees) for each flight track point
        lons= np.around( lon_arr[inds], decimals=2)  # Longitude (-180 to 180 degrees) for each flight track point
        
        # Decide what array to use as vertical coordinate 
        if len(pres_arr)> 0: # Pressure in hPa for each flight track point.
            if i==1: print('Using PRESSURE, not altitude.')
            pres= np.around( pres_arr[inds], decimals=2) #GC doesn't allow more than 2 decimals
            vert_header='PRESS '  # String for header
            vert_arr=pres # Set the vertical array to pressure
        else:
            if len(alt_arr)> 0:  # Altitude in m for each flight track point.
                if i==1: print('Using ALTITUDE, not pressure.')
                alt = np.around(alt_arr[inds], decimals=2) #GC doesn't allow more than 2 decimals
                vert_header='ALT '  # String for header
                vert_arr=alt  # Set the vertical array to altitude
            
        # Pack all this into a dataframe, because pandas writes tab delimited files nicely! 
        df= pd.DataFrame({'POINT': points, 'TYPE': typez, 'DD-MM-YYYY': day, 
                              'HH:MM':tms, 'LAT':lats, 'LON':lons, vert_header.strip(): vert_arr, 
                              'OBS':obs})
        # Append a line at the bottom of the DF that says its the end! 
        df= df.append({'POINT': 99999, 'TYPE':'{typ: >6}'.format(typ='END'), 'DD-MM-YYYY': '00-00-0000', 
                              'HH:MM': '00:00', 'LAT':0.00, 'LON':0.00, vert_header.strip(): 0.00, 
                              'OBS':0.000}, ignore_index=True)
        
        #==========================================================================
        # Preform Checks on Data used in this file. Check for out of bounds, dupes, Nans.
        #==========================================================================
        # Check to make sure passed values are valid/ in range of correct units. 
        if any(abs(lats)> 90):
            sys.exit('| Latitudes|  are > 90.')
        if any(abs(lons)> 180):
            sys.exit('|Longitudes| are > 180 ')
        if len(alt_arr)> 0: 
            if (any(alt< 0)) or (all(alt< 15)):
                sys.exit('Altitudes are either < 0 or in Kilometers, not meters.')
        if len(pres_arr)> 0: 
             if (any(pres< 0)) or (any(pres >1100)): # Check that 
                 sys.exit('|Pressures| are < 0 or > 1100 hPa. Check units.')
        
        # Check for & drop any rows that have NaNs/ Inf values/ tell people about them! 
        pd.options.mode.use_inf_as_na = True # Look fro Inf values too. 
        if df.isna().values.any() == True: 
            print('WARNING: Found NaNs in data. Dropping the following rows from the data set:')
            print(df[df.isna().any(axis=1)])
            df.dropna(axis=0, how='any', inplace=True) 
        
        # Check for any duplicate rows that drop them/ tell people.
        if df.duplicated().values.any() == True: 
            print('WARNING: Found duplicate rows in data.'+
                  'If you wish to drop them please set keyword drop_dupes=True.')
            print(df[df.duplicated()])
            if drop_dupes== True: 
                df=df.drop_duplicates(ignore_index=True) 
        
        #==========================================================================
        # Format the numerical strings so that they're in GEOS-Chem's expected format!
        #==========================================================================
        df.POINT=df.POINT.map('{: >5}'.format) # poitns are len 5 char strings only
        df.LAT=df.LAT.map('{:7.2f}'.format) # Got these vals used for each var from 
        df.LON=df.LON.map('{:7.2f}'.format) # planeflight_mod.f90 i/o checks.
        df.OBS=df.OBS.map('{:10.3f}'.format)
        if len(pres_arr)> 0:
            df.PRESS=df.PRESS.map('{:7.2f}'.format)
        elif len(alt_arr)> 0:
            df.ALT=df.ALT.map('{:7.2f}'.format)
        
        if overwrite==True: #Delete existing file with this name. Otherwise it appends...
            if os.path.isfile(outpath+filename):
                os.remove(outpath+filename)
        
        # Write the header lines that GEOS Chem expects to the file, considering format.
        header= '{strr: >6}'.format(strr='POINT ')+ \
            '{strr: >7}'.format(strr='TYPE ')+ \
            '{strr: >11}'.format(strr='DD-MM-YYYY ')+ \
            '{strr: >6}'.format(strr='HH:MM ')+\
            '{strr: >8}'.format(strr='LAT ')+\
            '{strr: >8}'.format(strr='LON ')+\
            '{strr: >8}'.format(strr=vert_header) +\
            '{strr: >10}'.format(strr='OBS')
        
        #======================================================================
        #            Begin writing the planeflight.dat text file
        #=====================================================================
        # This is just a list containing our header lines in the right order... 
        textList = [filename, username, today, spacer, ntracers, spacer] + \
                    tracer_list + [spacer, title, spacer, header]
        
        # Open the output file and write headers line by line. 
        outF = open(outpath+filename, "w")
        for line in textList:
            outF.write(line)
            outF.write("\n")
            
        # Write the lat/lon/time/ pres data a temp file, using ASCII encoding!
        # Pandas default in Python 3 uses UTF-8 encoding, which GEOS-Chem can't read.
        df.to_csv(outpath+filename+'_0', header=False, index=None, sep=' ', mode='a',
                  encoding='ascii') 
    
        # Annoyingly if we use a space as a delimiter, we get werid quotation m
        # marks around strings, soooo then we'll open the temp file, read line 
        # by line, and take out the quotation marks, and write that to the 
        # actual output file. 
        fin = open(outpath+filename+'_0', 'r')
        Lines = fin.readlines()
        for line in Lines:
            outF.write(line.replace('"','' ))

        outF.close() # Close the output file. 
        fin.close() # Close the tempororay file 
        os.remove(outpath+filename+'_0') # And delete the temp file. 

        print('Output saved at: '+ outpath + filename) # tell where output is saved.
