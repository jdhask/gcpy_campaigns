# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 10:40:31 2021

@author: jhask
"""

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
                            typestr='aSENEX',
                            datetimes=dt,
                            lat_arr=sn.GpsLat.values, 
                            lon_arr=sn.GpsLon.values, 
                            pres_arr=sn.StaticPrs.values,
                            username= 'jdh',
                            drop_dupes = False,
                            diags = ['all'],
                            diags_minus=minus,
                            print_diag_options= False)
