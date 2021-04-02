# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 11:27:21 2021

@author: Dr. Jessica D. Haskins 
"""

import os 

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
