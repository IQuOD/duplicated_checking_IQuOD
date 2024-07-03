#!/usr/bin/env python3

########################################################################
### @copyright Copyright (C) 2024 All Rights Reserved.
### @file  Create_DNA_Summary.py
### @brief 
### @version
###		Date	|	Author			|	Version		|	Description
### ------------|-------------------|---------------|---------------
### 2024-03-26	|                   |	1.0			|	Create
######################################################################


"""
This script used for pre-processing profile data (only supports WOD18 netCDF format files).
During the process:
--Numerical metadata are retained.
--String metadata are converted into numerical values by using the ASCII code to convert each letter (including spaces) and then add the ASCII code of each letter to obtain final values.
The results are stored in npz format file.

***For data that is not in WOD18 netCDF format, it needs to be rewritten as WOD18 format firstly.
The variables included in the WOD18 netCDF files can be checked in Table 3 in the README.md

we used the following metadata and secondly stat information to calculate the 'DNA' for each profiles
meta_names = ['WOD_unique_id','accession_number', 'dataset_id', 'lat', 'lon', 'year', 'month', 'day', 'probe type',
              'recorder', 'hour', 'minute', 'depth_number', 'maximum_depth', 'hasTemp', 'hasSalinity', 'hasOxygen',
              'hasChlonophyll', 'country_name', 'GMT_time', 'WMO_ID', 'dbase_orig', 'project_name', 'Platform',
              'ocean_vehicle', 'Institute', 'WOD_cruise_identifier', 'sum_temp', 'sum_salinity', 'sum_depth',
              'std_depth', 'std_temp', 'std_salinity', 'corr_temp_depth', 'corr_sal_depth']
The order in 'meta_names' CANNOT be modifed. Please strictly following this order. Keep NaN or set it as '' if this information is missing.

Read the Section 5 of README file for customize your own netCDf file.

Usage:
    Run this script and follow the prompt to enter the directory path containing netCDF files.
"""
import os
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
from datetime import datetime, timedelta
import warnings
import sys
import argparse

def validate_path(input_path):
    # Normalize the path
    normalized_path = os.path.normpath(input_path)

    # Check if the path exists and is a directory
    if not os.path.exists(normalized_path) or not os.path.isdir(normalized_path):
        return False

    return True


def read_netCDF_formatted_DNA_series(inputpath, outputpath):

    warnings.filterwarnings("ignore")

    # Path where the netCDF files are stored
    inputpath = os.path.normpath(os.path.abspath(inputpath))
    print(inputpath)


    # Get all file names in the directory that are not directories themselves
    filenames = [f for f in os.listdir(inputpath) if os.path.isfile(os.path.join(inputpath, f))]
    n_prof = len(filenames)

    # print(filenames)
    # Initialize data structures
    DNA_series = np.full((n_prof, 35), np.nan, dtype=np.float32)  # Similar to single(NaN(n_prof,35))
    meta_names = ['WOD_unique_id','accession_number', 'dataset_id', 'lat', 'lon', 'year', 'month', 'day', 'probe type',
                  'recorder', 'hour', 'minute', 'depth_number', 'maximum_depth', 'hasTemp', 'hasSalinity', 'hasOxygen',
                  'hasChlonophyll', 'country_name', 'GMT_time', 'WMO_ID', 'dbase_orig', 'project_name', 'Platform',
                  'ocean_vehicle', 'Institute', 'WOD_cruise_identifier', 'sum_temp', 'sum_salinity', 'sum_depth',
                  'std_depth', 'std_temp', 'std_salinity', 'corr_temp_depth', 'corr_sal_depth']
    txt = [[''] * 35 for _ in range(n_prof)]

    # Process each file
    for idx, filename in enumerate(filenames):
        # if(idx>=100):
        #     break
        print(f"Processing file {idx+1}/{n_prof}: {filename}")
        file_path = os.path.join(inputpath, filename)

        # Open the netCDF file
        try:
            f = Dataset(file_path, 'r')
        except Exception as e:
            print(f"Failed to open {filename}: {str(e)}")
            continue

        ### read 'depth' (z) attributes
        try:
            # Read 'z' (depth) data
            depth = f.variables['z'][:]
            # Apply conditions: set values outside the acceptable range to NaN
            depth = np.where((depth > 12000) | (depth < -10), np.nan, depth)

            # Calculate number of depth measurements
            depth_number = len(depth)  # Using len() since depth is a numpy array

            # Determine maximum depth, taking care only to consider valid (non-NaN) entries
            if np.any(~np.isnan(depth)):
                maximum_depth = np.nanmax(depth)
            else:
                maximum_depth = np.nan

            # Compute the sum and standard deviation of depth, rounding to four decimal places
            sum_depth = np.round(np.nansum(depth), 4)
            std_depth = np.round(np.nanstd(depth), 4)

            # Check if the calculated sum and standard deviation are NaN and handle if they are
            sum_depth = sum_depth if not np.isnan(sum_depth) else np.nan
            std_depth = std_depth if not np.isnan(std_depth) else np.nan
        except:
            # Handle the case where the 'z' variable is not found
            depth = np.nan
            sum_depth=np.nan
            std_depth=np.nan
            depth_number = np.nan
            maximum_depth = np.nan

        # Read 'Temperature' data
        try:
            temp = f.variables['Temperature'][:]

            # Apply conditions: mask values outside the acceptable range
            temp = np.where((temp > 50) | (temp < -2.5), np.nan, temp)

            # Filter out NaN values for further processing
            temp2 = temp[~np.isnan(temp)]
            depth2 = depth[~np.isnan(temp)]  # Assuming 'depth' is already defined and processed similarly

            # Check if there are valid temperature readings
            hasTemp = 1 if np.any(~np.isnan(temp2)) else 0

            # Compute the sum, standard deviation, and correlation coefficient, if applicable
            if hasTemp:
                sum_temp = np.round(np.nansum(temp2), 4)
                std_temp = np.round(np.nanstd(temp2), 4)
                # Calculate correlation if both arrays have non-NaN data and at least two data points
                if len(temp2) > 1 and len(depth2) > 1:
                    cor_temp_depth = np.round(np.corrcoef(temp2, depth2)[0, 1], 5)
                else:
                    cor_temp_depth = np.nan
            else:
                sum_temp = std_temp = cor_temp_depth = np.nan

            # Handle edge cases explicitly
            if np.isnan(cor_temp_depth):
                cor_temp_depth = np.nan
            if np.isnan(sum_temp) or sum_temp == 0.0:
                sum_temp = np.nan
            if np.isnan(std_temp) or std_temp == 0.0:
                std_temp = np.nan

        except:
            # Handle the case where the 'Temperature' variable is not available
            hasTemp = 0
            sum_temp = std_temp = cor_temp_depth = np.nan

       # Read 'Salinity' data
        try:
            sal = f.variables['Salinity'][:]

            # Apply conditions: mask values outside the acceptable range
            sal = np.where((sal > 45) | (sal < -1), np.nan, sal)

            # Filter out NaN values for further processing
            sal2 = sal[~np.isnan(sal)]
            depth2 = depth[~np.isnan(sal)]  # Assuming 'depth' is already defined and processed similarly

            # Check if there are valid salinity readings
            hasSalinity = 1 if np.any(~np.isnan(sal2)) else 0

            # Compute the sum, standard deviation, and correlation coefficient, if applicable
            if hasSalinity:
                sum_sal = np.round(np.nansum(sal2), 4)
                std_sal = np.round(np.nanstd(sal2), 4)
                # Calculate correlation if both arrays have non-NaN data and at least two data points
                if len(sal2) > 1 and len(depth2) > 1:
                    cor_sal_depth = np.round(np.corrcoef(sal2, depth2)[0, 1], 5)
                else:
                    cor_sal_depth = np.nan
            else:
                sum_sal = std_sal = cor_sal_depth = np.nan

            # Handle edge cases explicitly
            if np.isnan(cor_sal_depth):
                cor_temp_depth = np.nan
            if np.isnan(sum_sal) or sum_sal == 0.0:
                sum_sal = np.nan
            if np.isnan(std_temp) or std_temp == 0.0:
                std_sal = np.nan
        except:
            # Handle the case where the 'Temperature' variable is not available
            hasSalinity = 0
            sum_sal = std_sal = cor_sal_depth = np.nan

        hasOxygen = 1 if 'Oxygen' in f.variables else 0
        hasChlorophyll = 1 if 'Chlorophyll' in f.variables else 0

        # Try to read various attributes and variables
        try:
            wod_unique_id = f.variables['wod_unique_cast'][:]
        except:
            wod_unique_id = np.nan

        # Read geographic coordinates
        try:
            lat = np.round(f.variables['lat'][:], 4)
            lon = np.round(f.variables['lon'][:], 4)
        except:
            lat = lon = np.nan

        # Read and process date and time
        # try:
        time_var=f.variables['time']
        dtime = nc.num2date(time_var[-1],time_var.units)
        year = dtime.year
        month = dtime.month
        day = dtime.day
        hour = dtime.hour
        minute = dtime.minute

        # Read depth and temperature, handle invalid data
        try:
            depth = f.variables['z'][:]
            depth = np.where((depth > 12000) | (depth < -10), np.nan, depth)
            depth_number = len(depth)
            maximum_depth = depth[-1] if len(depth) > 0 else np.nan
            sum_depth = np.round(np.nansum(depth), 4)
            std_depth = np.round(np.nanstd(depth), 4)
        except:
            depth_number = maximum_depth = sum_depth = std_depth = np.nan

        try:
            temp = f.variables['Temperature'][:]
            temp = np.where((temp > 40) | (temp < -2.5), np.nan, temp)
            temp_filtered = temp[~np.isnan(temp)]
            depth_filtered = depth[~np.isnan(temp)]
            has_temp = 1 if len(temp_filtered) > 0 else 0
            sum_temp = np.round(np.nansum(temp_filtered), 4)
            std_temp = np.round(np.nanstd(temp_filtered), 4)
            cor_temp_depth = np.round(np.corrcoef(temp_filtered, depth_filtered)[0, 1], 5) if len(temp_filtered) > 1 else np.nan
        except:
            has_temp = sum_temp = std_temp = cor_temp_depth = np.nan

        # Other variables such as salinity can be processed similarly
        try:
            sal = f.variables['Salinity'][:]
            sal[(sal > 45) | (sal < -1)] = np.nan
            sal2 = sal[~np.isnan(sal)]
            depth2 = depth[~np.isnan(sal)]
            hasSalinity = 1 if sal2.size > 0 else 0
            sum_sal = np.round(np.nansum(sal2), 4)
            std_sal = np.round(np.nanstd(sal2), 4)
            cor_sal_depth = np.round(np.corrcoef(sal2, depth2)[0, 1], 5) if len(sal2) > 1 else np.nan
        except:
            hasSalinity = sum_sal = std_sal = cor_sal_depth = 0

        hasOxygen = 1 if 'Oxygen' in f.variables else 0
        hasChlonophyII = 1 if 'Chlorophyll' in f.variables else 0

        try:
            wod_unique_id = np.float64(f.variables['wod_unique_cast'][:])
        except:
            wod_unique_id = np.nan

        try:
            country_name = f.variables['country'][:]
            country_name = bytes(country_name[~country_name.mask]).decode('ascii')
        except:
            country_name = ''


        try:
            probe_type = f.variables['Temperature_Instrument'][:]
            probe_type = bytes(probe_type[~probe_type.mask]).decode('ascii')
        except:
            probe_type = ''

        try:
            need_z_fix = f.variables['need_z_fix'][:]
            need_z_fix = bytes(need_z_fix[~need_z_fix.mask]).decode('ascii')
        except:
            need_z_fix = ''

        try:
            recorder = f.variables['Recorder'][:]
            recorder = bytes(recorder[~recorder.mask]).decode('ascii')
        except:
            recorder = ''

        try:
            GMT_time = np.float64(f.variables['GMT_time'][:])
        except:
            GMT_time = np.nan

        try:
            WMO_id = np.float64(f.variables['WMO_ID'][:])
        except:
            WMO_id = np.nan

        # Reading and processing the 'dbase_orig' attribute
        try:
            dbase_orig = f.variables['dbase_orig'][:]
            dbase_orig = bytes(dbase_orig[~dbase_orig.mask]).decode('ascii')
        except:
            dbase_orig = ''

        # Reading and processing the 'Project' attribute
        try:
            project_name = f.variables['Project'][:]
            project_name = bytes(project_name[~project_name.mask]).decode('ascii')
        except:
            project_name = ''

        # Reading and processing the 'Platform' attribute
        try:
            platform = f.variables['platform'][:]
            platform = bytes(platform[~platform.mask]).decode('ascii')
            # print(platform)
        except:
            platform = ''

        # Reading and processing the 'Ocean_Vehicle' attribute
        try:
            ocean_vehicle = f.variables['Ocean_Vehicle'][:]
            ocean_vehicle = bytes(ocean_vehicle[~ocean_vehicle.mask]).decode('ascii')
            # print(ocean_vehicle)
        except:
            ocean_vehicle = ''

        # Reading 'Access_no'
        try:
            accession_number = f.variables['Access_no'][:]
        except:
            accession_number = np.nan

        # Reading and processing the 'Institute' attribute
        try:
            institute = f.variables['Institute'][:]
            institute = bytes(institute[~institute.mask]).decode('ascii') 
            # print('Institute ='+institute)       
        except:
            institute = ''

        # Reading and processing the 'WOD_cruise_identifier' attribute
        try:
            wod_cruise_identifier = f.variables['WOD_cruise_identifier'][:]
            wod_cruise_identifier = bytes(wod_cruise_identifier[~wod_cruise_identifier.mask]).decode('ascii')         
        except:
            wod_cruise_identifier = ''
        # print(wod_cruise_identifier)


        try:
            dataset_name = f.variables['dataset'][:]
            dataset_name = bytes(dataset_name[~dataset_name.mask]).decode('ascii')         
            # dataset_name = ''.join(chr(x) for x in f.variables['dataset'][:]).strip()
            dataset_name=dataset_name.lower()
            if any(sub in dataset_name for sub in ['bod', 'bottle', 'ocean station', 'osd', 'low-resolution', 'low resolution']):
                dataset_id = 1
            elif any(sub in dataset_name for sub in ['towed', 'uor', 'undulating']):
                dataset_id = 10
            elif any(sub in dataset_name for sub in ['ctd', 'xctd']):
                dataset_id = 2
            elif any(sub in dataset_name for sub in ['mbt', 'mechanica', 'mb']):
                dataset_id = 3
            elif any(sub in dataset_name for sub in ['xbt', 'xb', 'expendable']):
                dataset_id = 4
            elif 'sur' in dataset_name or 'surface' in dataset_name:
                dataset_id = 5
            elif any(sub in dataset_name for sub in ['apb', 'autonomous', 'animal']):
                dataset_id = 6
            elif any(sub in dataset_name for sub in ['mrb', 'moored', 'tao']):
                dataset_id = 7
            elif any(sub in dataset_name for sub in ['pfl', 'argo', 'profiling']):
                dataset_id = 8
            elif 'drb' in dataset_name or 'drifting' in dataset_name:
                dataset_id = 9
            elif 'gld' in dataset_name or 'glider' in dataset_name:
                dataset_id = 11
            elif 'dbt' in dataset_name:
                dataset_id = 12
            elif 'std' in dataset_name:
                dataset_id = 13
            elif 'microbt' in dataset_name:
                dataset_id = 14
            else:
                dataset_id = np.nan
        except:
            dataset_id = np.nan

        # Reading 'lat' and 'lon'
        try:
            latitude = np.round(f.variables['lat'][:], 4)
        except:
            latitude = np.nan
        try:
            longitude = np.round(f.variables['lon'][:], 4)
        except:
            ongitude = np.nan


        # Store the processed data
        ######  please make sure the 'position (order)' of each variables are consistent with the order of meta_names
        # txt[idx][0]=str(filename)
        txt[idx][8]=str(probe_type)
        txt[idx][9]=str(recorder)
        txt[idx][18]=str(country_name)
        txt[idx][21]=str(dbase_orig)
        txt[idx][22]=str(project_name)
        txt[idx][23]=str(platform)
        txt[idx][24]=str(ocean_vehicle)
        txt[idx][25]=str(institute)
        txt[idx][26]=str(wod_cruise_identifier)
        strings_columns_order=[8,9,18,21,22,23,24,25,26]

        DNA_series[idx,0:8]=[wod_unique_id, accession_number, dataset_id, latitude, longitude, year, month, day]
        DNA_series[idx,10:18]=[hour, minute, depth_number, maximum_depth, hasTemp, hasSalinity, hasOxygen,hasChlorophyll]
        DNA_series[idx,19:21]=[GMT_time,WMO_id]
        DNA_series[idx,27:35]=[sum_temp, sum_sal, sum_depth, std_depth, std_temp, std_sal, cor_temp_depth, cor_sal_depth]


        # Close the netCDF file
        f.close()


    # Delete WOD_unique_id
    DNA_series=DNA_series[:,1:]
    txt = [row[1:] for row in txt]
    del meta_names[0]   #Delete WOD_unique_id


    # Converts the string to the ASCILL code and sums
    variables_index_to_process = [x - 1 for x in strings_columns_order]
    for i in range(n_prof):
        for j in variables_index_to_process:
            if j < len(txt[i]): 
                # sum of all ACILL for each string varaible
                ASCII_sum = sum(ord(char) for char in txt[i][j] if char != ' ')
                DNA_series[i][j] = ASCII_sum


    ###### check output folders
    script_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
    # Get the parent directory of the current script's directory
    parent_directory = os.path.dirname(script_directory)
    # Define the path for the 'Input_files' directory within the parent directory

    # Check if the 'Input_files' directory exists, and create it if it does not
    if not os.path.exists(outputpath):
        os.makedirs(outputpath)

    # Save the data in .npz format
    output_filename=os.path.join(outputpath,'DNA_summary.npz')
    np.savez(output_filename, DNA_series=DNA_series, txt=txt,meta_names=meta_names,filenames=filenames)

    print('The DNA formatted file are output to current folder: '+outputpath+'\n')
    print('The DNA filename is: '+output_filename)


if __name__ == '__main__':

    OutputDir = os.path.dirname(os.path.abspath(__file__)) + "/../Input_files"
    InputDir = OutputDir + "/WOD18_sample_1995"

    parser = argparse.ArgumentParser(description='Create DNA Summary')
    parser.add_argument("-i", "--input", type=str, default=InputDir)
    parser.add_argument("-o", "--output", type=str, default=OutputDir)
    args = parser.parse_args()
    InputDir = args.input
    OutputDir = args.output

    # check Input/Output dir vaild
    isInputOK = validate_path(InputDir)
    isOutputOK = validate_path(OutputDir)
    if(not (isInputOK and isOutputOK)):
        print("The entered path is not valid. Please ensure the path is correct and try again.")
        raise Exception("Invalid InputDir or OutputDir!", InputDir, OutputDir)
    
    DNA_summary_filename = OutputDir + "/DNA_summary.npz"
    iAct = 1
    if os.path.exists(DNA_summary_filename):
        iAct = input("Update DNA Summary or not(1: Yes (default); 0: No): ")
        print(iAct)

    if (iAct == 1):
        read_netCDF_formatted_DNA_series(InputDir, OutputDir)
        print("DNA_summary.npz Complete !")
