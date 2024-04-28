#!/usr/bin/env python3
"""
This script is to run the crude screening check.
We set up 14 criteria (14 checks) to identify the maximum number of possible duplicates (see Figure 1 of Song et al., 2024).

Usage:
    Run this script and follow the prompt to enter the directory path containing netCDF files.
"""

import numpy as np
import pandas as pd
import os
import scipy.io as sio
from scipy.stats import zscore
try:
    import math_util_functions
except:
    from DC_OCEAN.util import math_util_functions as math_util_functions


def validate_file(input_path):
    # Normalize the path
    normalized_path = os.path.normpath(input_path)

    # Check if the fiile exists
    if not os.path.exists(normalized_path):
        return False

    return True




##### logical flow: see Figure 1 in the manuscrpt (standard 1)
def N02_1_DNA_check_standardize_line(DNA_series,filename_info):
    '''
    Standardize each line of data.
    Calculate the arithmetic average of each line and then compare which is closer.
    '''

    # Normalization: The data is normalized to data with a mean of 0 and a variance of 1
    DNA_series_copy = DNA_series.copy()

    ###### 判断WMO_ID在那一列？
    DNA_series_copy = np.delete(DNA_series_copy, 19, axis=1)  # Deleting the 20th column

    ###### standardizeData with mean as 0 and standard deviation as 1
    DNA_mapped = math_util_functions.standardizeData(np.transpose(DNA_series_copy), 0, 1)
    DNA_mapped=np.transpose(DNA_mapped)
    DNA_mapped[DNA_series_copy == 0] = 0

    # Calculate the arithmetic average of each line
    average_DNA = np.nanmean(DNA_mapped, axis=1)


    # Sort average_DNA in ascending order
    index = np.argsort(average_DNA)
    average_DNA = average_DNA[index]
    filename_info = filename_info[index][:]
    DNA_mapped = DNA_mapped[index, :]
    DNA_series = DNA_series[index, :]

    duplicate_filename_list=[]
    number_pairs = 0
    number_profiles = 0
    for i in range(len(average_DNA)):
        # print(i)
        number1 = average_DNA[i]
        difference = np.abs((number1 - average_DNA) / number1 * 100)
        difference[:i] = np.nan
        duplicate_number = np.sum(difference < 0.0001) #threshold value: 0.0001%
        if duplicate_number >= 2:
            difference[i] = np.nan
            id = np.array([i] + list(np.where(difference == np.nanmin(difference))[0]))
            DNA_series_small = DNA_series[id, :]
            
            if DNA_series[i, 1] == 7:  # MRB   #No.2 column denotes the instrument type
                continue
            
            # Depth or temperature or salinity are scaled or translation, go to output
            if np.any(np.abs(DNA_series_small[0, [32, 33]] - DNA_series_small[1, [32, 33]]) < 1e-4):
                duplicate_small_list=[]
                for m in id:
                    duplicate_small_list.append(filename_info[m])
                duplicate_filename_list.append(duplicate_small_list)

                number_pairs += 1
                number_profiles += duplicate_number
                continue
            
            fragment_same_number = np.nansum(np.abs(DNA_series_small[0, :] - DNA_series_small[1, :]) < 1e-5, axis=None)
            if fragment_same_number < 25:
                continue

            # Handling of specific probe types and conditions
            #If it is XBT CTD MBT BOT; the location difference is plus or minus 5 degrees within a month; the same probe----excludes navigation continuous observation
            #If type,platform, and vehicle are the same, but sum_temp,corr(temp,depth) are different, it is judged to be multiple observations on the same survey ship/platform on the same route
            if (DNA_series_small[0, 1] in {4, 2, 1, 3} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 22, 23, 25]] == DNA_series_small[1, [4, 5, 7, 22, 23, 25]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.099
                index3 = np.abs(DNA_series_small[0, 32] - DNA_series_small[1, 32]) > 0.001
                index4 = np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 5) and \
                         np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) > 1e-5)
                if index1 and index2 and index3 and index4:
                    continue
            
            # Exclude long-term continuous observation of fixed points/nearby points
            if (DNA_series_small[0, 1] in {1, 7, 5, 4} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 8, 21, 22, 23]] == DNA_series_small[1, [4, 5, 7, 8, 21, 22, 23]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.05
                index3 = np.abs(DNA_series_small[0, 28] - DNA_series_small[1, 28]) < 1e-5
                index4 = np.all(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 0.01)
                if index1 and index2 and index3 and index4:
                    continue

            duplicate_small_list=[]
            for m in id:
                duplicate_small_list.append(filename_info[m])
            duplicate_filename_list.append(duplicate_small_list)
            
            number_pairs += 1
            number_profiles += duplicate_number
    return duplicate_filename_list


def N02_2_DNA_check_standardize_depth_tem(DNA_series,filename_info):
    '''
    Focus on searching for potential duplicate pairs that are equal in sum_depth and sum_temperature.
    '''

    # Normalization: The data is normalized to data with a mean of 0 and a variance of 1
    DNA_series_copy = DNA_series.copy()
    DNA_series_temp_depth=DNA_series_copy[:,[26,28]]
    DNA_series_temp_depth[np.abs(DNA_series_temp_depth)>1e6]=np.nan


    # Calculate the arithmetic average of each line
    average_DNA = np.nanmean(DNA_series_temp_depth, axis=1)


    # Sort average_DNA in ascending order
    index = np.argsort(average_DNA)
    average_DNA = average_DNA[index]
    filename_info = filename_info[index][:]
    DNA_series = DNA_series[index, :]

    duplicate_filename_list=[]
    number_pairs = 0
    number_profiles = 0
    for i in range(len(average_DNA)):
        # print(i)
        number1 = average_DNA[i]
        difference = np.abs((number1 - average_DNA) / number1 * 100)
        difference[:i] = np.nan
        duplicate_number = np.sum(difference < 0.0001) #threshold value: 0.0001%
        if duplicate_number >= 2:
            difference[i] = np.nan
            id = np.array([i] + list(np.where(difference == np.nanmin(difference))[0]))
            DNA_series_small = DNA_series[id, :]
            
            if DNA_series[i, 1] == 7:  # MRB   #No.2 column denotes the instrument type
                continue
            
            if DNA_series[i, 11] <=3:  #depth number less than 3
                continue            
            
            fragment_same_number = np.nansum(np.abs(DNA_series_small[0, :] - DNA_series_small[1, :]) < 1e-5, axis=None)
            if fragment_same_number < 26:
                continue

            # Exclude sum_depth vary greatly
            if len(id) <= 1:
                continue
            valid_ids = []
            # Check each element to determine if it should be included based on your condition
            for m in range(1, len(id)):
                if np.abs(DNA_series_small[0, 28] - DNA_series_small[m, 28]) <= 1:
                    valid_ids.append(id[m])  # Add the index to valid_ids if the condition is met
            # Update 'id' to only include valid indices
            valid_ids.append(id[0])
            id = valid_ids

            # Check if the number of valid ids is one or less after filtering
            if len(id) <= 1:
                continue

            # Handling of specific probe types and conditions
            #If it is XBT CTD MBT BOT; the location difference is plus or minus 5 degrees within a month; the same probe----excludes navigation continuous observation
            #If type,platform, and vehicle are the same, but sum_temp,corr(temp,depth) are different, it is judged to be multiple observations on the same survey ship/platform on the same route
            if (DNA_series_small[0, 1] in {4, 2, 1, 3} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 22, 23, 25]] == DNA_series_small[1, [4, 5, 7, 22, 23, 25]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.099
                index3 = np.abs(DNA_series_small[0, 32] - DNA_series_small[1, 32]) > 0.001
                index4 = np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 5) and \
                         np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) > 1e-5)
                if index1 and index2 and index3 and index4:
                    continue
            
            # Exclude long-term continuous observation of fixed points/nearby points
            if (DNA_series_small[0, 1] in {1, 7, 5} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 8, 21, 22, 23]] == DNA_series_small[1, [4, 5, 7, 8, 21, 22, 23]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.05
                index3 = np.abs(DNA_series_small[0, 28] - DNA_series_small[1, 28]) < 1e-5
                index4 = np.all(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 0.01)
                if index1 and index2 and index3 and index4:
                    continue

            # Exclude sum_salinity unequal and not 0
            if (np.abs(DNA_series_small[0, 27] - DNA_series_small[1, 27]) > 1e-3) and (DNA_series_small[0, 27] > 1e-6):
                continue

            duplicate_small_list=[]
            for m in id:
                duplicate_small_list.append(filename_info[m])
            duplicate_filename_list.append(duplicate_small_list)
            
            number_pairs += 1
            number_profiles += duplicate_number
    return duplicate_filename_list



def N02_3_DNA_check_entropy_weight(DNA_series,filename_info):
    '''
    Normalize each line of data.
    Calculate the weighted average by entropy weight method and then compare which is closer.
    Use a portion of the metadata.
    '''
    # Normalization: The data is normalized to data with a mean of 0 and a variance of 1
    DNA_series_meta = DNA_series[:,0:26]
    ###### 判断WMO_ID在那一列
    DNA_series_meta = np.delete(DNA_series_meta, 19, axis=1)  # Deleting the 20th column

    ###### normalizeData to 0-1
    DNA_mapped = math_util_functions.normalizeData(np.transpose(DNA_series_meta), 0, 1)
    DNA_mapped=np.transpose(DNA_mapped)
    DNA_mapped[DNA_series_meta == 0] = 0

    #Use entropy weight method to calculate the weight
    [weight,_]=math_util_functions.entropy_weight(DNA_series_meta);

    # Calculate the weighted average
    average_DNA_single = np.full_like(DNA_mapped,np.nan,dtype=float)
    # Calculate the weighted value for each column
    for i in range(len(weight)):
        average_DNA_single[:, i] = DNA_mapped[:, i] * weight[i]
    average_DNA = np.nansum(average_DNA_single, axis=1)

    # print(average_DNA[:])
    # raise('here')

    # Sort average_DNA in ascending order
    index = np.argsort(average_DNA)
    average_DNA = average_DNA[index]
    filename_info = filename_info[index][:]
    DNA_mapped = DNA_mapped[index, :]
    DNA_series = DNA_series[index, :]
    DNA_series_meta=DNA_series_meta[index, :]

    duplicate_filename_list=[]
    number_pairs = 0
    number_profiles = 0
    for i in range(len(average_DNA)):
        # print(i)
        number1 = average_DNA[i]
        difference = np.abs((number1 - average_DNA) / number1 * 100)
        difference[:i] = np.nan
        duplicate_number = np.sum(difference < 0.0001) #threshold value: 0.0001%
        if duplicate_number >= 2:
            difference[i] = np.nan
            id = np.array([i] + list(np.where(difference == np.nanmin(difference))[0]))
            DNA_series_small = DNA_series[id, :]
            DNA_series_small_meta=DNA_series_meta[id, :]
            
            # Calculate how many similar fragments there are
            fragment_same_number = np.nansum(np.abs(DNA_series_small_meta[0, :] - DNA_series_small_meta[1, :]) < 1e-5)
            if DNA_series[i, 1] in [7,9,5]:  # MRB, SUR, MRB   #No.2 column denotes the instrument type
                if(fragment_same_number<23):
                    continue
            else:
                if(fragment_same_number<22):
                    continue
            
            # Handling of specific probe types and conditions
            #If it is XBT CTD MBT BOT; the location difference is plus or minus 5 degrees within a month; the same probe----excludes navigation continuous observation
            #If type,platform, and vehicle are the same, but sum_temp,corr(temp,depth) are different, it is judged to be multiple observations on the same survey ship/platform on the same route
            if (DNA_series_small[0, 1] in {4, 2, 1, 3} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 22, 23, 25]] == DNA_series_small[1, [4, 5, 7, 22, 23, 25]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.099 #sum_temp is different
                index3 = np.abs(DNA_series_small[0, 32] - DNA_series_small[1, 32]) > 0.001 #cor_temp_depth is different 
                index4 = np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 5) and \
                         np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) > 1e-5)
                if index1 and index2 and index3 and index4:
                    continue
            
            # Exclude long-term continuous observation of fixed points/nearby points
            if (DNA_series_small[0, 1] in {1, 7, 5} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 8, 21, 22, 23]] == DNA_series_small[1, [4, 5, 7, 8, 21, 22, 23]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.05
                index3 = np.abs(DNA_series_small[0, 28] - DNA_series_small[1, 28]) < 1e-5
                index4 = np.all(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 0.01)
                if index1 and index2 and index3 and index4:
                    continue

            duplicate_small_list=[]
            for m in id:
                duplicate_small_list.append(filename_info[m])
            duplicate_filename_list.append(duplicate_small_list)
            
            number_pairs += 1
            number_profiles += duplicate_number
    return duplicate_filename_list

def N02_4_DNA_check_entropy_weight_allinfo(DNA_series,filename_info):
    '''
    Normalize each line of data.
    Calculate the weighted average by entropy weight method and then compare which is closer.
    Use all the metadata.
    '''
    # Normalization: The data is normalized to data with a mean of 0 and a variance of 1
    DNA_series_copy = DNA_series.copy()

    DNA_series_copy = np.delete(DNA_series_copy, 19, axis=1)  # Deleting the 20th column (WMO_ID) as no useful information

    ###### normalizeData to 0-1
    DNA_mapped = math_util_functions.normalizeData(np.transpose(DNA_series_copy), 0, 1)
    DNA_mapped=np.transpose(DNA_mapped)
    DNA_mapped[DNA_series_copy == 0] = 0

    #Use entropy weight method to calculate the weight
    [weight,_]=math_util_functions.entropy_weight(DNA_series_copy);


    # Calculate the weighted average for each column
    average_DNA_single = np.full_like(DNA_mapped,np.nan,dtype=float)
    for i in range(len(weight)):
        average_DNA_single[:, i] = DNA_mapped[:, i] * weight[i]
    average_DNA = np.nansum(average_DNA_single, axis=1)

    # Sort average_DNA in ascending order
    index = np.argsort(average_DNA)
    average_DNA = average_DNA[index]
    filename_info = filename_info[index][:]
    DNA_mapped = DNA_mapped[index, :]
    DNA_series = DNA_series[index, :]

    duplicate_filename_list=[]
    number_pairs = 0
    number_profiles = 0
    for i in range(len(average_DNA)):
        # print(i)
        number1 = average_DNA[i]
        difference = np.abs((number1 - average_DNA) / number1 * 100)
        difference[:i] = np.nan
        duplicate_number = np.sum(difference < 0.0001) #threshold value: 0.0001%
        if duplicate_number >= 2:
            difference[i] = np.nan
            id = np.array([i] + list(np.where(difference == np.nanmin(difference))[0]))
            DNA_series_small = DNA_series[id, :]
            
            ## ignore MRB data
            if DNA_series[i, 1] == 7:  # MRB   #No.2 column denotes the instrument type
                continue
            
            if np.any(np.abs(DNA_series_small[0, [32, 33]] - DNA_series_small[1, [32, 33]]) < 1e-4):
                duplicate_small_list=[]
                for m in id:
                    duplicate_small_list.append(filename_info[m])
                duplicate_filename_list.append(duplicate_small_list)

                number_pairs += 1
                number_profiles += duplicate_number
                continue
            
            fragment_same_number = np.nansum(np.abs(DNA_series_small[0, :] - DNA_series_small[1, :]) < 1e-5, axis=None)
            if fragment_same_number < 27:
                continue

            
            # Handling of specific probe types and conditions
            #If it is XBT CTD MBT BOT; the location difference is plus or minus 5 degrees within a month; the same probe----excludes navigation continuous observation
            #If type,platform, and vehicle are the same, but sum_temp,corr(temp,depth) are different, it is judged to be multiple observations on the same survey ship/platform on the same route
            if (DNA_series_small[0, 1] in {4, 2, 1, 3} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 22, 23, 25]] == DNA_series_small[1, [4, 5, 7, 22, 23, 25]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.099 #sum_temp is different
                index3 = np.abs(DNA_series_small[0, 32] - DNA_series_small[1, 32]) > 0.001 #cor_temp_depth is different 
                index4 = np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 5) and \
                         np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) > 1e-5)
                if index1 and index2 and index3 and index4:
                    continue
            
            # Exclude long-term continuous observation of fixed points/nearby points
            if (DNA_series_small[0, 1] in {1, 7, 5} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 8, 21, 22, 23]] == DNA_series_small[1, [4, 5, 7, 8, 21, 22, 23]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.05
                index3 = np.abs(DNA_series_small[0, 28] - DNA_series_small[1, 28]) < 1e-5
                index4 = np.all(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 0.01)
                if index1 and index2 and index3 and index4:
                    continue

            duplicate_small_list=[]
            for m in id:
                duplicate_small_list.append(filename_info[m])
            duplicate_filename_list.append(duplicate_small_list)
            
            number_pairs += 1
            number_profiles += duplicate_number
    return duplicate_filename_list

def N02_5_DNA_check_entropy_weight_noLATLON(DNA_series,filename_info):
    '''
    Normalize each line of data.
    Calculate the weighted average by entropy weight method and then compare which is closer.
    Ignore latitude and longitude information, the goal is to find duplicate pairs that may have been manipulated in latitude and longitude.
    '''

    # Normalization: The data is normalized to data with a mean of 0 and a variance of 1
    # delete latitude, longitude and WMO ID information
    indices = np.r_[0:2, 4:19, 20:34]
    DNA_series_copy = DNA_series.copy()
    DNA_series_copy = DNA_series_copy[:,indices]

    ###### normalizeData to 0-1
    DNA_mapped = math_util_functions.normalizeData(np.transpose(DNA_series_copy), 0, 1)
    DNA_mapped=np.transpose(DNA_mapped)
    DNA_mapped[DNA_series_copy == 0] = 0

    #Use entropy weight method to calculate the weight
    [weight,_]=math_util_functions.entropy_weight(DNA_series_copy);


    # Calculate the weighted average for each column
    average_DNA_single = np.full_like(DNA_mapped,np.nan,dtype=float)
    for i in range(len(weight)):
        average_DNA_single[:, i] = DNA_mapped[:, i] * weight[i]
    average_DNA = np.nansum(average_DNA_single, axis=1)

    # Sort average_DNA in ascending order
    index = np.argsort(average_DNA)
    average_DNA = average_DNA[index]
    filename_info = filename_info[index][:]
    DNA_mapped = DNA_mapped[index, :]
    DNA_series = DNA_series[index, :]
    DNA_series_copy=DNA_series_copy[index, :]

    duplicate_filename_list=[]
    number_pairs = 0
    number_profiles = 0
    for i in range(len(average_DNA)):
        # print(i)
        number1 = average_DNA[i]
        difference = np.abs((number1 - average_DNA) / number1 * 100)
        difference[:i] = np.nan
        duplicate_number = np.sum(difference < 0.0001) #threshold value: 0.0001%
        if duplicate_number >= 2:
            difference[i] = np.nan
            id = np.array([i] + list(np.where(difference == np.nanmin(difference))[0]))
            DNA_series_small = DNA_series[id, :]
            
            ## ignore MRB data
            if DNA_series[i, 1] == 7:  # MRB   #No.2 column denotes the instrument type
                continue
            
            fragment_same_number = np.nansum(np.abs(DNA_series_small[0, :] - DNA_series_small[1, :]) < 1e-5, axis=None)
            if fragment_same_number < 26:
                continue

            
            # Handling of specific probe types and conditions
            #If it is XBT CTD MBT BOT; the location difference is plus or minus 5 degrees within a month; the same probe----excludes navigation continuous observation
            #If type,platform, and vehicle are the same, but sum_temp,corr(temp,depth) are different, it is judged to be multiple observations on the same survey ship/platform on the same route
            if (DNA_series_small[0, 1] in {4, 2, 1, 3} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 22, 23, 25]] == DNA_series_small[1, [4, 5, 7, 22, 23, 25]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.099 #sum_temp is different
                index3 = np.abs(DNA_series_small[0, 32] - DNA_series_small[1, 32]) > 0.001 #cor_temp_depth is different 
                index4 = np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 5) and \
                         np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) > 1e-5)
                if index1 and index2 and index3 and index4:
                    continue
            
            # Exclude long-term continuous observation of fixed points/nearby points
            if (DNA_series_small[0, 1] in {1, 7, 5} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 8, 21, 22, 23]] == DNA_series_small[1, [4, 5, 7, 8, 21, 22, 23]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.05
                index3 = np.abs(DNA_series_small[0, 28] - DNA_series_small[1, 28]) < 1e-5
                index4 = np.all(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 0.01)
                if index1 and index2 and index3 and index4:
                    continue

            duplicate_small_list=[]
            for m in id:
                duplicate_small_list.append(filename_info[m])
            duplicate_filename_list.append(duplicate_small_list)
            
            number_pairs += 1
            number_profiles += duplicate_number
    return duplicate_filename_list

def N02_6_DNA_check_entropy_weight_noDepth(DNA_series,filename_info):
    '''
    Normalize each line of data.
    Calculate the weighted average by entropy weight method and then compare which is closer.
    Ignore depth number and maximum depth information, the goal is to find duplicate pairs that may have been manipulated in depth number and maximum depth.
    '''

    # delete depth_number, maximum depth and WMO ID
    DNA_series_copy = DNA_series.copy()
    DNA_series_meta = DNA_series_copy[:,0:26]
    DNA_series_meta = np.delete(DNA_series_meta, [11,12,19], axis=1) 

    ###### normalizeData to 0-1
    DNA_mapped = math_util_functions.normalizeData(np.transpose(DNA_series_meta), 0, 1)
    DNA_mapped=np.transpose(DNA_mapped)
    DNA_mapped[DNA_series_meta == 0] = 0

    #Use entropy weight method to calculate the weight
    [weight,_]=math_util_functions.entropy_weight(DNA_series_meta);


    # Calculate the weighted average for each column
    average_DNA_single = np.full_like(DNA_mapped,np.nan,dtype=float)
    for i in range(len(weight)):
        average_DNA_single[:, i] = DNA_mapped[:, i] * weight[i]
    average_DNA = np.nansum(average_DNA_single, axis=1)

    # Sort average_DNA in ascending order
    index = np.argsort(average_DNA)
    average_DNA = average_DNA[index]
    filename_info = filename_info[index][:]
    DNA_mapped = DNA_mapped[index, :]
    DNA_series = DNA_series[index, :]
    DNA_series_meta=DNA_series_meta[index, :]

    duplicate_filename_list=[]
    number_pairs = 0
    number_profiles = 0
    for i in range(len(average_DNA)):
        number1 = average_DNA[i]
        difference = np.abs((number1 - average_DNA) / number1 * 100)
        difference[:i] = np.nan
        duplicate_number = np.sum(difference < 0.0001) #threshold value: 0.0001%
        if duplicate_number >= 2:
            difference[i] = np.nan
            id = np.array([i] + list(np.where(difference == np.nanmin(difference))[0]))
            DNA_series_small = DNA_series[id, :]
            DNA_series_small_meta=DNA_series_meta[id, :]
            
            
            fragment_same_number = np.nansum(np.abs(DNA_series_small_meta[0, :] - DNA_series_small_meta[1, :]) < 1e-5, axis=None)
            if DNA_series[i, 1] in [7,9,5]:  #MRB SUR DRB
                if fragment_same_number < 23:
                    continue
            else:
                if fragment_same_number < 20:
                    continue
            
            # Handling of specific probe types and conditions
            # Exclude long-term continuous observation of fixed points/nearby points
            if (DNA_series_small[0, 1] in {1, 7, 5} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 8, 21, 22, 23]] == DNA_series_small[1, [4, 5, 7, 8, 21, 22, 23]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.05
                index3 = np.abs(DNA_series_small[0, 28] - DNA_series_small[1, 28]) < 1e-5
                index4 = np.all(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 0.01)
                if index1 and index2 and index3 and index4:
                    continue

            duplicate_small_list=[]
            for m in id:
                duplicate_small_list.append(filename_info[m])
            duplicate_filename_list.append(duplicate_small_list)
            
            number_pairs += 1
            number_profiles += duplicate_number
    return duplicate_filename_list


def N02_7_DNA_check_normalize_column(DNA_series,filename_info):
    '''
    Normalize each column of data.
    Calculate the arithmetic average of each line and then compare which is closer.
    '''

    # Normalization: The data is normalized to data with a mean of 0 and a variance of 1
    DNA_series_copy = DNA_series.copy()

    ###### 判断WMO_ID在那一列？
    DNA_series_copy = np.delete(DNA_series_copy, 19, axis=1)  # Deleting the 20th column

    ###### normalizeData with mean as 0 and standard deviation as 1
    DNA_mapped = math_util_functions.normalizeData(DNA_series_copy, 0, 1)
    DNA_mapped[:,4]=0
    DNA_mapped[DNA_series_copy == 0] = 0

    # Calculate the arithmetic average of each line
    average_DNA = np.nanmean(DNA_mapped, axis=1)


    # Sort average_DNA in ascending order
    index = np.argsort(average_DNA)
    average_DNA = average_DNA[index]
    filename_info = filename_info[index][:]
    DNA_mapped = DNA_mapped[index, :]
    DNA_series = DNA_series[index, :]

    duplicate_filename_list=[]
    number_pairs = 0
    number_profiles = 0
    for i in range(len(average_DNA)):
        # print(i)
        number1 = average_DNA[i]
        difference = np.abs((number1 - average_DNA) / number1 * 100)
        difference[:i] = np.nan
        duplicate_number = np.sum(difference < 0.0001) #threshold value: 0.0001%
        if duplicate_number >= 2:
            difference[i] = np.nan
            id = np.array([i] + list(np.where(difference == np.nanmin(difference))[0]))
            DNA_series_small = DNA_series[id, :]
            
            if DNA_series[i, 1] in [5,7,9]:  # MRB SUR DRB instruments   #No.2 column denotes the instrument type
                continue
            
            if np.any(np.abs(DNA_series_small[0, [32, 33]] - DNA_series_small[1, [32, 33]]) < 1e-4):
                duplicate_small_list=[]
                for m in id:
                    duplicate_small_list.append(filename_info[m])
                duplicate_filename_list.append(duplicate_small_list)

                number_pairs += 1
                number_profiles += duplicate_number
                continue
            
            fragment_same_number = np.nansum(np.abs(DNA_series_small[0, :] - DNA_series_small[1, :]) < 1e-5, axis=None)
            if fragment_same_number < 25:
                continue

            # Handling of specific probe types and conditions
            #If it is XBT CTD MBT BOT; the location difference is plus or minus 5 degrees within a month; the same probe----excludes navigation continuous observation
            #If type,platform, and vehicle are the same, but sum_temp,corr(temp,depth) are different, it is judged to be multiple observations on the same survey ship/platform on the same route
            if (DNA_series_small[0, 1] in {4, 2, 1, 3} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 22, 23, 25]] == DNA_series_small[1, [4, 5, 7, 22, 23, 25]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.099
                index3 = np.abs(DNA_series_small[0, 32] - DNA_series_small[1, 32]) > 0.001
                index4 = np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 5) and \
                         np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) > 1e-5)
                if index1 and index2 and index3 and index4:
                    continue
            
            # Exclude long-term continuous observation of fixed points/nearby points
            if (DNA_series_small[0, 1] in {1, 7, 5} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 8, 21, 22, 23]] == DNA_series_small[1, [4, 5, 7, 8, 21, 22, 23]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.05
                index3 = np.abs(DNA_series_small[0, 28] - DNA_series_small[1, 28]) < 1e-5
                index4 = np.all(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 0.01)
                if index1 and index2 and index3 and index4:
                    continue

            duplicate_small_list=[]
            for m in id:
                duplicate_small_list.append(filename_info[m])
            duplicate_filename_list.append(duplicate_small_list)
            
            number_pairs += 1
            number_profiles += duplicate_number
    return duplicate_filename_list


def N02_8_DNA_check_normalize_depth_tem(DNA_series,filename_info):
    '''
    Focus on searching for potential duplicate pairs that are equal in sum_depth and sum_temperature.
    Normalize each column of data.
    '''
    # Normalization: The data is normalized to data with a mean of 0 and a variance of 1
    DNA_series_copy = DNA_series.copy()

    DNA_series_temp_depth=DNA_series_copy[:,[26,28]]

    DNA_series_temp_depth[np.abs(DNA_series_temp_depth)>1e6]=np.nan


    #using mapminmax normalize each column of data
    DNA_series_temp_depth=math_util_functions.normalizeData(DNA_series_temp_depth,0,1);

    # Calculate the arithmetic average of each line
    average_DNA = np.nanmean(DNA_series_temp_depth, axis=1)


    # Sort average_DNA in ascending order
    index = np.argsort(average_DNA)
    average_DNA = average_DNA[index]
    filename_info = filename_info[index][:]
    DNA_series = DNA_series[index, :]

    duplicate_filename_list=[]
    number_pairs = 0
    number_profiles = 0
    for i in range(len(average_DNA)):
        # print(i)
        number1 = average_DNA[i]
        difference = np.abs((number1 - average_DNA) / number1 * 100)
        difference[:i] = np.nan
        duplicate_number = np.sum(difference < 0.0001) #threshold value: 0.0001%
        if duplicate_number >= 2:
            difference[i] = np.nan
            id = np.array([i] + list(np.where(difference == np.nanmin(difference))[0]))
            DNA_series_small = DNA_series[id, :]
            
            if DNA_series[i, 1] in [5,7,9]:  # MRB SUR DRB  #No.2 column denotes the instrument type
                continue
            
            if DNA_series[i, 11] <=3:  #depth number less than 3
                continue            
            
            fragment_same_number = np.nansum(np.abs(DNA_series_small[0, :] - DNA_series_small[1, :]) < 1e-5, axis=None)
            if fragment_same_number < 26:
                continue

            # Exclude sum_depth vary greatly
            if len(id) <= 1:
                continue
            valid_ids = []
            # Check each element to determine if it should be included based on your condition
            for m in range(1, len(id)):
                if np.abs(DNA_series_small[0, 28] - DNA_series_small[m, 28]) <= 1:
                    valid_ids.append(id[m])  # Add the index to valid_ids if the condition is met
            # Update 'id' to only include valid indices
            valid_ids.append(id[0])
            id = valid_ids

            # Check if the number of valid ids is one or less after filtering
            if len(id) <= 1:
                continue

            # Handling of specific probe types and conditions
            #If it is XBT CTD MBT BOT; the location difference is plus or minus 5 degrees within a month; the same probe----excludes navigation continuous observation
            #If type,platform, and vehicle are the same, but sum_temp,corr(temp,depth) are different, it is judged to be multiple observations on the same survey ship/platform on the same route
            if (DNA_series_small[0, 1] in {4, 2, 1, 3} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 22, 23, 25]] == DNA_series_small[1, [4, 5, 7, 22, 23, 25]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.099
                index3 = np.abs(DNA_series_small[0, 32] - DNA_series_small[1, 32]) > 0.001
                index4 = np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 5) and \
                         np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) > 1e-5)
                if index1 and index2 and index3 and index4:
                    continue
            
            # Exclude long-term continuous observation of fixed points/nearby points
            if (DNA_series_small[0, 1] in {1, 7, 5} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 8, 21, 22, 23]] == DNA_series_small[1, [4, 5, 7, 8, 21, 22, 23]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.05
                index3 = np.abs(DNA_series_small[0, 28] - DNA_series_small[1, 28]) < 1e-5
                index4 = np.all(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 0.01)
                if index1 and index2 and index3 and index4:
                    continue

            # Exclude sum_salinity unequal and not 0
            if (np.abs(DNA_series_small[0, 27] - DNA_series_small[1, 27]) > 1e-3) and (DNA_series_small[0, 27] > 1e-6):
                continue

            duplicate_small_list=[]
            for m in id:
                duplicate_small_list.append(filename_info[m])
            duplicate_filename_list.append(duplicate_small_list)
            
            number_pairs += 1
            number_profiles += duplicate_number
    return duplicate_filename_list



def N02_9_DNA_check_entropy_normalizeData_column(DNA_series,filename_info):
    '''
    Normalize each column of data.
    Calculate the weighted average by entropy weight method and then compare which is closer.
    Use a portion of the metadata.
    '''
    # Normalization: The data is normalized to data with a mean of 0 and a variance of 1
    DNA_series_meta = DNA_series[:,0:26]

    DNA_series_meta = np.delete(DNA_series_meta, 19, axis=1)  # Deleting the 20th column

    ###### normalizeData to 0-1
    DNA_mapped = math_util_functions.normalizeData(DNA_series_meta, 0, 1)
    DNA_mapped[:,4]=0; #delete year column
    DNA_mapped[DNA_series_meta == 0] = 0

    #Use entropy weight method to calculate the weight
    [weight,_]=math_util_functions.entropy_weight(DNA_series_meta);

    # Calculate the weighted average
    average_DNA_single = np.full_like(DNA_mapped,np.nan,dtype=float)
    # Calculate the weighted value for each column
    for i in range(len(weight)):
        average_DNA_single[:, i] = DNA_mapped[:, i] * weight[i]
    average_DNA = np.nansum(average_DNA_single, axis=1)

    # print(average_DNA[:])
    # raise('here')

    # Sort average_DNA in ascending order
    index = np.argsort(average_DNA)
    average_DNA = average_DNA[index]
    filename_info = filename_info[index][:]
    DNA_mapped = DNA_mapped[index, :]
    DNA_series = DNA_series[index, :]
    DNA_series_meta=DNA_series_meta[index, :]

    duplicate_filename_list=[]
    number_pairs = 0
    number_profiles = 0
    for i in range(len(average_DNA)):
        # print(i)
        number1 = average_DNA[i]
        difference = np.abs((number1 - average_DNA) / number1 * 100)
        difference[:i] = np.nan
        duplicate_number = np.sum(difference < 0.0001) #threshold value: 0.0001%
        if duplicate_number >= 2:
            difference[i] = np.nan
            id = np.array([i] + list(np.where(difference == np.nanmin(difference))[0]))
            DNA_series_small = DNA_series[id, :]
            DNA_series_small_meta=DNA_series_meta[id, :]
            
            # Calculate how many similar fragments there are
            fragment_same_number = np.nansum(np.abs(DNA_series_small_meta[0, :] - DNA_series_small_meta[1, :]) < 1e-5)
            if DNA_series[i, 1] in [7,9,5]:  # MRB, SUR, MRB   #No.2 column denotes the instrument type
                if(fragment_same_number<23):
                    continue
            else:
                if(fragment_same_number<22):
                    continue
            
            # Handling of specific probe types and conditions
            #If it is XBT CTD MBT BOT; the location difference is plus or minus 5 degrees within a month; the same probe----excludes navigation continuous observation
            #If type,platform, and vehicle are the same, but sum_temp,corr(temp,depth) are different, it is judged to be multiple observations on the same survey ship/platform on the same route
            if (DNA_series_small[0, 1] in {4, 2, 1, 3} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 22, 23, 25]] == DNA_series_small[1, [4, 5, 7, 22, 23, 25]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.099 #sum_temp is different
                index3 = np.abs(DNA_series_small[0, 32] - DNA_series_small[1, 32]) > 0.001 #cor_temp_depth is different 
                index4 = np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 5) and \
                         np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) > 1e-5)
                if index1 and index2 and index3 and index4:
                    continue
            
            # Exclude long-term continuous observation of fixed points/nearby points
            if (DNA_series_small[0, 1] in {1, 7, 5} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 8, 21, 22, 23]] == DNA_series_small[1, [4, 5, 7, 8, 21, 22, 23]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.05
                index3 = np.abs(DNA_series_small[0, 28] - DNA_series_small[1, 28]) < 1e-5
                index4 = np.all(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 0.01)
                if index1 and index2 and index3 and index4:
                    continue

            duplicate_small_list=[]
            for m in id:
                duplicate_small_list.append(filename_info[m])
            duplicate_filename_list.append(duplicate_small_list)
            
            number_pairs += 1
            number_profiles += duplicate_number
    return duplicate_filename_list

def N02_10_DNA_check_normalizeData_weight_column_allinfo(DNA_series,filename_info):
    '''
    Normalize each column of data.
    Calculate the weighted average by entropy weight method and then compare which is closer.
    Use all the metadata.
    '''
    # Normalization: The data is normalized to data with a mean of 0 and a variance of 1
    DNA_series_copy = DNA_series.copy()
    ###### 判断WMO_ID在那一列
    DNA_series_copy = np.delete(DNA_series_copy, 19, axis=1)  # Deleting the 20th column

    ###### normalizeData to 0-1
    DNA_mapped = math_util_functions.normalizeData(DNA_series_copy, 0, 1)
    DNA_mapped[:,4]=0
    DNA_mapped[DNA_series_copy == 0] = 0

    #Use entropy weight method to calculate the weight
    [weight,_]=math_util_functions.entropy_weight(DNA_series_copy);


    # Calculate the weighted average for each column
    average_DNA_single = np.full_like(DNA_mapped,np.nan,dtype=float)
    for i in range(len(weight)):
        average_DNA_single[:, i] = DNA_mapped[:, i] * weight[i]
    average_DNA = np.nansum(average_DNA_single, axis=1)

    # Sort average_DNA in ascending order
    index = np.argsort(average_DNA)
    average_DNA = average_DNA[index]
    filename_info = filename_info[index][:]
    DNA_mapped = DNA_mapped[index, :]
    DNA_series = DNA_series[index, :]

    duplicate_filename_list=[]
    number_pairs = 0
    number_profiles = 0
    for i in range(len(average_DNA)):
        # print(i)
        number1 = average_DNA[i]
        difference = np.abs((number1 - average_DNA) / number1 * 100)
        difference[:i] = np.nan
        duplicate_number = np.sum(difference < 0.0001) #threshold value: 0.0001%
        if duplicate_number >= 2:
            difference[i] = np.nan
            id = np.array([i] + list(np.where(difference == np.nanmin(difference))[0]))
            DNA_series_small = DNA_series[id, :]
            
            ## ignore MRB data
            if DNA_series[i, 1] in [5,7,9]:  # SUR DRB MRB   #No.2 column denotes the instrument type
                continue
            
            if np.any(np.abs(DNA_series_small[0, [32, 33]] - DNA_series_small[1, [32, 33]]) < 1e-4):
                duplicate_small_list=[]
                for m in id:
                    duplicate_small_list.append(filename_info[m])
                duplicate_filename_list.append(duplicate_small_list)

                number_pairs += 1
                number_profiles += duplicate_number
                continue
            
            fragment_same_number = np.nansum(np.abs(DNA_series_small[0, :] - DNA_series_small[1, :]) < 1e-5, axis=None)
            if fragment_same_number < 27:
                continue

            
            # Handling of specific probe types and conditions
            #If it is XBT CTD MBT BOT; the location difference is plus or minus 5 degrees within a month; the same probe----excludes navigation continuous observation
            #If type,platform, and vehicle are the same, but sum_temp,corr(temp,depth) are different, it is judged to be multiple observations on the same survey ship/platform on the same route
            if (DNA_series_small[0, 1] in {4, 2, 1, 3} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 22, 23, 25]] == DNA_series_small[1, [4, 5, 7, 22, 23, 25]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.099 #sum_temp is different
                index3 = np.abs(DNA_series_small[0, 32] - DNA_series_small[1, 32]) > 0.001 #cor_temp_depth is different 
                index4 = np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 5) and \
                         np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) > 1e-5)
                if index1 and index2 and index3 and index4:
                    continue
            
            # Exclude long-term continuous observation of fixed points/nearby points
            if (DNA_series_small[0, 1] in {1, 7, 5} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 8, 21, 22, 23]] == DNA_series_small[1, [4, 5, 7, 8, 21, 22, 23]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.05
                index3 = np.abs(DNA_series_small[0, 28] - DNA_series_small[1, 28]) < 1e-5
                index4 = np.all(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 0.01)
                if index1 and index2 and index3 and index4:
                    continue

            duplicate_small_list=[]
            for m in id:
                duplicate_small_list.append(filename_info[m])
            duplicate_filename_list.append(duplicate_small_list)
            
            number_pairs += 1
            number_profiles += duplicate_number
    return duplicate_filename_list



def N02_11_DNA_check_normalize_column_entropy_weight(DNA_series,filename_info):
    '''
    Normalize each column of data.
    Calculate the weighted average by entropy weight method and then compare which is closer.
    Ignore latitude and longitude information, the goal is to find duplicate pairs that may have been manipulated in latitude and longitude.
    '''

    # Normalization: The data is normalized to data with a mean of 0 and a variance of 1
    # delete latitude, longitude and WMO ID information
    indices = np.r_[0:2, 4:19, 20:34]
    DNA_series_copy = DNA_series.copy()
    DNA_series_copy = DNA_series_copy[:,indices]

    #Use entropy weight method to calculate the weight
    [weight,_]=math_util_functions.entropy_weight(DNA_series_copy)

    # normalizeData to 0-1
    DNA_mapped = math_util_functions.normalizeData(DNA_series_copy, 0, 1)
    DNA_mapped[:,4]=0 
    DNA_mapped[DNA_series_copy == 0] = 0

    # Calculate the weighted average for each column
    average_DNA_single = np.full_like(DNA_mapped,np.nan,dtype=float)
    for i in range(len(weight)):
        average_DNA_single[:, i] = DNA_mapped[:, i] * weight[i]
    average_DNA = np.nansum(average_DNA_single, axis=1)

    # Sort average_DNA in ascending order
    index = np.argsort(average_DNA)
    average_DNA = average_DNA[index]
    filename_info = filename_info[index][:]
    DNA_mapped = DNA_mapped[index, :]
    DNA_series = DNA_series[index, :]
    DNA_series_copy=DNA_series_copy[index, :]

    duplicate_filename_list=[]
    number_pairs = 0
    number_profiles = 0
    for i in range(len(average_DNA)):
        # print(i)
        number1 = average_DNA[i]
        difference = np.abs((number1 - average_DNA) / number1 * 100)
        difference[:i] = np.nan
        duplicate_number = np.sum(difference < 0.0001) #threshold value: 0.0001%
        if duplicate_number >= 2:
            difference[i] = np.nan
            id = np.array([i] + list(np.where(difference == np.nanmin(difference))[0]))
            DNA_series_small = DNA_series[id, :]
            
            ## ignore MRB data
            if DNA_series[i, 1] in [5,7,9]:  # MRB   #No.2 column denotes the instrument type
                continue
            
            fragment_same_number = np.nansum(np.abs(DNA_series_small[0, :] - DNA_series_small[1, :]) < 1e-5, axis=None)
            if fragment_same_number < 26:
                continue

            
            # Handling of specific probe types and conditions
            #If it is XBT CTD MBT BOT; the location difference is plus or minus 5 degrees within a month; the same probe----excludes navigation continuous observation
            #If type,platform, and vehicle are the same, but sum_temp,corr(temp,depth) are different, it is judged to be multiple observations on the same survey ship/platform on the same route
            if (DNA_series_small[0, 1] in {4, 2, 1, 3} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 22, 23, 25]] == DNA_series_small[1, [4, 5, 7, 22, 23, 25]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.099 #sum_temp is different
                index3 = np.abs(DNA_series_small[0, 32] - DNA_series_small[1, 32]) > 0.001 #cor_temp_depth is different 
                index4 = np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 5) and \
                         np.any(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) > 1e-5)
                if index1 and index2 and index3 and index4:
                    continue
            
            # Exclude long-term continuous observation of fixed points/nearby points
            if (DNA_series_small[0, 1] in {1, 7, 5} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 8, 21, 22, 23]] == DNA_series_small[1, [4, 5, 7, 8, 21, 22, 23]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.05
                index3 = np.abs(DNA_series_small[0, 28] - DNA_series_small[1, 28]) < 1e-5
                index4 = np.all(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 0.01)
                if index1 and index2 and index3 and index4:
                    continue

            duplicate_small_list=[]
            for m in id:
                duplicate_small_list.append(filename_info[m])
            duplicate_filename_list.append(duplicate_small_list)
            
            number_pairs += 1
            number_profiles += duplicate_number
    return duplicate_filename_list


def N02_12_DNA_check_entropy_weight_column_noDepth(DNA_series,filename_info):
    '''
    Normalize each column of data.
    Calculate the weighted average by entropy weight method and then compare which is closer.
    Ignore depth number and maximum depth information, the goal is to find duplicate pairs that may have been manipulated in depth number and maximum depth.
    '''

    # delete depth_number, maximum depth and WMO ID
    DNA_series_copy = DNA_series.copy()
    DNA_series_meta = DNA_series_copy[:,0:26]
    DNA_series_meta = np.delete(DNA_series_meta, [11,12,19], axis=1) 

    ###### normalizeData to 0-1
    DNA_mapped = math_util_functions.normalizeData(DNA_series_meta, 0, 1)
    DNA_mapped[:,4]=0;
    DNA_mapped[DNA_series_meta == 0] = 0

    #Use entropy weight method to calculate the weight
    [weight,_]=math_util_functions.entropy_weight(DNA_series_meta);


    # Calculate the weighted average for each column
    average_DNA_single = np.full_like(DNA_mapped,np.nan,dtype=float)
    for i in range(len(weight)):
        average_DNA_single[:, i] = DNA_mapped[:, i] * weight[i]
    average_DNA = np.nansum(average_DNA_single, axis=1)

    # Sort average_DNA in ascending order
    index = np.argsort(average_DNA)
    average_DNA = average_DNA[index]
    filename_info = filename_info[index][:]
    DNA_mapped = DNA_mapped[index, :]
    DNA_series = DNA_series[index, :]
    DNA_series_meta=DNA_series_meta[index, :]

    duplicate_filename_list=[]
    number_pairs = 0
    number_profiles = 0
    for i in range(len(average_DNA)):
        number1 = average_DNA[i]
        difference = np.abs((number1 - average_DNA) / number1 * 100)
        difference[:i] = np.nan
        duplicate_number = np.sum(difference < 0.0001) #threshold value: 0.0001%
        if duplicate_number >= 2:
            difference[i] = np.nan
            id = np.array([i] + list(np.where(difference == np.nanmin(difference))[0]))
            DNA_series_small = DNA_series[id, :]
            DNA_series_small_meta=DNA_series_meta[id, :]
        
            fragment_same_number = np.nansum(np.abs(DNA_series_small_meta[0, :] - DNA_series_small_meta[1, :]) < 1e-5, axis=None)
            if DNA_series[i, 1] in [7,9,5]:  #MRB SUR DRB
                if fragment_same_number < 23:
                    continue
            else:
                if fragment_same_number < 20:
                    continue
            
            # Handling of specific probe types and conditions
            # Exclude long-term continuous observation of fixed points/nearby points
            if (DNA_series_small[0, 1] in {1, 7, 5} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 8, 21, 22, 23]] == DNA_series_small[1, [4, 5, 7, 8, 21, 22, 23]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.05
                index3 = np.abs(DNA_series_small[0, 28] - DNA_series_small[1, 28]) < 1e-5
                index4 = np.all(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 0.01)
                if index1 and index2 and index3 and index4:
                    continue

            duplicate_small_list=[]
            for m in id:
                duplicate_small_list.append(filename_info[m])
            duplicate_filename_list.append(duplicate_small_list)
            
            number_pairs += 1
            number_profiles += duplicate_number
    return duplicate_filename_list


def N02_13_DNA_check_PCA_90_allinfo(DNA_series,filename_info):
    '''
    Use principal component analysis to calculate the DNA of each profile
    The principal component contribution threshold is set to 90%
    '''
    # delete depth_number, maximum depth and WMO ID
    DNA_series_copy = DNA_series.copy()
    DNA_series_copy = np.delete(DNA_series_copy, 19, axis=1) 
    DNA_series_copy[np.isnan(DNA_series_copy)]=0

    ###### Use principal component analysis(PCA) to calculate the DNA of each profile
    DNA_mapped=math_util_functions.PCA_DNA_profiles(DNA_series_copy,0.9)

    #Use entropy weight method to calculate the weight
    [weight,_]=math_util_functions.entropy_weight(DNA_mapped);


    # Calculate the weighted average for each column
    average_DNA_single = np.full_like(DNA_mapped,np.nan,dtype=float)
    for i in range(len(weight)):
        average_DNA_single[:, i] = DNA_mapped[:, i] * weight[i]
    average_DNA = np.nansum(average_DNA_single, axis=1)

    # Sort average_DNA in ascending order
    index = np.argsort(average_DNA)
    average_DNA = average_DNA[index]
    filename_info = filename_info[index][:]
    DNA_mapped = DNA_mapped[index, :]
    DNA_series = DNA_series[index, :]

    duplicate_filename_list=[]
    number_pairs = 0
    number_profiles = 0
    for i in range(len(average_DNA)):
        number1 = average_DNA[i]
        difference = np.abs((number1 - average_DNA) / number1 * 100)
        difference[:i] = np.nan
        duplicate_number = np.sum(difference < 0.001) #threshold value: 0.001%
        if duplicate_number >= 2:
            difference[i] = np.nan
            id = np.array([i] + list(np.where(difference == np.nanmin(difference))[0]))
            DNA_series_small = DNA_series[id, :]
        
            if DNA_series[i, 1] in [7,9,5]:  #MRB SUR DRB
                continue  

            # Depth or temperature or salinity are scaled or translation, go to output
            if np.any(np.abs(DNA_series_small[0, [32, 33]] - DNA_series_small[1, [32, 33]]) < 1e-4):
                duplicate_small_list=[]
                for m in id:
                    duplicate_small_list.append(filename_info[m])
                duplicate_filename_list.append(duplicate_small_list)

                number_pairs += 1
                number_profiles += duplicate_number
                continue

            #Calculate how many similar fragments there are    
            fragment_same_number = np.nansum(np.abs(DNA_series_small[0, :] - DNA_series_small[1, :]) < 1e-5, axis=None)
            if fragment_same_number < 27:
                continue

            
            # Handling of specific probe types and conditions
            # Exclude long-term continuous observation of fixed points/nearby points
            if (DNA_series_small[0, 1] in {1, 7, 5} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 8, 21, 22, 23]] == DNA_series_small[1, [4, 5, 7, 8, 21, 22, 23]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.05
                index3 = np.abs(DNA_series_small[0, 28] - DNA_series_small[1, 28]) < 1e-5
                index4 = np.all(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 0.01)
                if index1 and index2 and index3 and index4:
                    continue

            duplicate_small_list=[]
            for m in id:
                duplicate_small_list.append(filename_info[m])
            duplicate_filename_list.append(duplicate_small_list)
            
            number_pairs += 1
            number_profiles += duplicate_number
    return duplicate_filename_list




def N02_14_DNA_check_PCA_95_allinfo(DNA_series,filename_info):
    '''
    Use principal component analysis to calculate the DNA of each profile
    The principal component contribution threshold is set to 95%
    '''
    # delete depth_number, maximum depth and WMO ID
    DNA_series_copy = DNA_series.copy()
    DNA_series_copy = np.delete(DNA_series_copy, 19, axis=1) 
    DNA_series_copy[np.isnan(DNA_series_copy)]=0

    ###### Use principal component analysis(PCA) to calculate the DNA of each profile
    DNA_mapped=math_util_functions.PCA_DNA_profiles(DNA_series_copy,0.95)

    #Use entropy weight method to calculate the weight
    [weight,_]=math_util_functions.entropy_weight(DNA_mapped);


    # Calculate the weighted average for each column
    average_DNA_single = np.full_like(DNA_mapped,np.nan,dtype=float)
    for i in range(len(weight)):
        average_DNA_single[:, i] = DNA_mapped[:, i] * weight[i]
    average_DNA = np.nansum(average_DNA_single, axis=1)

    # Sort average_DNA in ascending order
    index = np.argsort(average_DNA)
    average_DNA = average_DNA[index]
    filename_info = filename_info[index][:]
    DNA_mapped = DNA_mapped[index, :]
    DNA_series = DNA_series[index, :]

    duplicate_filename_list=[]
    number_pairs = 0
    number_profiles = 0
    for i in range(len(average_DNA)):
        number1 = average_DNA[i]
        difference = np.abs((number1 - average_DNA) / number1 * 100)
        difference[:i] = np.nan
        duplicate_number = np.sum(difference < 0.001) #threshold value: 0.001%
        if duplicate_number >= 2:
            difference[i] = np.nan
            id = np.array([i] + list(np.where(difference == np.nanmin(difference))[0]))
            DNA_series_small = DNA_series[id, :]
        
            if DNA_series[i, 1] in [7,9,5]:  #MRB SUR DRB
                continue  

            # Depth or temperature or salinity are scaled or translation, go to output
            if np.any(np.abs(DNA_series_small[0, [32, 33]] - DNA_series_small[1, [32, 33]]) < 1e-4):
                duplicate_small_list=[]
                for m in id:
                    duplicate_small_list.append(filename_info[m])
                duplicate_filename_list.append(duplicate_small_list)

                number_pairs += 1
                number_profiles += duplicate_number
                continue

            #Calculate how many similar fragments there are    
            fragment_same_number = np.nansum(np.abs(DNA_series_small[0, :] - DNA_series_small[1, :]) < 1e-5, axis=None)
            if fragment_same_number < 27:
                continue

            
            # Handling of specific probe types and conditions
            # Exclude long-term continuous observation of fixed points/nearby points
            if (DNA_series_small[0, 1] in {1, 7, 5} and DNA_series_small[1, 1] == DNA_series_small[0, 1]):
                index1 = np.all(DNA_series_small[0, [4, 5, 7, 8, 21, 22, 23]] == DNA_series_small[1, [4, 5, 7, 8, 21, 22, 23]])
                index2 = np.abs(DNA_series_small[0, 26] - DNA_series_small[1, 26]) > 0.05
                index3 = np.abs(DNA_series_small[0, 28] - DNA_series_small[1, 28]) < 1e-5
                index4 = np.all(np.abs(DNA_series_small[0, [2, 3]] - DNA_series_small[1, [2, 3]]) < 0.01)
                if index1 and index2 and index3 and index4:
                    continue

            duplicate_small_list=[]
            for m in id:
                duplicate_small_list.append(filename_info[m])
            duplicate_filename_list.append(duplicate_small_list)
            
            number_pairs += 1
            number_profiles += duplicate_number
    return duplicate_filename_list


def main(DNA_summary_filename):
    # Load *.npz data
    print('loading the DNA summary files....')
    data = np.load(DNA_summary_filename)
    DNA_series = data['DNA_series']
    meta_names = data['meta_names']
    filename_info = data['filenames']


    ####### N02_1: run Crude Scrren" Arithmetic mean: standardlizatioin check
    print('Running the Crude Screen check: the No.1 criteria check...')
    duplicate_filename_list_1=N02_1_DNA_check_standardize_line(DNA_series,filename_info)

    ##### run Crude Scrren" Arithmetic mean: standardlizatioin check
    print('Running the Crude Screen check: the No.2 criteria check...')
    duplicate_filename_list_2=N02_2_DNA_check_standardize_depth_tem(DNA_series,filename_info)

    ##### run Crude Scrren" Arithmetic mean: standardlizatioin check
    print('Running the Crude Screen check: the No.3 criteria check...')
    duplicate_filename_list_3=N02_3_DNA_check_entropy_weight(DNA_series,filename_info)

    print('Running the Crude Screen check: the No.4 criteria check...')
    duplicate_filename_list_4=N02_4_DNA_check_entropy_weight_allinfo(DNA_series,filename_info)

    print('Running the Crude Screen check: the No.5 criteria check...')
    duplicate_filename_list_5=N02_5_DNA_check_entropy_weight_noLATLON(DNA_series,filename_info)

    print('Running the Crude Screen check: the No.6 criteria check...')
    duplicate_filename_list_6=N02_6_DNA_check_entropy_weight_noDepth(DNA_series,filename_info)

    print('Running the Crude Screen check: the No.7 criteria check...')
    duplicate_filename_list_7=N02_7_DNA_check_normalize_column(DNA_series,filename_info)

    print('Running the Crude Screen check: the No.8 criteria check...')
    duplicate_filename_list_8=N02_8_DNA_check_normalize_depth_tem(DNA_series,filename_info)

    print('Running the Crude Screen check: the No.9 criteria check...')
    duplicate_filename_list_9=N02_9_DNA_check_entropy_normalizeData_column(DNA_series,filename_info)

    print('Running the Crude Screen check: the No.10 criteria check...')
    duplicate_filename_list_10=N02_10_DNA_check_normalizeData_weight_column_allinfo(DNA_series,filename_info)
    
    print('Running the Crude Screen check: the No.11 criteria check...')
    duplicate_filename_list_11=N02_11_DNA_check_normalize_column_entropy_weight(DNA_series,filename_info)

    print('Running the Crude Screen check: the No.12 criteria check...')
    duplicate_filename_list_12=N02_12_DNA_check_entropy_weight_column_noDepth(DNA_series,filename_info)

    print('Running the Crude Screen check: the No.13 criteria check...')
    duplicate_filename_list_13=N02_13_DNA_check_PCA_90_allinfo(DNA_series,filename_info)

    print('Running the Crude Screen check: the No.14 criteria check...')
    duplicate_filename_list_14=N02_14_DNA_check_PCA_95_allinfo(DNA_series,filename_info)

    #combine all possible duplicate list from No.1 checks to No.14 checks
    All_possible_duplicate_list=duplicate_filename_list_1+duplicate_filename_list_2+duplicate_filename_list_3+duplicate_filename_list_4+duplicate_filename_list_5+duplicate_filename_list_6+duplicate_filename_list_7+duplicate_filename_list_8+duplicate_filename_list_9+duplicate_filename_list_10+duplicate_filename_list_11+duplicate_filename_list_12+duplicate_filename_list_13+duplicate_filename_list_14

    # 'unique' the pair
    All_possible_duplicate_list = math_util_functions.process_file_pairs_generic(All_possible_duplicate_list)

    print('The number of the possible duplicates pairs are:')
    for i in All_possible_duplicate_list:
        print(i)

    print('The number of the possible duplicates pairs are: '+str(len(All_possible_duplicate_list)))

    return All_possible_duplicate_list


def save_txt_duplicate_list(All_possible_duplicate_list,DNA_summary_filename):
    ### write the possible_duplicate_list to text file
    [script_directory,_]=os.path.split(DNA_summary_filename)
    if not os.path.exists(script_directory):
        os.makedirs(script_directory)

    output_file=os.path.join(script_directory,'sorted_unique_pairs_generic.txt')
    with open(output_file, 'w') as file:
        for pair in All_possible_duplicate_list:
            file.write(f"{pair[0]} {pair[1]}\n")

    #output
    print('The possible duplicates pair list is stored in: '+output_file)
    print('Then, please run the M01/M02 files to determine whether the potential duplicate pairs are exact/possible/no duplicates or not')



if __name__ == '__main__':
    # filepath = '../Input_files/DNA_summary.npz'
    DNA_summary_filename = input("Please enter the path to your DNA summary files (*.npz): ").lstrip().rstrip()
    
    if validate_file(DNA_summary_filename):
        All_possible_duplicate_list=main(DNA_summary_filename)
        save_txt_duplicate_list(All_possible_duplicate_list,DNA_summary_filename)
        print('SUCCESSFULLY run the crude screen check!!')
    else:
        print("The entered path is not valid. Please ensure the path is correct and try again.")


