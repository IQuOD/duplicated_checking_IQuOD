#!/usr/bin/env python3

########################################################################
### @copyright Copyright (C) 2024 All Rights Reserved.
### @file  M00_Duplicate_Check_MAIN.py
### @brief 
### @version
###		Date	|	Author			|	Version		|	Description
### ------------|-------------------|---------------|---------------
### 2024-03-26	|                   |	1.1			|	Create
### 2024-06-01	|                   |	1.2			|	Create
### 2024-07-04	|                   |	1.3			|	Create
######################################################################

"""
    Main DC_OCEAN program (See Section 5 in README.md)
    By calling the Duplicate_Checker.py, this program will be run in two modes for duplicate checking.
    The lgoical flow of this program can be obtained https://github.com/IQuOD/duplicated_checking_IQuOD/blob/main/flowchart.pdf
    
    mode = 1:  DuplicateCheckeManual
    Manually check whether the potential duplicates are exact duplicates based on some criterias
    This function is manually ehck one by one pair
    input data: Filenames of the possible duplicates
    output: whether it is exact duplicated, possible duplicate or non-duplicates. (Screen output)
    
    mode = 0:  DuplicateCheckeList
    Similar with mode = 1, but check automatically with the possible duplicate list.
    input data: the txt file output from the ./support/N01_Possible_Duplicate_Check.py
    output: two txt files: the duplicated list and the non-duplicated list. These two files can be opened by using Excel etc.
"""


import os
import argparse
import netCDF4 as nc
import numpy as np
import math
import Duplicate_Checker
from util import country_table as t_country
from util import compair_main as compair_main
import warnings
warnings.filterwarnings('ignore')
warnings

def DuplicateCheckeManual(checker, InputDir, OutputDir):
    if checker.validate_file(InputDir):
        checker.duplicate_checke_manual(InputDir)
    else:
        print("The entered path of netCDF files is not valid. Please ensure the path is correct and try again.")

def DuplicateCheckeList(checker, InputDir, OutputDir):
    # input the path with filename (*.txt) of the potential duplicated list output from N01_possible_duplicates.py
    potential_txt_path = OutputDir + "/sorted_unique_pairs_generic.txt"
    if checker.validate_file(potential_txt_path):
        netCDF_filepath = InputDir
        if checker.validate_file(netCDF_filepath):
            checker.duplicate_checke_multiple(netCDF_filepath, potential_txt_path)
        else:
            print('The entered path of netCDF files is not valid. Please try again.')
    else:
        print("The entered path of potential duplicated list is not valid. Please ensure the path is correct and try again.")

if __name__ == '__main__':

    ######### optional input parameter
    OutputDir = os.path.dirname(os.path.abspath(__file__)) + "/Input_files"
    InputDir = OutputDir + "/WOD18_sample_1995"
    PSSDir = OutputDir
    iMode = 0   # 0: List ; 1: Manual
    parser = argparse.ArgumentParser(description='Duplicate Check')
    parser.add_argument("-i", "--input", type=str, default=InputDir)
    parser.add_argument("-o", "--output", type=str, default=OutputDir)
    parser.add_argument("-d", "--PSSfolder", type=str, default=PSSDir)
    parser.add_argument("-m", "--mode", type=int, default=0)
    args = parser.parse_args()
    InputDir = args.input
    OutputDir = args.output
    PSSDir = args.PSSfolder           
    iMode = args.mode
    
    #initialization the code environment (class)
    oChecker = Duplicate_Checker.DuplicateChecker()  
    oChecker.InitEnvironment(InputDir, OutputDir)  

    if(iMode == 0):
        #mode = 0:  DuplicateCheckeList
        DuplicateCheckeList(oChecker, InputDir, OutputDir)
    elif (iMode == 1):
        # mode = 1:  DuplicateCheckeManual
        DuplicateCheckeManual(oChecker, InputDir, OutputDir)
