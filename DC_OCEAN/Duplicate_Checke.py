#!/usr/bin/env python3

########################################################################
### @copyright Copyright (C) 2024 All Rights Reserved.
### @file  Duplicate_Checke.py
### @brief 
### @version
###		Date	|	Author			|	Version		|	Description
### ------------|-------------------|---------------|---------------
### 2024-03-26	| Huifeng, Zhetao   |	1.0			|	Create
######################################################################

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

    OutputDir = os.path.dirname(os.path.abspath(__file__)) + "/Input_files"
    InputDir = OutputDir + "/WOD18_sample_1995"
    DNADir = OutputDir
    iMode = 0   # 0: List ; 1: Manual
    parser = argparse.ArgumentParser(description='Duplicate Check')
    parser.add_argument("-i", "--input", type=str, default=InputDir)
    parser.add_argument("-o", "--output", type=str, default=OutputDir)
    parser.add_argument("-d", "--dna", type=str, default=DNADir)
    parser.add_argument("-m", "--mode", type=int, default=0)
    args = parser.parse_args()
    InputDir = args.input
    OutputDir = args.output
    DNADir = args.dna
    iMode = args.mode

    oChecker = Duplicate_Checker.DuplicateChecker()
    oChecker.InitEnvironment(InputDir, OutputDir)

    if(iMode == 0):
        DuplicateCheckeList(oChecker, InputDir, OutputDir)
    elif (iMode == 1):
        DuplicateCheckeManual(oChecker, InputDir, OutputDir)
