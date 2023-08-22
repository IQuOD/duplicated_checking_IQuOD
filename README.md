# duplicated_checking_IQuOD
Author: Zhetao Tan (IAP/CAS), Xinyi Song (IAP/CAS)

Contributor: International Quality-controlled Ocean Database (*IQuOD*) members



## Overview

**This algorithm aims at detecting the duplicate ocean *in-situ* profiles by reducing the computational intensity with a cost-effective way.**

It used a so called 'DNA' method, with defining a DNA for each profiles by using the primary metadata (e.g., latitude, longitude, instrument types etc.) and the secondary data (e.g., sum of depth, sum of temperature, standard deviation of temperature in the profile).

The assumption of this checking is: if it is a duplicate pair, most of the metadata and the observational data are identical.

The duplicate checking algorithm is contributed to the IQuOD group.

The codes need to be runned with MATLAB and Python 3.

A scentific paper to introduce this algorithm is in prepration.



**Below is the simply instruction of the algorithm. Further information will be provided in the scentific paper.** 

**Running orders:**

(1) N01_read_excel_1.m

(2) N02_1_duplicate_check_main_2_mapstd_average.m

(3) N02_2_duplicate_check_main_2_onlydepthtemp_compair.m

(4) N02_3_duplicate_check_main_2_weight_meta.m

(5) N02_4_duplicate_check_main_2_weighted_allinfo.m

(6) N02_5_duplicate_check_weight_meta_noLATLON.m

(7) N02_6_duplicate_check_main_2_weight_meta_noDepthinfo.m

(8) N02_7_duplicate_check_main_2_mapminmax_average.m

(9) N02_8_duplicate_check_main_2_mapminmax_onlydepthtemp_compair.m

(10) N02_9_duplicate_check_main_2_mapminmax_weight_meta.m

(11) N02_10_duplicate_check_main_2_mapminmax_weighted_allinfo.m

(12) N02_11_duplicate_check_mapminmax_weight_meta_noLATLON.m

(13) N02_12_duplicate_check_main_2_mapminmax_weight_meta_noDepthinfo.m

(14) N02_13_duplicate_check_main_2_PAC_90_allinfo.m

(15) N02_14_duplicate_check_main_2_PAC_95_allinfo.m

(16) N03_potential_duplicate_unique_list_3.m

(17) N04_check_nc_duplicate_list.py

(18) N04_check_nc_duplicate.py



Input: A test file `CODCv1_1995.xlsx.`

Output: duplicated and non-duplicated list files`DuplicateList_potential_duplicate_ALL_1995_unique.txt` and `Unduplicatelist_potential_duplicate_ALL_1995_unique.txt`

The above two output files could be opened by using Excel.



Here, `N01_read_excel_1.m` is to read the input test file.

`N02***.m` and `N03 `are to calculate the 'DNA' for each profiles in different weights and then detect the potential duplicated pairs.

`N04_check_nc_duplicate_list.py` is to automatically check whether the potential duplicated pairs output from N03 is the real duplicated or not.

`N04_check_nc_duplicate.py` is to manually check whether the potential duplicated pairs output from N03 is the real duplicated or not. It needs to input the filename of the potential duplicated pairs manually.



## Update logs

(1) January 15th, 2023: updated the N04 program with adding minor revision.
(2) February 3rd,2023: expanded the N02 series of procedures, at present, the N02_1** to N02_6** programs are based on the normalization of data by row; the N02_7** to N02_12** programs are based on the normalization of data by column; the N02_13** and N02_14** program are based on the principal component analysis method.
(3) March 29th,2023:updated the N04 program with adding minor revision;Added only output duplicate data file name and accession number program to facilitate sensitivity check; Added a program to output non-duplicate data for manual inspection;Added procedures for checking sensitivity.

(4) August 22rd, 2023: Finalized the first version of the duplicate checking alogrithim (v1.0)



## Citation

[1] Xinyi Song, Zhetao Tan,  Lijing Cheng et al, 2023: An open source algorithm of duplicate checking for ocean *in-situ* profiles. (In prepration)

[2] Xinyi Song, Zhetao Tan, Lijing Cheng. 2023, A benchmark dataset for ocean profiles duplicate checking. http://dx.doi.org/10.12157/IOCAS.20230821.001



