# duplicated_checking_IQuOD (DC_OCEAN)
Author: Zhetao Tan (IAP/CAS), Xinyi Song (IAP/CAS)

Contributor: International Quality-controlled Ocean Database (*IQuOD*) members

Zenodo software DOI: https://doi.org/10.5281/zenodo.10689637


## Overview

**This algorithm aims at detecting the duplicate ocean *in-situ* profiles by reducing the computational intensity with a cost-effective way.**

It used a so called 'DNA' method, with defining a DNA for each profiles by using the primary metadata (e.g., latitude, longitude, instrument types etc.) and the secondary data (e.g., sum of depth, sum of temperature, standard deviation of temperature in the profile).

The assumption of this checking is: if it is a duplicate pair, most of the metadata and the observational data are identical.

The duplicate checking algorithm is contributed to the IQuOD group.

The codes need to be runned with MATLAB and Python 3.

A scentific paper to introduce this algorithm is in prepration.



**A detailed user manual can be found in:** `DC_OCEAN_user_document-v1.0.md`



**Below is the simply instruction of the algorithm. Further information will be provided in the scentific paper.** 

**Running orders:**

(1) ./support/N00_read_metadata.m

(2) ./support/N01_formatted_metadata.m 

(3) ./support/N02_1_duplicate_check_main_2_mapstd_average.m

(4) ./support/N02_2_duplicate_check_main_2_onlydepthtemp_compair.m

(5) ./support/N02_3_duplicate_check_main_2_weight_meta.m

(6) ./support/N02_4_duplicate_check_main_2_weighted_allinfo.m

(7) ./support/N02_5_duplicate_check_weight_meta_noLATLON.m

(8) ./support/N02_6_duplicate_check_main_2_weight_meta_noDepthinfo.m

(9) ./support/N02_7_duplicate_check_main_2_mapminmax_average.m

(10) ./support/N02_8_duplicate_check_main_2_mapminmax_onlydepthtemp_compair.m

(11) ./support/N02_9_duplicate_check_main_2_mapminmax_weight_meta.m

(12) ./support/N02_10_duplicate_check_main_2_mapminmax_weighted_allinfo.m

(13) ./support/N02_11_duplicate_check_mapminmax_weight_meta_noLATLON.m

(14) ./support/N02_12_duplicate_check_main_2_mapminmax_weight_meta_noDepthinfo.m

(15) ./support/N02_13_duplicate_check_main_2_PAC_90_allinfo.m

(16) ./support/N02_14_duplicate_check_main_2_PAC_95_allinfo.m

(17) ./support/N03_potential_duplicate_unique_list_3.m

(18) M01_MAIN_check_nc_duplicate_manual.py

(19) M02_MAIN_check_nc_duplicate_list.py



Input: A test file folder `./input_files/WOD18_sample_1995`

Output: duplicated and non-duplicated list files`DuplicateList_potential_duplicate_ALL_1995_unique.txt` and `Unduplicatelist_potential_duplicate_ALL_1995_unique.txt`

The above two output files could be opened by using Excel.



Here, `./supports/N00_read_metadata.m` is to read the input test (netCDF format) file.

`N02***.m` and `N03 `are to calculate the 'DNA' for each profile in different weights and then detect the potential duplicated pairs.

`M02_MAIN_check_nc_duplicate_list.py` is to automatically check whether the potential duplicated pairs output from N03 is the real duplicate or not.

`M01_MAIN_check_nc_duplicate_manual.py` is to manually check whether the potential duplicated pairs output from N03 is the real duplicated or not. It needs to input the filename of the potential duplicated pairs manually.



## Update logs

(1) January 15th, 2023: updated the N04 program with adding minor revision.
(2) February 3rd,2023: expanded the N02 series of procedures, at present, the N02_1** to N02_6** programs are based on the normalization of data by row; the N02_7** to N02_12** programs are based on the normalization of data by column; the N02_13** and N02_14** program are based on the principal component analysis method.
(3) March 29th,2023:updated the N04 program with adding minor revision; Added only output duplicate data file name and accession number program to facilitate sensitivity check; Added a program to output non-duplicate data for manual inspection;Added procedures for checking sensitivity.

(4) August 22rd, 2023: Finalized the first version of the duplicate checking alogrithim (v1.0)

(5) January 5, 2024:Improved version  (v1.1)

## Citation

[1] Xinyi Song, Zhetao Tan,  Lijing Cheng et al, 2024: An open source algorithm of duplicate checking for ocean *in-situ* profiles. (In prepration)

[2] Xinyi Song, Zhetao Tan, Lijing Cheng. 2023, A benchmark dataset for ocean profiles duplicate checking. http://dx.doi.org/10.12157/IOCAS.20230821.001



