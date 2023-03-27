# duplicated_checking_IQuOD
The duplicate checking algorithms contributed to IQuOD group

Running orders:
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
(17) N04_check_nc_duplicate_list4.py
(18) N05_make_stat_4.m

Here, N02***.m is to calculate the 'DNA' for each profiles in different weights.

Input test file: WOD_1995_test.xlsx.

Output potential duplicated list file: DuplicateList_potential_duplicate_ALL_1995.txt


If updates, please write down the update logs here (如有更新，请在此简要说明)

(1) January 15th, 2023: updated the N04 program with adding minor revision.
(2) February 3rd,2023: expanded the N02 series of procedures, at present, the N02_1** to N02_6** programs are based on the normalization of data by row; the N02_7** to N02_12** programs are based on the normalization of data by column; the N02_13** and N02_14** program are based on the principal component analysis method.
(3) March 27th,2023:updated the N04 program with adding minor revision;Added only output duplicate data file name and accession number program to facilitate sensitivity check; Added a program to output non-duplicate data for manual inspection;Added procedures for checking sensitivity.
