# Duplicated_checking_IQuOD (DC_OCEAN)

Release v1.1

Author: Zhetao Tan (IAP/CAS), Xinyi Song (IAP/CAS), Lijing Cheng (IAP/CAS)

Contributors: International Quality-controlled Ocean Database (*IQuOD*) members

Institute of Atmospheric Physics, Chinese Academy of Sciences (IAP/CAS)

 

## Important Updates:

Now the 'DC_OCEAN' package is uploaded to pypi (https://pypi.org/project/DC-OCEAN/1.1). For those of you interests, you can easily and freely accessed via 'pip' or 'conda' with the following commands: 

```shell
pip install DC_OCEAN
```

Please make sure  **PIP** fits your version of >=python3.7



## 1. Overview

**This algorithm, namely DC_OCEAN, aims at detecting the ocean *in-situ* duplicate profiles by reducing the computational intensity in a cost-effective way.**

It utilizes a 'DNA' method, which assigns a 'DNA' to each profile using primary metadata (e.g., latitude, longitude, instrument types) and secondary data (e.g., sum of depth, sum of temperature, standard deviation of temperature in the profile). This approach is inspired by the similarity between profile data characteristics and the structure of DNA in biology, where each profile represents a 'DNA' and different metadata information corresponds to distinct segments on this 'DNA.' 

The core assumption of this checking method is that if it's a duplicate pair, most of the metadata and observational data will be identical.

The duplicate checking algorithm is contributed to the IQuOD group.

The codes need to be run with MATLAB and Python.

A scientific paper introducing this algorithm can be referenced from Song et al., 2024, FMS.

In this manual, Section 2 presents the concept and composition of DC_OCEAN, while Section 4 explains how to install DC_OCEAN in Python and get started with it. 

## 2. Introduction of DC_OCEAN

DC_OCEAN is an open-source Python library designed for detecting duplicate profiles in ocean *in-situ* observations, such as temperature and salinity profiles. It was developed to identify duplicate ocean profiles, label them, and simultaneously reduce computational demands and human workload associated with manual quality control.

> #### Why DC_OCEAN

* **DC_OCEAN** represents the first open-source software package for checking duplicate ocean observation profiles.

* The performance and robustness of **DC_OCEAN** have been meticulously analyzed and evaluated in a scientific peer-reviewed journal (refer to Song et al., 2023; in preparation).
* As part of the contribution to the International Quality-controlled Ocean Database (IQuOD ) Task Team, specifically the Duplicate Checking Task Team, the **DC_OCEAN** has been adopted and recommended by IQuOD.
* **DC_OCEAN** utilizes 'DNA' algorithms, which employ mathematical, statistical methods like the entropy weight method and principal component analysis to establish a unique 'DNA' for each profile. This approach offers greater flexibility in comparisons and significantly reduces the time complexity of the screening process, all while ensuring screening accuracy.

#### 2.1 Composition for DC_OECAN

The DC_OCEAN is composed of two main components:

* **The first component** involves the preprocessing of metadata by calculating their corresponding 'DNA' for each profile. These files are stored in the 'support' folder.

* **The second component** is the core program of DC_OCEAN, designed to determine whether potential duplicate pairs are real duplicates or not.

**The first component includes a total of 17 scripts:**

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

**The second component consists of two files:**

(1) M01_MAIN_check_nc_duplicate_manual.py

(2) M02_MAIN_check_nc_duplicate_list.py



In short, there are 6 steps to run the DC_OCEAN (see Table 1).

<center> Table 1. The composition of DC_OCEAN <center/>

| **Order** |                  **Filename**                   |                         **Comments**                         |
| :-------: | :---------------------------------------------: | :----------------------------------------------------------: |
|     1     |           support/N00_read_metadata.m           | Read the metadata and records from the WOD18 netCDF file. (Some example files are provided) |
|     2     |        support/N01_formatted_metadata.m         | Formatted the metadata and records from the data output from N00_read_metadata.m |
|     3     |                support/N02***.m                 | Utilize fourteen distinct screening criteria to calculate the 'DNA' and identify potential duplicate pairs (14 files in total order from N02_1 to N02_14). |
|     4     | support/N03_potential_duplicate_unique_list_3.m |          Combine and unique the result of N02***.m.          |
|     5     |      M01_MAIN_check_nc_duplicate_manual.py      | Determine  whether the potential duplicates from N03 is the real  duplicates or not by manually check. |
|     6     |       M02_MAIN_check_nc_duplicate_list.py       |       The same as Order 4, but in automatically check.       |

 For more deatils and interpreation of the codes above, please refer to Song et al., 2024 (FMS, under review) [1].

 

## 3. Installation

#### 3.1 Requirement packages

* Python 3 (>=3.7)
* numpy (>= 1.19.1)
* timezonefinder (>= 5.2.0)
* netCDF4 (>= 1.5.5.1)

> Computer Memory: >8GB is obligatory.
>
> MATLAB is needed to run the support codes.

#### 3.2 Installing DC_OCEAN

**Step1: Using pip to quickly install**

If you don’t already have **PIP** running on your machine, first you need to **install pip**, then you can run:

```shell
pip install DC_OCEAN
```

Please make sure  **PIP** fits your version of python3.X. In some machines, you should use pip3 to install DC_OCEAN because **"pip"** may be linked to python2.X 

If you use the [conda](http://conda.io/) package manager `conda install -c conda-forge DC_OCEAN`

Then, you will wait for **several seconds** to install the package.

 

**Step2: Making a first and easiest duplicate checking test.** 

Here, we provide a **demo**. Now, you can make a first and most effortless test to check whether the DC_OCEAN package works well.

Launch your Python 3

Then go to the `tests` folder, run the Example file: `Example1_check_nc_duplicate_demo.py`

```shell
cd <DC_OCEAN path>/tests
python3 Example1_check_nc_duplicate_demo.py
```

Then, the following information is output:

```shell
---------Please input two netCDF files which are potential duplicates--------
The first netCDF file name is: 
```

Then, input the following two NetCDF filename:

```shell
WOD18_19750602_00177_OSD.nc
WOD18_19750603_00105_OSD.nc
```

```shell
Output profile information or not(1: Yes; 0: No) : 1
```

If return the following result, congratulations!! The DC_OCEAN package works well. 

```shell
              WOD_id:             8891325 ,              8891307
            Acess_no:                 416 ,                  416
             Dataset: bottle/rossette/net ,  bottle/rossette/net
                 Lat:             39.5000 ,              39.0167
                Long:            136.7667 ,             136.7333
                Year:                1975 ,                 1975
               Month:                   6 ,                    6
                 Day:                   3 ,                    2
                Hour:                   0 ,                   20
              Minute:                  30 ,                   54
           sum_depth:            935.0000 ,             935.0000
            sum_temp:             74.5400 ,              74.5400
        sum_salinity:            341.0000 ,             341.0000
          Probe type:                     ,                     
            Recorder:                     ,                     
        Depth_number:                  10 ,                   10
       Maximum Depth:             300.000 ,              300.000
             hasTemp:                   1 ,                    1
         hasSalinity:                   1 ,                    1
           hasOxygen:                   0 ,                    0
      hasChlonophyII:                   0 ,                    0
             Country:                  28 ,                   28
            GMT_time:               0.500 ,               20.920
     Database_origin:                     ,                     
        Project_name:   POD (1963 - 1996) ,    POD (1963 - 1996)
            Platform:       TATEYAMA-MARU ,        TATEYAMA-MARU
             Vehicle:                     ,                     
           Institute:TOYAMA PREFECTURAL FISHERIES EXPERIMENTAL STATION , TOYAMA PREFECTURAL FISHERIES EXPERIMENTAL STATION
WOD_cruise_identifier:            JP024640 ,             JP024640
      Wind_Direction:55 DEGREES - 64 DEGREES ,                     
          Wind_Speed:            999.0000 ,             999.0000
           Std_depth:             91.9796 ,              91.9796
            Std_temp:              5.4069 ,               5.4069
        Std_salinity:              0.0647 ,               0.0647
     Corr_temp&depth:             -0.8597 ,              -0.8597
      Corr_sal&depth:             -0.7677 ,              -0.7677


Spatial-temporal checks--Simultaneously but at different location: 0
Spatial-temporal checks--Simultaneously and co-located: 0
Correlation check: 0
Truncation check: 0
Layer by layer check - wrong location: 1
Layer by layer check - wrong date: 1
Layer by layer check - wrong time: 1
Layer by layer check - wrong country: 0
Layer by layer check - wrong instrument: 0
Exact duplicates check: 0
Interpolation (missing data) check: 0
CTD double data check: 0

    depth1     depth2 depth_diff      temp1      temp2  temp_diff   sal_diff
     0.000      0.000      0.000    17.1000    17.1000     0.0000     0.0000
    10.000     10.000      0.000    14.1300    14.1300     0.0000     0.0000
    20.000     20.000      0.000    12.6000    12.6000     0.0000     0.0000
    30.000     30.000      0.000     9.6700     9.6700     0.0000     0.0000
    50.000     50.000      0.000     6.9900     6.9900     0.0000     0.0000
    75.000     75.000      0.000     5.3100     5.3100     0.0000     0.0000
   100.000    100.000      0.000     4.2600     4.2600     0.0000     0.0000
   150.000    150.000      0.000     2.5000     2.5000     0.0000     0.0000
   200.000    200.000      0.000     1.3200     1.3200     0.0000     0.0000
   300.000    300.000      0.000     0.6600     0.6600     0.0000     0.0000
Duplicate result is: Near Duplicate
```

Now, you can get started with DC_OCEAN!



## 4. Getting Started with DC_OCEAN

Here, we will use the *in-situ* observational profiles in 1995 downloaded from the World Ocean Database (WOD18) to run the DC_OCEAN, aiming to detect the potential duplicate profiles within this dataset. 

#### 4.1 Run support files 

Begin by using `CODCv1_1995.xlsx` as the input file to execute `support/N01_read_excel_1.m`. This script preprocesses the profile data and metadata, employing ASCII to transform character (string) variables into numerical variables. This process generates the `DNA_summary_1995.mat` file, saved in the `Input_files` folder. You'll find three variables in this file: `DNA_series`, `filename_info`, and `variable_name`.

Next, use the `DNA_summary_1995.mat` file as input for the sequential execution of `support/N02_1_duplicate_check_main_2_mapstd_average.m` through `N02_14_duplicate_check_main_2_PAC_95_allinfo.m`. These scripts apply various 'DNA' algorithms to uncover as many potential duplicates as possible. Each `N02` script generates a corresponding `.txt` file stored in the year-specific folder within `support/potential_duplicates_output`.

Finally, execute `support/N03_potential_duplicate_unique_list.m` to combine and unique the result of the `N02***.m` . The mixed results are saved in the file `support/potential_duplicate_ALL_1995_unique_1119.txt`, which can be easily opened using Excel for further examination.

#### 4.2 Run main files

This program aims to use the knowledge of physical oceanography and the expert experiences in the field tests to determine the authenticity of potential duplicates identified in Section 4.1.  A total of 7 criteria are set.

* Simultaneously but at a different location

* Simultaneously and co-located

* Correlation check

* Truncation check

*  Layer-by-layer check (wrong location、wrong date、wrong time、wrong country、wrong instrument)

* Exact duplicates check (i.e., depth-by-depth check)

* Interpolation (missing data) check

**Two input and output methods are provided:**

1. Manually inputting the file name of potential duplicate pairs using `M01_MAIN_check_nc_duplicate_manual.py`, which checks one pair at a time (refer to 4.2.1 for details).
2. Directly using the `.txt` file output in Section 4.1 as input for the main program (`M02_MAIN_check_nc_duplicate_list.py`) to automatically check all potential duplicate data pairs (see Section 4.2.2 for details).

 ##### 4.2.1 Manual check: M01_MAIN_check_nc_duplicate_manual.py

Using `M01_MAIN_check_nc_duplicate_manual.py` enables a manual check, providing a side-by-side comparison of metadata information between potential duplicate and unduplicated profile data pairs. This facilitates a more precise determination of their duplicate status.

The `M01_MAIN_check_nc_duplicate_manual.py` is storage at the <DC_OCEAN> main folder.

```python
########## check whether the potential duplicates are 'true' duplicates based on some criterias
######## This check is manually, one by one pair
#### The automatic check is in the N04_check_nc_duplicate_list.py
##### input data: the potential pair
#### output: whether it is real duplicated or not. (Screen output)

import DC_OCEAN
import netCDF4 as nc
import numpy as np
import math
import os
from DC_OCEAN.util import country_table as t_country
from DC_OCEAN.util import compair_main as compair_main
import warnings
warnings.filterwarnings('ignore')

Class Duplicate_check(object):
	def __init__(self):
    pass
    
    def run(self):
        while True:
            print('---------Please input two netCDF files which are potential duplicates--------')
            file1=input('The first netCDF file name is: ').rstrip().lstrip()
            file2=input('The second netCDF file name is: ').rstrip().lstrip()
            isOutput_detail = input("Output profile information or not(1: Yes; 0: No)")
            ......
            ### Read the first netCDF file data
            content1=self.read_nc_data(filepath1)   # content1 is a dictionary
            ### Read the second netCDF file data
            content2=self.read_nc_data(filepath2)

            ### Output the information of two netCDF files
            self.output_info_pairs(content1,content2)

            ### Determine whether it is really repeated
            isDuplicated=compair_main.compair(content1,content2)

            if(isOutput_detail=='1'):
                self.output_detail(content1,content2)

            if(isDuplicated==1):
                print('Duplicate result is: Exact Duplicate')
            elif(isDuplicated==2):
                print('Duplicate result is: Near Duplicate')
            else:
                print('Duplicate result is: Not Duplicate')
                
	def output_detail(self,content1,content2):
  	# Output data information and secondary processing information of two profiles
  	......
    def output_info_pairs(self,content1,content2):
  	# Output detail metadata information of two profiles
  	......
  	def read_nc_data(self,file):
  	# Read netCDF file
  	......
  	def find_id_country(self,country_name):
  	......
  	def add_parameters(self,params, **kwargs):
  	......
  	def find_order_dataset(self,dataset_name):
  	......
    
def main(netCDF_filepath):
    dc=Duplicate_check()
    dc.run(netCDF_filepath)

if __name__ == '__main__':
    ### You should revise this path to fit your own purposes
    netCDF_filepath='/Users/zqtzt/Downloads/WODdata'  
    main(netCDF_filepath)
```

> Please update the ***netCDF_filepath*** to suit your specific case. We've provided a demo using WOD18 data for all of 1995 in netCDF format. You can download the compressed file [here](www.ocean.iap.ac.cn/) and then extract it to your local directory.

Run the file  `M01_MAIN_check_nc_duplicate_manual.py` and enter the file name of the profile pair to be checked according to the prompts on the screen output:

```shell
---------Please input two netCDF files which are potential duplicates--------
The first netCDF file name is: CASv1_T_S_19950308_00729_CTD.nc
The second netCDF file name is: CASv1_T_S_19950308_00730_CTD.nc
```

Input 1 or 0 to determine whether to output metadata information according to your needs (yes is 1, no is 0); in this case, 1 is used.

```shell
Output profile information or not(1: Yes; 0: No)1
```

```shell
              WOD_id:            13363774 ,             13363775
            Acess_no:               68956 ,                68956
             Dataset:                 CTD ,                  CTD
      Cast_Direction:                     ,                     
                 Lat:             38.9167 ,              38.9317
                Long:            141.7333 ,             141.7400
                Year:                1995 ,                 1995
               Month:                   3 ,                    3
                 Day:                   8 ,                    8
                Hour:                   1 ,                    1
              Minute:                  27 ,                   27
           sum_depth:            110.0000 ,             110.0000
            sum_temp:             48.2300 ,              37.9420
        sum_salinity:            170.8860 ,             136.7540
          Probe type:                     ,                     
            Recorder:                     ,                     
        Depth_number:                   5 ,                    4
       Maximum Depth:              50.000 ,               50.000
             hasTemp:                   1 ,                    1
         hasSalinity:                   1 ,                    1
           hasOxygen:                   0 ,                    0
      hasChlonophyII:                   0 ,                    0
             Country:                  28 ,                   28
            GMT_time:               1.450 ,                1.450
     Database_origin:                     ,                     
        Project_name:   POD (1963 - 1996) ,    POD (1963 - 1996)
            Platform:          IWATE MARU ,           IWATE MARU
             Vehicle:                     ,                     
           Institute:IWATE PREFECTURAL FISHERIES EXPERIMENTAL STATION , IWATE PREFECTURAL FISHERIES EXPERIMENTAL STATION
WOD_cruise_identifier:            JP035710 ,             JP035710
      Wind_Direction:                     ,                     
          Wind_Speed:            999.0000 ,             999.0000
           Std_depth:             17.2047 ,              14.7902
            Std_temp:              0.3426 ,               0.1142
        Std_salinity:              0.0406 ,               0.0355
     Corr_temp&depth:             -0.4092 ,               0.8784
      Corr_sal&depth:              0.9490 ,               0.9125


Spatial-temporal checks--Simultaneously but at different location: 0
Spatial-temporal checks--Simultaneously and co-located: 0
Correlation check: 0
Truncation check: 0
Layer by layer check - wrong location: 0
Layer by layer check - wrong date: 0
Layer by layer check - wrong time: 0
Layer by layer check - wrong country: 0
Layer by layer check - wrong instrument: 0
Exact duplicates check: 0
Interpolation (missing data) check: 1
CTD double data check: 0


    depth1      temp1  salinity1      depth2      temp2  salinity2
     0.000    10.3000    34.1270      10.000     9.3250    34.1390
    10.000     9.3200    34.1400      20.000     9.4300    34.1710
    20.000     9.4300    34.1720      30.000     9.5930    34.2180
    30.000     9.5900    34.2190      50.000     9.5940    34.2260
    50.000     9.5900    34.2280                                  
Duplicate result is: Near Duplicate
```

According to the running results, the two profile data are near duplicate, with the duplication type attributed to interpolation (missing data).

##### 4.2.2 M02_MAIN_check_nc_duplicate_list.py

The theory is consistent with Section 4.2.1, with the only difference being the modification of input and output formats.

It should be noted that **the input of this code is sourced from the output in 4.1**

```python
###### Only the parts that are different from the code in 4.2.1 are shown here.
#### This program is used to determine whether the potential duplicate pairs quickly identified in the N02 step are actually duplicated, and if so, output
##### input data: the txt file output from the N03_potential_duplicate_unique_list_3.m
#### output: two txt files: the duplicated list and the non-duplicated list. These two files can be opened by using Excel
import DC_OCEAN
import netCDF4 as nc
import numpy as np
import math
import os
from DC_OCEAN.util import country_table as t_country
from DC_OCEAN.util import compair_main as compair_main
import warnings
warnings.filterwarnings('ignore')
warnings

class Duplicate_check(object):
        def __init__(self):
        pass

    def read_potential_txt(self,txt_path):
        data=[]
        with open(txt_path,'r') as f:
            for line in f.readlines():
                ss=line.split()
                data.append(ss)

        return data

    def run(self,netCDF_filepath,potential_txt_path):

        ### Read potential_files_txt
        potential_files_list=self.read_potential_txt(potential_txt_path)

        potential_output_path='DuplicateList_'+potential_txt_path
        duplicate_number=0
        fid_duplicate_list=open(potential_output_path,'w+')
        print('filename1, filename2, unique_id_cast1, unique_id_cast2, same_moment_diff_loc_cruise, diff_records_in_same_Moment&Loc_cruise, scaled_records, rounded_truncate, wrong_location, wrong_date, wrong_moments, wrong_country, wrong_instru_types, identical_info, interpolated_pairs, CTD multiple observations, ',end='',file=fid_duplicate_list)
        print('Instrument_cast1, Instrument_cast2, Accession_cast1, Accession_cast2, lat_cast1, lat_cast2, lon_cast1, lon_cast2, year_cast1, year_cast2, month_cast1, month_cast2, day_cast1, day_cast2, hour_cast1, hour_cast2, minute_cast1, minute_cast2,',end='',file=fid_duplicate_list)
        print('probe_type_cast1, probe_type_cast2, recorder_cast1, recorder_cast2, depth_number_cast1, depth_number_cast2, maximum_depth_cast1, maximum_depth_cast2, country_cast1, country_cast2, GMT_time_cast1, GMT_time_cast2, dbase_orig_cast1, dbase_orig_cast2,',end='',file=fid_duplicate_list)
        print('project_cast1, project_cast2, Platform_cast1, Platform_cast2, ocean_vehicle_cast1, ocean_vehicle_cast2,WOD_cruise_identifier1,WOD_cruise_identifier2,Institute1,Institute2,need_z_fix1,need_z_fix2,sum_depth_cast1, sum_depth_cast2, sum_temp_cast1, sum_temp_cast2, sum_salinity_cast1, sum_salinity_cast2',file=fid_duplicate_list)

        ### Output a txt file containing nonduplicated profiles
        potential_output_unduplicate_path = 'Unduplicatelist_' + potential_txt_path
        fid_unduplicate_list = open(potential_output_unduplicate_path, 'w+')

        for i,potential_pairs in enumerate(potential_files_list):
            file1=potential_pairs[0].rstrip().lstrip()
            for i in range(1,len(potential_pairs)):
                file2=potential_pairs[i].rstrip().lstrip()
                # isOutput_detail = input("Output profile information or not(1: Yes; 0: No)")
                isOutput_detail='0'
                ......
                
                ### Read the first netCDF file data
                try:
                    content1=self.read_nc_data(filepath1)   # content1 is a dictionary
                    ### Read the second netCDF file data
                    content2=self.read_nc_data(filepath2)
                except:
                    print('Failed reading: '+file1+' and '+file2)
                    continue

                ### Compare the data
                isDuplicated,duplicate_multimodels=compair_main.compair(content1,content2)

                ### Output nonduplicated profile pair information
                if (isDuplicated == False):
                    self.output_UnduplicateList_txt(fid_unduplicate_list,content1,content2,file1,file2)

                elif(isDuplicated==1 or isDuplicated==2):
                    ### Output pair information

                    print(file1,file2)

                    duplicate_number=duplicate_number+1

                    ### Output metadata information and duplicate type of duplicate profile pairs
                    self.output_DuplicateList_txt(fid_duplicate_list,content1,content2,duplicate_multimodels,file1,file2)

                    if(isOutput_detail=='1'):
                        self.output_detail(content1,content2)

                    if(isDuplicated==1):
                        print(file1+' v.s. '+file2+': Exact Duplicate')
                    elif(isDuplicated==2):
                        print(file1+' v.s. '+file2+': Near Duplicate')
                    else:
                        print(file1+' v.s. '+file2+': No Duplicate')

                del isDuplicated

        print("duplicate_number: " + str(duplicate_number))  
        print("Two files output: "+potential_output_unduplicate_path +' and '+potential_output_path)
        print("Finished!")
        ......
        
def main(netCDF_filepath):
    #potential_txt_path='./*.txt'
    dc=Duplicate_check()
    potential_txt_path=input('Please input the filename of the potential duplicated list filename (e.g., potential_duplicate_ALL_1995_unique.txt):')
    dc.run(netCDF_filepath,potential_txt_path)


if __name__ == '__main__':
    #please change this path to suit your own netCDF fiels path
    netCDF_filepath='/Users/zqtzt/Downloads/WODdata'  
    main(netCDF_filepath)
```

Run the code and input the result of `N03*.m` (`potential_duplicate_ALL_1995_unique_1119.txt`) according to the output of the prompt on the screen:

```shell
Please input the path to potential txt:potential_duplicate_ALL_1995_unique_1119.txt
```

Subsequently, two text files are generated:

* `DuplicateList_potential_duplicate_ALL_1995_unique.txt`：Contains filenames of duplicate data and their corresponding metadata.

* `Unduplicatelist_potential_duplicate_ALL_1995_unique.txt`: Contains filenames of unduplicate data and their corresponding metadata.

Table 2 presents the variables found in the `*.txt` files and their corresponding metadata.

<center> Table 2. The output metadata list in M02_MAIN_check_nc_duplicate_list.py <center/>

|        **Variable**        | **Corresponding metadata fullname** |
| :------------------------: | :---------------------------------: |
|          filename          |              filename               |
|         unique_id          |            WOD unique id            |
|      Instrument_cast       |               dataset               |
|       Accession_cast       |          accession number           |
|          lat_cast          |              latitude               |
|          lon_cast          |              longitude              |
|         year_cast          |                year                 |
|         month_cast         |                month                |
|          day_cast          |                 day                 |
|         hour_cast          |                hour                 |
|        minute_cast         |               minute                |
|      probe_type_cast       |             probe type              |
|       recorder_cast        |              recorder               |
|     depth_number_cast      |            depth number             |
|     maximum_depth_cast     |            maximum depth            |
|        country_cast        |               country               |
|       GMT_time_cast        |              GMT time               |
|      dbase_orig_cast       |            dbase origin             |
|        project_cast        |               project               |
|       Platform_cast        |              platform               |
|     ocean_vehicle_cast     |            ocean vehicle            |
| WOD_cruise_identifier_cast |        WOD cruise identifier        |
|       Institute_cast       |              institute              |
|      need_z_fix_cast       |             need z_fix              |
|       sum_depth_cast       |        sum of depth records         |
|       sum_temp_cast        |     sum of temperature records      |
|     sum_salinity_cast      |       sum of salinity records       |

The above data files can be viewed and modified using Excel.

  

## 5. References

For more information about the DC_OCEAN, please refer to the documents or links below:

DC_OCEAN Github Project: https://github.com/IQuOD/duplicated_checking_IQuOD

IQuOD project and Task Team Duplicates: https://www.iquod.org/about.html

**For more information about the DC_OCEAN (performance evaluation, scientific application)**, please refer to our scientific paper:

>  X. Song, Z. Tan, R. Locarnini, S. Simoncelli, R. Cowley, S.i Kizu, T. Boyer, F. Reseghetti, G. Castelao, V. Gouretski, L. Cheng, 2024: An open-source algorithm for identification of duplicates in ocean database. *Frontier in Marine Science*
>

## 6. License

**DC_OCEAN** is licensed under the [Apache-2.0 License]( https://github.com/IQuOD/duplicated_checking_IQuOD/blob/main/LICENSE).

## 7. Citation

Please remember to cite this study if you use DC_OCEAN for any purposes:

[1] X. Song, Z. Tan, R. Locarnini, S. Simoncelli, R. Cowley, S.i Kizu, T. Boyer, F. Reseghetti, G. Castelao, V. Gouretski, L. Cheng, 2024: An open-source algorithm for identification of duplicates in ocean database. *Frontier in Marine Science*

[2] Tan, Z. and Song, X. (2024) IQuOD/duplicated_checking_IQuOD: DC_OCEAN (DC_OCEAN). Zenodo. https://doi.org/10.5281/zenodo.10689637

## 8. Acknowledgment

This study is supported by the Strategic Priority Research Program of the Chinese Academy of Sciences (Grant no. XDB42040402). We extend our thanks to all the IQuOD members who contributed to the manual checks of potential duplicates. We are grateful for the support of the International Oceanographic Data and Information Exchange (IODE) program. Special thanks to Edward King from CSIRO for providing valuable insights and reference materials on duplicate checking codes.

## 9. Questions and feedback

We warmly welcome feedback, questions, forks, pull requests, and improvements for the DC_OCEAN project within the IQuOD community!!

If you have any questions, suggestions, or come across any bugs in the program, or if you're interested in debugging or enhancing the DC_OCEAN project, please don't hesitate to get in touch:

·    [Create an issue](https://github.com/IQuOD/duplicated_checking_IQuOD/issues) in the Github community

·    [Pull requests](https://github.com/IQuOD/duplicated_checking_IQuOD/pulls) your debugged/improved codes in the GitHub community

·    Send us an email at: **tanzhetao@mail.iap.ac.cn** or **songxinyi231@mails.ucas.ac.cn**

## 10. Update logs

* January 15, 2023: updated the `N04` program with adding minor revisions.
* February 3, 2023: expanded the `N02` series of procedures. At present, the `N02_1**` to `N02_6**` programs are based on the normalization of data by row; the `N02_7**` to `N02_12**` programs are based on the normalization of data by column; the `N02_13**` and `N02_14**` program are based on the principal component analysis method.
* March 29, 2023: updated the `N04` program with minor revision; Added only output duplicate data file name and accession number program to facilitate sensitivity check; Added a program to output non-duplicate data for manual inspection; Added procedures for checking sensitivity.
* August 22, 2023: Finalized the first version of the duplicate checking algorithm (v1.0)
* November 2023: Issued DC_OCEAN Python package (v1.0).
* March 2024:  Issued DC_OCEAN Python package (v1.1) and linked the package to the Zenodo.
