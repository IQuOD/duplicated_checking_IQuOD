# DC_OCEAN: An open-source algorithm for identification of duplicates in ocean databases

Release v1.3.3

Author: Zhetao Tan (IAP/CAS), Xinyi Song (IAP/CAS), Lijing Cheng (IAP/CAS), Rebecca Cowley (CSIRO), Huifeng Yuan (CNIC/CAS), Guilherme Castelao (SIO), Simona Simoncelli (INGV), Shoichi Kizu (Tohoku University), Ricardo Locarnini (NOAA/NCEI), Tim Boyer (NOAA/NCEI), Franco Reseghetti (INGV), Viktor Gouretski (IAP/CAS)

International Quality-controlled Ocean Database (*IQuOD*)



*Author affiliations:*

Institute of Atmospheric Physics, Chinese Academy of Sciences (IAP/CAS)

Climate Science Centre, Environment, Commonwealth Scientific and Industrial Research Organisation (CSIRO), Australia.

Computer Network Information Center (CNIC), Chinese Academy of Sciences

Scripps Institution of Oceanography (SIO), University of California, United States.

Istituto Nazionale di Geofisica e Vulcanologia (INGV), Italy.

Department of Geophysics, Tohoku University, Japan.

NOAA National Centers for Environmental Information, United States.



## 1. Overview

**This algorithm, namely  DC_OCEAN, aims at detecting the ocean *in-situ* duplicate profiles by reducing the computational intensity in a cost-effective way.**

It utilizes a 'profile summary score (PSS)' method, which assigns a numerical value to each profile using primary metadata (e.g., latitude, longitude, instrument types) and secondary data (e.g., sum of depth, sum of temperature, standard deviation of temperature in the profile). 

The core assumption of this algorithm is that if it's a duplicate pair, most of the metadata and observational data will be identical.

The duplicate checking algorithm can support various groups including IQuoD, IAP/CAS, WOD/NCEI, CODC etc.

The codes need to be run with Python 3.

**A scientific paper introducing this algorithm could be found in Song et al., 2024.**

## 2. Introduction of DC_OCEAN

DC_OCEAN is an open-source Python library designed for detecting duplicate profiles in ocean *in-situ* observations, such as temperature and salinity profiles. It was developed to identify duplicate ocean profiles, label them, and simultaneously reduce computational demands and human workload associated with manual quality control.

> #### Why DC_OCEAN

* **DC_OCEAN** represents the first open-source software package for checking duplicate ocean observation profiles.

* The performance and robustness of **DC_OCEAN** have been thoroughly analyzed and evaluated in a scientific peer-reviewed journal (refer to Song et al., 2024; FMS).
* This software also partly contributes to the International Quality-controlled Ocean Database (IQuOD ) Task Team, specifically the Duplicate Checking Task Team, the **DC_OCEAN** has been adopted and recommended by IQuOD.
* **DC_OCEAN** calculates the Profile Summary Score (PSS), which employ mathematical, statistical methods like the entropy weight method and principal component analysis (PCA) to establish a unique numerical value for each profile. This approach offers greater flexibility in comparisons and significantly reduces the time complexity of the screening process, all while ensuring screening accuracy.

#### 2.1 Composition for DC_OECAN

**The DC_OCEAN is composed of two main components:**

* **The first component** involves the processing of metadata by calculating their corresponding PSS for each profile. These files are stored in the 'support' folder.

* **The second component** is the core program of DC_OCEAN, designed to determine whether potential duplicate pairs are real duplicates or not.



**The first component includes a total of 2 scripts:**

(1) N00_Create_Profile_Summary_Score.py

This script aims at reading the metadata and secondary information from the original netCDF files (we use the WOD18 single netCDF format) and then processing the metadata to create a data-metadata full list for later estimating the Profile Summary Score (PSS).

(2) N01_Possible_Duplicate_Check.py

This script aims at utilizing 14 distinct screening criteria to calculate the Profile Summary Score (PSS) and identify potential duplicate pairs. The output is a potential duplicate pair list file (*.txt).



**The second component consists of two files:**

(1) Duplicate_Checker.py

(2) M00_Duplicate_Check_MAIN.py



In short, there are 4 steps to run the DC_OCEAN (see Table 1).

<center> Table 1. The composition of DC_OCEAN <center/>

| **Order** |             **Filename**              |                         **Comments**                         |
| :-------: | :-----------------------------------: | :----------------------------------------------------------: |
|     1     |   support/N00_Create_Profile_Summary_Score.py   |                   Preprocess the metadata                    |
|     2     | support/N01_Possible_Duplicate_Check.py | Utilize fourteen distinct screening criteria to calculate the Profile Summary Score and identify potential duplicate pairs. |
|     3     | Duplicate_Checker.py | Determine whether the potential duplicates in the results of Possible_Duplicate_check.py are the real duplicates or not by manually checking and automatic check. |
|     4     |  M00_Duplicate_Check_MAIN.py  |         The overall flow of duplicate check. **It is the entry for the duplicate check program.**         |

 For more details and interpreation of the codes above, please refer to Song et al., 2024, Frontier in Marine Science.

## 3. Installation

#### 3.1 Requirement packages

* Python 3 (>=3.8)
* numpy (= 1.19.1)
* timezonefinder (= 6.0.1)
* netCDF4 (= 1.5.5.1)
* pandas (= 1.0.3)
* scipy (=1.7.3)
* argparse (>=1.4.0)

> Computer Memory: >8GB is obligatory.
>

#### 3.2 Installing DC_OCEAN

Now the 'DC_OCEAN' package is uploaded to pypi (https://pypi.org/project/DC-OCEAN/1.3.2). For those of you interested, you can easily and freely access via 'pip' with the following steps:

**Step1: Using pip to quickly install**

If you don’t already have **PIP** running on your machine, first you need to **install pip**, then you can run:

```shell
pip install DC_OCEAN
```

Please make sure  **PIP** fits your version of python3.X. In some machines, you should use pip3 to install DC_OCEAN because **"pip"** may be linked to python2.X 

Then, you will wait for **several seconds** to install the package.

If you fail this step, you can manually install the package with the `DC_OCEAN-1.3.2-py3-none-any.whl`file:

```shell
pip install DC_OCEAN-1.3.2-py3-none-any
```



**Step 2: Make a first and easiest QC test.** 

Here, we provide a **demo**. Now, you can make a first try with this test to check whether the DC_OCEAN package works well.

Launch your Python 3

Then go to the `tests` folder, and run the Example file: `Example1_check_nc_duplicate_demo.py`

```shell
cd <DC_OCEAN path>/tests
python3 Example1_check_nc_duplicate_demo.py
```

Then, the following information outputs:

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

##4. Logical flow of DC_OCEAN

#### 4.1  support files to calculate the Profile Summary Score and potential duplicates list 

In this step, we will pre-processes the profiles and metadata, using ASCII to convert character (string) variables into numerical values. This process generates the `../Input_files/Profile_Summary_Score_list.npz` file. You'll find three variables in this npz file: `PSS_series`, `filename_info`, and `variable_name`.

```shell
cd <DC_ocean>/support
python python N00_Create_Profile_Summary_Score.py -i <input_path> -o <output_path>
```

The acceptable input parameters are as follows:

```
-i or --input: the path that storge all netCDF file. The default value is <DC_ocean>/Input_files/WOD18_sample_1995
-o or --output: the path that storge Profile_Summary_Score_list.npz. The default value is <DC_ocean>/Input_files
```

Here, we try the following command:

```shell
python N00_Create_Profile_Summary_Score.py -i <DC_ocean>/Input_files/WOD18_sample_1995 -o <DC_ocean>/Input_files
```

> **Please replace the '<DC_ocean>' to the DC_OCEAN installed path.**
>
> Due to the upload limitation in the GitHub repository, the <DC_ocean>/Input_files/WOD18_sample_1995 only contains less than 200 netCDF files. To test as much as possible netCDF files, **you can download the compressed files of the full data from** [here]([https://github.com/IQuOD/duplicated_checking_IQuOD/blob/main/WOD18_sample_1995.zip), and uncompress it to the `<DC_ocean>/Input_files/WOD18_sample_1995` path

If return the following result, congratulations!! The first step works well.

```
Processing file 1/1883: wod_007274912O.nc
....
....
Processing file 1880/1883: wod_007275383O.nc
Processing file 1881/1883: wod_007274923O.nc
Processing file 1882/1883: wod_007274937O.nc
Processing file 1883/1883: wod_007275397O.nc
******************************************************

The Profile Summary Score formatted file are output to current folder: <DC_ocean>\Input_files


The Profile Summary Score filename is: <DC_ocean>\Input_files\Profile_Summary_Score_list.npz
Profile_Summary_Score_list.npz Complete !
```

If a file with the same name already exists in the output path, the following question will go out, please give an answer.
```
Update Profile Summary Score list or not(1: Yes (default); 0: No): 
```

Next, use the `Profile_Summary_Score_list.npz` file as input for the sequential execution of `support/N01_Possible_Duplicate_Check.py`. These scripts apply various strategies algorithms (with 14 strategies in total) to uncover as many potential duplicates as possible. This process generates the potential duplicate pairs list:

```shell
cd <DC_ocean>/support
python N01_Possible_Duplicate_Check.py -d <DC_ocean>/Input_files/Profile_Summary_Score_list.npz
```

If return the following result, congratulations!! The first step works well.

```
loading the Profile Summary Score list files....
Running the Crude Screen check: the No.1 criteria check...
Running the Crude Screen check: the No.2 criteria check...
Running the Crude Screen check: the No.3 criteria check...
Running the Crude Screen check: the No.4 criteria check...
Running the Crude Screen check: the No.5 criteria check...
Running the Crude Screen check: the No.6 criteria check...
Running the Crude Screen check: the No.7 criteria check...
Running the Crude Screen check: the No.8 criteria check...
Running the Crude Screen check: the No.9 criteria check...
Running the Crude Screen check: the No.10 criteria check...
Running the Crude Screen check: the No.11 criteria check...
Running the Crude Screen check: the No.12 criteria check...
Running the Crude Screen check: the No.13 criteria check...
Running the Crude Screen check: the No.14 criteria check...
The number of the potential duplicates pairs are:
('wod_007276168O.nc', 'wod_007276473O.nc')
...
...
('wod_007275019O.nc', 'wod_007275021O.nc')
('wod_007275041O.nc', 'wod_007276232O.nc')
The number of the possible duplicates pairs are: 258

*************FINISHED****************
The possible duplicates list is stored in: <DC_OCEAN>\Input_files\sorted_unique_pairs_generic.txt

Then, please run the MAIN files (M00_Duplicate_Check_MAIN.py) to determine whether the potential duplicate pairs are exact/possible/no duplicates or not
SUCCESSFULLY run the crude screen check!!
```

Here, the potential duplicate pairs list is saved as `./Input_files/sorted_unique_pairs_generic.txt`, which can be easily opened using Excel for further examination.

You can now go to the Section 4.2.

#### 4.2 Run main files (M00_Duplicate_Check_MAIN.py)

This program aims to use the knowledge of physical oceanography and the expert experiences in the field tests to determine the authenticity of potential duplicates identified in Section 4.1.  A total of 7 criteria are set.

* Simultaneously but at a different location

* Simultaneously and co-located

* Correlation check

* Truncation check

*  Layer-by-layer check (wrong location, wrong date, wrong time, wrong country, wrong instrument)

* Exact duplicates check (i.e., depth-by-depth check)

* Interpolation (missing data) check

**Two input and output methods are provided:**

1. Manually inputting the file name of potential duplicate pairs by calling `DuplicateCheckeManual` in `M00_Duplicate_Check_MAIN.py` , which checks the potential duplicates one by one  (refer to 4.2.1 for details).
2. Automatically checking all the potential duplicate list by calling `DuplicateCheckeManual` in `M00_Duplicate_Check_MAIN.py` by using the `.txt` file output in Section 4.1 as input for the main program (see Section 4.2.2 for details).

##### 4.2.1 Manual check: `DuplicateCheckeManual`

Using `DuplicateCheckeManual` in  `Duplicate_Checker.py` enables a manual check, providing a side-by-side comparison of metadata information between potential duplicate and unduplicated profile data pairs. This facilitates a more precise determination of duplicates.

The manual check codes are storage in the <DC_OCEAN> main folder (`Duplicate_Checker.py`) at Line 79-131.

```python
'''
	This program is used to determine whether the potential duplicate pairs quickly identified in the N02 step are actually duplicates, if so, output
	the data: the txt file output from the ./support/N01_Possible_Duplicate_Check.py
	output: two txt files: the duplicated list and the non-duplicated list. These two files can be opened by using Excel etc.
'''
def duplicate_checke_manual(self, netCDF_filepath):
    while True:
        print('---------Please input two netCDF files which are potential duplicates--------')
        file1=input('The first netCDF file name is: ').rstrip().lstrip()
        file2=input('The second netCDF file name is: ').rstrip().lstrip()
        isOutput_detail = input("Output profile information or not(1: Yes; 0: No)")

        # index_str=file1.rfind('_')
        # date1=file1[index_str-14:index_str-6]
        # year1=date1[0:4]
        # month1=date1[4:6]
        # path1=os.path.join(netCDF_filepath,year1,month1)

        # index_str=file2.rfind('_')
        # date2=file2[index_str-14:index_str-6]
        # year2=date2[0:4]
        # month2=date2[4:6]
        # path2=os.path.join(netCDF_filepath,year2,month2)

        filepath1=os.path.join(netCDF_filepath,file1)
        filepath2=os.path.join(netCDF_filepath,file2)

        ### Read the first netCDF file data
        content1=self.read_nc_data(filepath1)   # content1 is a dictionary
        ### Read the second netCDF file data
        content2=self.read_nc_data(filepath2)

        ### Output the information of two netCDF files
        self.output_info_pairs(content1,content2)

        ### Determine whether it is really repeated
        isDuplicated,duplicate_multimodels=compair_main.compair(content1,content2)

        if(isOutput_detail=='1'):
            self.output_detail(content1,content2)

            if(isDuplicated==1):
                print('Duplicate result is: Exact Duplicate')
                elif(isDuplicated==2):
                    print('Duplicate result is: Possible Duplicate')
                else:
                    print('Duplicate result is: Not Duplicate')
```

> Please specify the ***netCDF_filepath*** to suit your specific case. We've provided a demo using WOD18 data in 1995 with netCDF format. You can download the compressed file [here](www.ocean.iap.ac.cn/) and then extract it to your local directory.



##### 4.2.2 Atuomatically check (DuplicateCheckeList)

The logical flow is consistent with Section 4.2.1, with the only difference being the modification of input and output formats.

It should be noted that **the input of this code is sourced from the output in 4.1**

```python
"""
    This program is used to determine whether the potential duplicate pairs quickly identified in the N01 step are actually duplicated, and if so, output
    input data: the txt file output from the ./support/N01_Possible_Duplicate_Check.py
    output: two txt files: the duplicated list and the non-duplicated list. These two files can be opened by using Excel etc.
"""
def duplicate_checke_multiple(self,netCDF_filepath,potential_txt_path):

    ### Read potential_files_txt
    potential_files_list=self.read_potential_txt(potential_txt_path)


    # script_directory = os.path.dirname(potential_txt_path)
    script_directory, _filename = os.path.split(potential_txt_path)

    potential_output_path=os.path.join(script_directory,'DuplicateList_'+_filename)

    duplicate_number=0
    fid_duplicate_list=open(potential_output_path,'w+')
    print('filename1, filename2, unique_id_cast1, unique_id_cast2, same_moment_diff_loc_cruise, diff_records_in_same_Moment&Loc_cruise, scaled_records, rounded_truncate, wrong_location, wrong_date, wrong_moments, wrong_country, wrong_instru_types, identical_info, interpolated_pairs, CTD multiple observations, ',end='',file=fid_duplicate_list)
    print('Instrument_cast1, Instrument_cast2, Accession_cast1, Accession_cast2, lat_cast1, lat_cast2, lon_cast1, lon_cast2, year_cast1, year_cast2, month_cast1, month_cast2, day_cast1, day_cast2, hour_cast1, hour_cast2, minute_cast1, minute_cast2,',end='',file=fid_duplicate_list)
    print('probe_type_cast1, probe_type_cast2, recorder_cast1, recorder_cast2, depth_number_cast1, depth_number_cast2, maximum_depth_cast1, maximum_depth_cast2, country_cast1, country_cast2, GMT_time_cast1, GMT_time_cast2, dbase_orig_cast1, dbase_orig_cast2,',end='',file=fid_duplicate_list)
    print('project_cast1, project_cast2, Platform_cast1, Platform_cast2, ocean_vehicle_cast1, ocean_vehicle_cast2, WOD_cruise_identifier1,WOD_cruise_identifier2,Institute1,Institute2,need_z_fix1,need_z_fix2,sum_depth_cast1, sum_depth_cast2, sum_temp_cast1, sum_temp_cast2, sum_salinity_cast1, sum_salinity_cast2',file=fid_duplicate_list)

    ### Output a txt file containing nonduplicated profiles
    potential_output_unduplicate_path = os.path.join(script_directory,'Unduplicatelist_' + _filename)
    fid_unduplicate_list = open(potential_output_unduplicate_path, 'w+')

    for i,potential_pairs in enumerate(potential_files_list):
        file1=potential_pairs[0].rstrip().lstrip()
        for i in range(1,len(potential_pairs)):
            file2=potential_pairs[i].rstrip().lstrip()
            # isOutput_detail = input("Output profile information or not(1: Yes; 0: No)")
            isOutput_detail='0'
            
            # index_str=file1.rfind('_')
            # date1=file1[index_str-14:index_str-6]
            # year1=date1[0:4]
            # month1=date1[4:6]
            # day1=date1[6:8]
            # path1=os.path.join(netCDF_filepath,year1,month1)

            # index_str=file2.rfind('_')
            # date2=file2[index_str-14:index_str-6]
            # year2=date2[0:4]
            # month2=date2[4:6]
            # day2=date2[6:8]
            # path2=os.path.join(netCDF_filepath,year2,month2)
            
            filepath1=os.path.join(netCDF_filepath,file1)
            filepath2=os.path.join(netCDF_filepath,file2)
            #print(filepath1)

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

            ### Output non-duplicate profile pair information
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
                    print(file1+' v.s. '+file2+': Possible Duplicate')
                else:
                    print(file1+' v.s. '+file2+': No Duplicate')

            del isDuplicated
    print('\n\n')
    print('***************FINISHED********************')
    print("duplicate_number: " + str(duplicate_number))
    print('\n')
    print("Two files output: "+potential_output_unduplicate_path +' and '+potential_output_path)
    print("Finished!")
```

Subsequently, two text files are generated:

* `duplicatelist_sorted_unique_pairs_generic.txt`：Contains filenames of duplicate data and their corresponding metadata.

* `Unduplicatelist_sorted_unique_pairs_generic.txt`: Contains filenames of non-duplicate data and their corresponding metadata.

Table 2 presents the variables saved in the `*.txt` files and their corresponding metadata.

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



## 5. Getting Started with DC_OCEAN

Here, we will use some *in-situ* observational profiles in 1995 downloaded from the World Ocean Database (WOD18) to run the DC_OCEAN, aiming to detect the potential duplicate profiles within this dataset. These netCDF files are stored in `<DC_ocean>/Input_files/WOD18_sample_1995`

**Please read the following instructions before you run the demo**:

Ensure all necessary dependencies are installed. You can run the MAIN files `M00_Duplicate_Check_MAIN.py` using the following command:

```shell
python M00_Duplicate_Check_MAIN.py [-i <input_directory>] [-o <output_directory>] [-d <PSS_directory>] [-m <mode>]
```

**Parameter Description:**

- `-i, --input`: Path to the input directory (default: `Input_files/WOD18_sample_1995`).
- `-o, --output`: Path to the output directory (default: `Input_files`).
- `-d, --PSSfolder`: Path to the Profile Summary Score directory (default: `Input_files`).
- `-m, --mode`: Execuation mode (default: 0)
  - `0`: Automatically duplicate checking mode (List)
  - `1`: Manual duplicate checking mode

**Example #1:** To run the code with the default parameters:

```
cd <DC_OCEAN>
python M00_Duplicate_Check_MAIN.py
Update Profile Summary Score list or not(1: Yes (default); 0: No): 1
```

Then, it will automatically running all supporting files (`N00_Create_Profile_Summary_Score.py` and `N01_Possible_Duplicate_Check.py`). The output files are the duplicate and non-duplicate lists (test files) saved in the default output directory.



**Example #2:** To customize mode parameters with different input and output path:

```shell
python M00_Duplicate_Check_MAIN.py -m 0
```

or perform manual duplicate check : 

```shell
python M00_Duplicate_Check_MAIN.py -m 1
```

within this mode, you will be required to input the file names:

```
---------Please input two netCDF files which are potential duplicates--------
The first netCDF file name is:wod_007274512O.nc
The second netCDF file name is: wod_007274809O.nc
```

choose whether you want to output the comparasion information (1-yes; 0: no):

```
Output profile information or not(1: Yes; 0: No):  1
```

Then, it will output the results similar as Section 3.



**Example #3:** One can customize parameters as needed:

```
python M00_Duplicate_Check_MAIN.py -i /path/to/input -o /path/to/output -m 0
```

In this example, the input directory and output directory **are entered via argument passing.**




## 6. Notes for WOD18 netCDF format

In this algorithm, the input data and the data format **is** WOD18 (World Ocean Database 2018) single netCDF file format. The format can be referenced [here](https://www.ncei.noaa.gov/access/world-ocean-database/wod-codes.html#second). The variables we used are shown in Table 3. 

Here, we also provide a *.cdl file [here](https://github.com/IQuOD/duplicated_checking_IQuOD/blob/main/ocean_data_netCDF_format.cdl) for user to build an input file. <u>**Therefore, if you need to use your custom format rather than using WOD18 format, please follow the *cdl file to customize your input netCDF files, otherwise the program will report errors.**</u>

> Note: If a specific data field does not have information (i.e., the data on your side doesn't contain the information/value of this specific data field), please set it to its default value.

<center>Table 3. The input WOD18 data format list for ./support/N00_read_data_metadata.py <center/>

|     Variable name      |                           Comment                            | Data type |
| :--------------------: | :----------------------------------------------------------: | :-------: |
|       Access_no        |  NODC accession number (used to find original data at NODC)  |    int    |
|        country         |                           country                            |   char    |
|        dataset         |                         WOD dataset                          |   char    |
|          lat           |                           latitude                           |   float   |
|          lon           |                          longitude                           |   float   |
|        Project         |                         Project name                         |   char    |
|      Temperature       |                         Temperature                          |   float   |
|          time          |                             time                             |  double   |
| WOD_cruise_identifier  | two byte country code + WOD cruise number (unique to country code) |   char    |
|    wod_unique_cast     |                       wod unique cast                        |    int    |
|           z            |                    depth below sea level                     |   float   |
|        Salinity        |                           Salinity                           |   float   |
|         Oxygen         |                            Oxygen                            |   float   |
|      Chlorophyll       |                         Chlorophyll                          |   float   |
| Temperature_Instrument |                 Device used for measurement                  |   char    |
|       need_z_fix       |                instruction for fixing depths                 |   char    |
|        Recorder        |              Device which recorded measurement               |   char    |
|        GMT_time        |                           GMT time                           |   float   |
|         WMO_ID         |                   WMO identification code                    |    int    |
|       dbase_orig       |           Database from which data were extracted            |   char    |
|        platform        |     Name of platform from which measurements were taken      |   char    |
|     Ocean_Vehicle      |                        Ocean vehicle                         |   char    |
|       Institute        |            name of institute which collected data            |   char    |



The *.cdl file is attached here:

```shell
netcdf sample {
dimensions:
    depth = unlimited;   // Or specify an upper limit if known

variables:
    float z(depth);              // Depth below sea level in meters
        z:units = "meters";
        z:long_name = "depth_below_sea_level";
        z:valid_range = 0.0f, 11000.0f;
        z:positive = "down"	
        z:_FillValue = NaN;

    float Temperature(depth);    // Temperature in degrees Celsius
        Temperature:units = "degrees_Celsius";
        Temperature:long_name = "sea_water_temperature";
        Temperature:valid_range = -2.0f, 50.0f;
        Temperature:_FillValue = NaN;

    float Salinity(depth);       // Salinity in PSU (Practical Salinity Units)
        Salinity:units = "PSU";
        Salinity:long_name = "sea_water_salinity";
        Salinity:valid_range = 0.0f, 45.0f;
        Salinity:_FillValue = NaN;

    float Oxygen(depth);                // Oxygen concentration in mL/L
        Oxygen:units = "Micromole per kilogram (µmol/kg)";
        Oxygen:long_name = "Sea water dissolved oxygen";
        Oxygen:valid_range = 0.0f, 1000.0f;
        Oxygen:_FillValue = NaN;

    float Chlorophyll(depth);           // Chlorophyll concentration in mg/m^3
        Chlorophyll:units = "Microgram per liter (µg l-1)";
        Chlorophyll:long_name = "Chlorophyll";
        Chlorophyll:valid_range = 0.0f, 50.0f;
        Chlorophyll:_FillValue = NaN;

    int wod_unique_cast;         // WOD unique cast identifier
        wod_unique_cast:long_name = "WOD Unique Cast Identifier";
        wod_unique_cast:valid_range = unlimited; 
        wod_unique_cast:_FillValue = -9999;

    float lat;                   // Latitude in decimal degrees
        lat:units = "degrees_north";
        lat:long_name = "Latitude";
        lat:valid_range = -90.0f, 90.0f;
        lat:_FillValue =NaN;

    float lon;                   // Longitude in decimal degrees
        lon:units = "degrees_east";
        lon:long_name = "Longitude";
        lon:valid_range = -180.0f, 180.0f;
        lon:_FillValue = NaN;

    double time;                 // Time of data collection in seconds since 1970-01-01 00:00:00 UTC
        time:units = days since 1770-01-01 00:00:00 UTC";
        time:long_name = "Time";
        time:valid_range = unlimited; 
        time:_FillValue = NaN;

    string country;              // Country name where data was collected
        country:long_name = "Country";
        country:comment = "https://www.ncei.noaa.gov/access/world-ocean-database/CODES/country-list.html"
        country:_FillValue = "";

    string Temperature_Instrument;  // Device used for temperature measurement
        Temperature_Instrument:long_name = "Temperature Instrument";
        Temperature_Instrument:comment = "https://www.ncei.noaa.gov/access/world-ocean-database/CODES/s_29_probe_type.html"
        Temperature_Instrument:_FillValue = "";

    string need_z_fix;           // Instruction for fixing depths
        need_z_fix:long_name = "Need Z Fix for XBT bias correction";
        need_z_fix:comment = "https://www.ncei.noaa.gov/access/world-ocean-database/CODES/s_54_needs_depth_fix.html"
        need_z_fix:_FillValue = "";

    string Recorder;             // Device which recorded measurement
        Recorder:long_name = "Recorder";
        Recorder:comment = "https://www.ncei.noaa.gov/access/world-ocean-database/CODES/s_32_recorder.html"
        Recorder:_FillValue = "";

    float GMT_time;              // GMT time of data collection
        GMT_time:units = "hours";
        GMT_time:long_name = "GMT Time";
        GMT_time:valid_range = 0.0f, 24.0f;
        GMT_time:_FillValue = NaN;

    int WMO_ID;                  // WMO identification code
        WMO_ID:long_name = "WMO_identification_code";
        WMO_ID:valid_range = 0, 999999;
        WMO_ID:_FillValue = NaN;

    string dbase_orig;           // Database from which data were extracted
        dbase_orig:long_name = "Database of Origin";
        dbase_orig:_FillValue = "";

    string Project;              // Name of the project
        Project:long_name = "Project";
        Project:comment = "name or acronym of project under which data were measured"
        Project:_FillValue = "";

    string platform;             // Name of platform from which measurements were taken
        platform:long_name = "Platform";
        platform:comment = "https://www.ncei.noaa.gov/access/world-ocean-database/CODES/s_3_platform.html"
        platform:_FillValue = "";

    string Ocean_Vehicle;        // Ocean vehicle used for data collection
        Ocean_Vehicle:long_name = "Ocean Vehicle";
        Ocean_Vehicle:comment = "https://www.ncei.noaa.gov/access/world-ocean-database/CODES/s_74_ocean_vehicle.html"
        Ocean_Vehicle:_FillValue = "";

    int Access_no;               // NODC accession number (used to find original data at NODC)
        Access_no:long_name = "NODC_accession_number";
        Access_no:comment = "https://www.ncei.noaa.gov/access/world-ocean-database/CODES/s_1_accession.html"
        Access_no:_FillValue =NaN;

    string Institute;            // Name of the institute which collected data
        Institute:long_name = "Institute";
        Institute:_FillValue = "";

    string WOD_cruise_identifier;  // Two byte country code + WOD cruise number (unique to country code)
        WOD_cruise_identifier:long_name = "WOD Cruise Identifier";
        WOD_cruise_identifier:comment = "Two byte country code + WOD cruise number (unique to country code)"
        WOD_cruise_identifier:_FillValue = "";

    string dataset;              // WOD dataset identifier
        dataset:long_name = "WOD_dataset";
        dataset:comment = "https://www.ncei.noaa.gov/access/world-ocean-database/CODES/wod-datasets.html"
        dataset:_FillValue = "";

// Note: If a specific data field does not have information, please set it to its default value.

```



## 7. References

For more information about the DC_OCEAN, please refer to the documents or links below:

DC_OCEAN Github Project: https://github.com/IQuOD/duplicated_checking_IQuOD

IQuOD project and Task Team Duplicates: https://www.iquod.org/about.html

**For more information about the DC_OCEAN (performance evaluation, scientific application)**, please refer to:

> X. Song, Z. Tan, R. Locarnini, S. Simoncelli, R. Cowley, S.i Kizu, T. Boyer, F. Reseghetti, G. Castelao, V. Gouretski, L. Cheng*, 2024: DC_OCEAN: An open-source algorithm for identification of duplicates in ocean database. *Frontier in Marine Science*. [https://doi.org/10.3389/fmars.2024.1403175](https://doi.org/10.3389/fmars.2024.1403175)

## 8. License

**DC_OCEAN** is licensed under the [Apache-2.0 License]( https://github.com/IQuOD/duplicated_checking_IQuOD/blob/main/LICENSE).

## 9. Citation

Please **REMEMBER** to cite this study if you use DC_OCEAN for any purposes:

**[1]** X. Song, Z. Tan, R. Locarnini, S. Simoncelli, R. Cowley, S.i Kizu, T. Boyer, F. Reseghetti, H. Yuan, G. Castelao, V. Gouretski, L. Cheng*, 2024: DC_OCEAN: An open-source algorithm for identification of duplicates in ocean database. *Frontier in Marine Science*. [https://doi.org/10.3389/fmars.2024.1403175](https://doi.org/10.3389/fmars.2024.1403175)

**[2]** Zhetao Tan, Xinyi Song, Huifeng Yuan, & Rebecca Cowley. (2024). IQuOD/duplicated_checking_IQuOD: v1.3.1 (v1.3.1). Zenodo. https://doi.org/10.5281/zenodo.12697240

## 10. Acknowledgment

This study is supported by the Strategic Priority Research Program of the Chinese Academy of Sciences (Grant no. XDB42040402). We extend our thanks to all the IQuOD members who contributed to the manual checks of potential duplicates. We are grateful for the support of the International Oceanographic Data and Information Exchange (IODE) program. Special thanks to Edward King from CSIRO for providing valuable insights and reference materials on duplicate checking codes.

## 11. Questions and feedback

We warmly welcome feedback, questions, requests for the DC_OCEAN !!

If you have any questions, suggestions, find any bugs, or you're interested in further developing the DC_OCEAN software, please contact us:

* [Create an issue](https://github.com/IQuOD/duplicated_checking_IQuOD/issues) in the GitHub community
* [Pull requests](https://github.com/IQuOD/duplicated_checking_IQuOD/pulls) your debugged/improved codes in the GitHub community.
* Send us an email at: **tanzhetao19@mails.ucas.ac.cn** or **songxinyi231@mails.ucas.ac.cn**

## 12. Update logs

* January 15, 2023: updated the `N04` program with adding minor revisions.
* February 3, 2023: expanded the `N02` series of procedures. At present, the `N02_1**` to `N02_6**` programs are based on the normalization of data by row; the `N02_7**` to `N02_12**` programs are based on the normalization of data by column; the `N02_13**` and `N02_14**` program are based on the principal component analysis method.
* March 29, 2023: updated the `N04` program with minor revision; Added only output duplicate data file name and accession number program to facilitate sensitivity check; Added a program to output non-duplicate data for manual inspection; Added procedures for checking sensitivity.
* August 22, 2023: Finalized the first version of the duplicate checking algorithm (v1.0)
* November 2023: Issued DC_OCEAN Python package (v1.0).
* March 2024: Issued DC_OCEAN Python package (v1.1) and linked the package to the Zenodo.
* May 2024: Issued DC_OCEAN Python package (v1.2)
* July 2024: Code encapsulation and framework enhancement, issue DC_OCEAN Python package (v1.3.1)
* September 2024: Paper (Song et al., 2024) published in the Frontier in Marine Sciences; Issue the DC_OCEAN Python package (v.1.3.3)
