import numpy as np
import math
from timezonefinder import TimezoneFinder
from datetime import datetime

##################### Further improvements can be made based on known duplications! #####################
def compair(content1,content2):
    temp1=content1['temp']
    depth1=content1['depth']
    salinity1=content1['sal']

    temp2=content2['temp']
    depth2=content2['depth']
    salinity2=content2['sal']

    # Nonduplicated：0 ; Exact duplicated: 1 ; Near duplicated: 2

    isDuplicate=False

    ### 2023.1.2 added flag
    near_duplicated=False
    exact_duplicated=False

    duplicate_check1=False  #Spatial-temporal checks: Simultaneously but at different location
    duplicate_check2=False  #Spatial-temporal checks: Simultaneously and co-located
    duplicate_check3=False  #Data records checks: Correlation check
    duplicate_check4=False  #Data records checks: Truncation check
    duplicate_check5=False  #Data records checks: Layer by layer check-wrong location
    duplicate_check6=False  #Data records checks: Layer by layer check-wrong date
    duplicate_check7=False  #Data records checks: Layer by layer check-wrong time
    duplicate_check8=False  #Data records checks: Layer by layer check-wrong country
    duplicate_check9=False  #Data records checks: Layer by layer check-wrong instrument
    duplicate_check10=False #Data records checks: Accurate repeat check (depth-by-depth check)
    duplicate_check11=False #Data records checks: Interpolation check
    duplicate_check12=False #Data records checks: CTD double data check
    
    #### First check: Simultaneously but at different location
    #### Navigation observation--On the same voyage two observations were made on the same ship at the same time in different locations(unreasonable)
    if(content1['dataset_id'] in [1,2,3,4] and content2['dataset_id'] in [1,2,3,4]):  # instrument type: XBT MBT BOTTLE CTD
        index1=content1['year'] == content2['year'] and content1['month'] == content2['month'] and content1['day'] == content2['day'] and content1['hour'] == content2['hour'] and content1['minute'] == content2['minute']
        index2=content1['hour'] != 0 and content2['hour']!=0 # hour is not 0
        index3=content1['ocean_vehicle']==content2['ocean_vehicle'] and content1['Platform']==content2['Platform'] and content1['Institute']==content2['Institute']
        if(index1 and index2 and index3):
            ### Plantform code not missing but without unique identifier(MMSI; IMO; HULL; call-sign etc.)--nonduplicated
            if(content1['Platform']!=''):
                if ('(' not in content1['Platform']!=''):
                    duplicate_check1 = False
            else:  # otherwise continue to judge
                if (content1['ocean_vehicle'] == '' or content1['WOD_cruise_identifier'] == ''):
                    duplicate_check1=False
                else: # same time(accurate to minute)
                    if(content1['minute'] != 0 and content2['minute']!=0):
                        if(np.abs(content1['latitude']-content2['latitude'])>0.01 or np.abs(content1['longitude']-content2['longitude'])>0.01):  # different locations(threshold value: 0.01)
                            duplicate_check1=True  #Near duplicated
                        else:  # same time(accurate to hour)
                            juli=distance(content1['latitude'],content2['latitude'],content1['longitude'],content2['longitude'])
                            if(juli>30):   # threshold value: 30km ( With a maximum speed of 30km per hour, the maximum distance between two positions within an hour is unlikely to exceed 30km, and if it does, it is unlikely to happen, and there is reason to believe that it may be repeated.)
                               duplicate_check1=True  #Near duplicated

    #### Second check: Simultaneously and co-located
    #### (excluding buoy data) Multiple observations occur at the same time (accruate to minute) and at the same location
    if (content1['dataset_id'] not in [5,7,9] and content2['dataset_id'] not in [5,7,9]):
        index1=content1['year'] == content2['year'] and content1['month'] == content2['month'] and content1['day'] == content2['day'] and content1['hour'] == content2['hour'] and content1['minute'] == content2['minute']
        index2=(content1['hour'] != 0 and content2['hour']!=0 and content1['minute'] != 0 and content2['minute'] != 0) or (content1['hour'] == 0 and content2['hour']==0 and content1['minute'] != 0 and content2['minute'] != 0)
        index3=content1['ocean_vehicle']==content2['ocean_vehicle'] and content1['Platform']==content2['Platform']
        juli=distance(content1['latitude'],content2['latitude'],content1['longitude'],content2['longitude'])
        # distance less than 1km (Old instruments do not record latitude and longitude with as much precision and may only be retained to two decimal places, resolution is about 0.01 degrees.)
        index4=juli<=1
        index5=(content1['depth_number']>1) and (content2['depth_number']>1)
        if(index1 and index2 and index3 and index4 and index5):
            duplicate_check2=True

        ### Excluding OSD and CTD、glider、Argo、APB (High precision, as long as the time, minutes, are different, can be judged as "different time"), plantform codes not missing, same location, same time (within 1 hour)--near duplicated
        if (content1['dataset_id'] not in [1,2,6,8,11] and content2['dataset_id'] not in [1,2,6,8,11]) and (content1['Platform'] != '' and content2['Platform'] != '') and (content1['year'] == content2['year']):
            index1 = juli <= 0.01
            # Use the datetime libraries to calculate the time difference
            date1 = datetime(content1['year'], content1['month'], content1['day'], content1['hour'],content1['minute'])
            date2 = datetime(content2['year'], content2['month'], content2['day'], content2['hour'],content2['minute'])
            duration=date1-date2
            if duration.days < 0:
                duration = date2 - date1
            days=abs(duration.days)
            hours=abs(duration.seconds/3600)
            index2 = days ==0
            index3 = hours <= 1
            index4 = (content1['hour'] != 0 and content2['hour'] != 0 and content1['minute'] != 0 and content2['minute'] != 0)
            # Same ship--plantform codes has unique identifier(MMSI; IMO; HULL; call-sign etc.)
            index5 = ('(' in content1['Platform'] and '(' in content2['Platform'] and content1['Platform'] == content2['Platform'])
            index6 = content1['ocean_vehicle'] == content2['ocean_vehicle'] and content1['Platform'] == content2['Platform']
            index7 = content1['depth_number'] > 1 and content2['depth_number'] > 1
            if (index1 and index2 and index4 and index3 and index5 and index6 and index7):
                duplicate_check2 = True

        ### Add meteorological wind speed and direction information for judgment--has wind speed and direction information but not same--nonduplicated
        if((content1['Wind_Direction']!='' and content2['Wind_Direction']!='') or (content1['Wind_Speed']!=999 and content2['Wind_Speed']!=999)):
            if(np.abs(content1['Wind_Speed']-content2['Wind_Speed'])>1e-3 or content1['Wind_Direction'] != content2['Wind_Direction']):
                duplicate_check2=False

        ### The up and down data of instruments glider, CTD, Argo(PFL) are not duplicates
        if  (content1['dataset_id']  in [2,11,8] and content2['dataset_id']  in [2,11,8]):
            if (np.abs(content1['latitude'] - content2['latitude']) < 0.1 and np.abs(content1['longitude'] - content2['longitude']) < 0.1):  # Approximate latitude and longitude
                if(content1['Cast_Direction']!='' and content2['Cast_Direction']!='' and content1['Cast_Direction']!=content2['Cast_Direction']): # has Cast_Direction information; one is UP the other is DOWN
                    duplicate_check2 = False

        ### The GPS accuracy of CTD, APB, Argo and glider is very high, and if the longitude and latitude are different after three or four decimal places, it can not be determined as the same place
        if (content1['dataset_id'] in [2,6,8,11] and content2['dataset_id'] in [2,6,8,11] and duplicate_check2 == True):
            index1=(np.abs(content1['latitude'] - content2['latitude']) != 0) or (np.abs(content1['longitude'] - content2['longitude']) != 0) # Latitude and longitude are not exactly the same
            if (index1):
                duplicate_check2 = False

        ### Exclude surface data
        if(len(temp1[~np.isnan(temp1)])<2 or len(temp2[~np.isnan(temp2)])<2):
            duplicate_check2 = False

    #### Third check: Correlation check
    #### Depth or temperature or salinity are scaled or translation
    if(np.abs(content1['cor_temp_depth']-content2['cor_temp_depth'])<1e-3 or (np.abs(content1['cor_sal_depth']-content2['cor_sal_depth'])<1e-4) and content1['cor_sal_depth']!=999 and content2['cor_sal_depth']!=999 and  content1['cor_temp_depth']!=999 and content2['cor_temp_depth']!=999):
        if(content1['dataset_id'] == content2['dataset_id'] and content1['depth_number']==content2['depth_number']):
            juli_distance = distance(content1['latitude'], content2['latitude'], content1['longitude'], content2['longitude'])
            if(juli_distance<0.5):    # distance less than 500m
                ### After removing nan data from temperature data, depth number still more than 2
                if(content1['depth_number']>2 and content2['depth_number']>2 and len(temp1[~np.isnan(temp1)])>2 and len(temp2[~np.isnan(temp2)])>2):
                    isMode_depth=check_Mode(depth1,depth2,0.85)
                    isMode_temp=check_Mode(temp1,temp2,0.85)
                    if(isMode_temp or isMode_depth):
                        duplicate_check3=True

                    ### After removing nan data from salinity data, depth number still more than 2
                        if(content1['cor_sal_depth']!=999 and content2['cor_sal_depth']!=999 and len(salinity1[~np.isnan(salinity1)])>2 and len(salinity2[~np.isnan(salinity2)])>2):
                            isMode_salinity=check_Mode(salinity1,salinity2,0.85)
                            if(isMode_salinity):
                                duplicate_check3=True

        ### Exclude fixed-point observations from coastal stations
        if (content1['Platform']==content2['Platform']):
            index1 = (np.abs(content1['latitude'] - content2['latitude'] ) < 0.01 and np.abs( content1['longitude'] - content2['longitude'] ) < 0.01) # fixed-point observation
            index2 = (np.sum(temp1) != np.sum(temp2)) # temperature are not the same
            ### It is very shallow along the coast, and the temperature measured at each layer may be the same
            t_ave1 = np.sum(temp1) / content1['depth_number']
            t_ave2 = np.sum(temp2) / content2['depth_number']
            index3 = (np.abs(np.sum(temp1 - t_ave1)) < 0.1)
            index4 = (np.abs(np.sum(temp2 - t_ave2)) < 0.1)
            if (index1 and index2 and index3 and index4):
                duplicate_check3 = False

        ### The up and down data of instruments GLD are not duplicates
        if (content1['dataset_id'] == 11 and content2['dataset_id'] == 11):
            index1 = content1['year'] == content2['year'] and content1['month'] == content2['month'] and content1['day'] == content2['day'] and content1['hour'] == content2['hour'] and content1['minute'] == content2['minute']
            index2 = (content1['hour'] != 0 and content2['hour'] != 0 and content1['minute'] != 0 and content2['minute'] != 0) or (content1['hour'] == 0 and content2['hour'] == 0 and content1['minute'] != 0 and content2['minute'] != 0) or (content1['hour'] != 0 and content2['hour'] != 0 and content1['minute'] == 0 and content2['minute'] == 0)
            index3 = (content1['ocean_vehicle'] == content2['ocean_vehicle'] and content1['ocean_vehicle'] != '') or (content1['WOD_cruise_identifier'] == content2['WOD_cruise_identifier'] and content1['WOD_cruise_identifier'] != '')
            index4 = content1['access_no'] == content2['access_no']
            index5 = (np.abs(content1['latitude'] - content2['latitude']) < 0.1 and np.abs(content1['longitude'] - content2['longitude']) < 0.1)  # Approximate latitude and longitude
            index6 = (content1['Cast_Direction'] != '' and content2['Cast_Direction'] != '' and content1['Cast_Direction'] !=content2['Cast_Direction'])  # has Cast_Direction information; one is UP the other is DOWN
            if (index1 and index2 and index3 and index4 and index5 and index6):
                duplicate_check3 = False

    #### Fourth check: Truncation check
    #### Temperature or salinity records are rounded or truncated
    ### Same depth number
    ### Qualification: The depth number cannot be too small (more than 3)
    if(content1['depth_number']==content2['depth_number'] and np.abs(content1['sum_depth']-content2['sum_depth'])<1 and content1['depth_number']>3 and content2['depth_number']>3):
        level=content1['depth_number']
        if(np.all(np.abs(temp1-temp2)<1e-5)):  # exclude exact duplicates
            duplicate_check4=False
        else:
            index1=np.sum(np.abs(np.round(temp1,2)-temp2)<1e-6)/level>0.85 or np.sum(np.abs(np.round(temp2,2)-temp1)<1e-6)/level>0.85
            index2=np.sum(np.abs(np.round(temp1,1)-temp2)<1e-6)/level>0.85 or np.sum(np.abs(np.round(temp2,1)-temp1)<1e-6)/level>0.85
            index3=np.sum(np.abs(np.round(temp1,0)-temp2)<1e-6)/level>0.85 or np.sum(np.abs(np.round(temp2,0)-temp1)<1e-6)/level>0.85

            index4= np.sum(np.abs(truncate(temp1,2)-temp2)<1e-6)/level>0.85 or np.sum(np.abs(truncate(temp2,2)-temp1)<1e-6)/level>0.85
            index5= np.sum(np.abs(truncate(temp1,1)-temp2)<1e-6)/level>0.85 or np.sum(np.abs(truncate(temp2,1)-temp1)<1e-6)/level>0.85
            index6= np.sum(np.abs(truncate(temp1,0)-temp2)<1e-6)/level>0.85 or np.sum(np.abs(truncate(temp2,0)-temp1)<1e-6)/level>0.85

            if(index1 or index2 or index3 or index4 or index5 or index6):
                duplicate_check4=1

    ### Different depth number
    ### Qualification: The depth number cannot be too small (more than 3)
    elif (content1['depth_number']!=content2['depth_number'] and content1['depth_number']>3 and content2['depth_number']>3):
        if (content1['depth_number'] > content2['depth_number']):
            # depth2 interpolates to the depth of depth1
            temp2_interp = np.interp(depth1, depth2, temp2, left=np.nan, right=np.nan, period=None)
            level1=content2['depth_number']
            index1 = np.sum(np.abs(np.round(temp1, 2) - temp2_interp) < 1e-6) / level1 > 0.85 or np.sum(np.abs(np.round(temp2_interp, 2) - temp1) < 1e-6) / level1 > 0.85
            index2 = np.sum(np.abs(np.round(temp1, 1) - temp2_interp) < 1e-6) / level1 > 0.85 or np.sum(np.abs(np.round(temp2_interp, 1) - temp1) < 1e-6) / level1 > 0.85
            index3 = np.sum(np.abs(np.round(temp1, 0) - temp2_interp) < 1e-6) / level1 > 0.85 or np.sum(np.abs(np.round(temp2_interp, 0) - temp1) < 1e-6) / level1 > 0.85

            index4 = np.sum(np.abs(truncate(temp1, 2) - temp2_interp) < 1e-6) / level1 > 0.85 or np.sum(np.abs(truncate(temp2_interp, 2) - temp1) < 1e-6) / level1 > 0.85
            index5 = np.sum(np.abs(truncate(temp1, 1) - temp2_interp) < 1e-6) / level1 > 0.85 or np.sum(np.abs(truncate(temp2_interp, 1) - temp1) < 1e-6) / level1 > 0.85
            index6 = np.sum(np.abs(truncate(temp1, 0) - temp2_interp) < 1e-6) / level1 > 0.85 or np.sum(np.abs(truncate(temp2_interp, 0) - temp1) < 1e-6) / level1 > 0.85
        else:
            # depth1 interpolates to the depth of depth2
            temp1_interp=np.interp(depth2, depth1, temp1, left=np.nan, right=np.nan, period=None)
            level2=content1['depth_number']
            index1 = np.sum(np.abs(np.round(temp1_interp, 2) - temp2) < 1e-6) / level2 > 0.85 or np.sum(np.abs(np.round(temp2, 2) - temp1_interp) < 1e-6) / level2 > 0.85
            index2 = np.sum(np.abs(np.round(temp1_interp, 1) - temp2) < 1e-6) / level2 > 0.85 or np.sum(np.abs(np.round(temp2, 1) - temp1_interp) < 1e-6) / level2 > 0.85
            index3 = np.sum(np.abs(np.round(temp1_interp, 0) - temp2) < 1e-6) / level2 > 0.85 or np.sum(np.abs(np.round(temp2, 0) - temp1_interp) < 1e-6) / level2 > 0.85

            index4 = np.sum(np.abs(truncate(temp1_interp, 2) - temp2) < 1e-6) / level2 > 0.85 or np.sum(np.abs(truncate(temp2, 2) - temp1_interp) < 1e-6) / level2 > 0.85
            index5 = np.sum(np.abs(truncate(temp1_interp, 1) - temp2) < 1e-6) / level2 > 0.85 or np.sum(np.abs(truncate(temp2, 1) - temp1_interp) < 1e-6) / level2 > 0.85
            index6 = np.sum(np.abs(truncate(temp1_interp, 0) - temp2) < 1e-6) / level2 > 0.85 or np.sum(np.abs(truncate(temp2, 0) - temp1_interp) < 1e-6) / level2 > 0.85
        if (index1 or index2 or index3 or index4 or index5 or index6):
            duplicate_check4 = 1
            ### distance less than 500m
            juli_distance = distance(content1['latitude'], content2['latitude'], content1['longitude'],content2['longitude'])
            if (juli_distance > 0.5):
                duplicate_check4 = False

    #### Fifth check: Layer by layer check
    #### Compare deph1,depth2,depth_diff,temp1,temp2,temp_diff,sal_diff
    duplicate_check_combined=False
    if(content1['depth_number']==content2['depth_number']):
        dz = np.sum(np.abs(depth1-depth2))
        dt = np.sum(np.abs(temp1-temp2))
        ds = np.sum(np.abs(salinity1-salinity2))
        if(np.all(np.isnan(ds))):
            duplicate1 = dz < 0.1 and dt < 0.1
        else:
            duplicate1 = dz < 0.1 and dt < 0.1 and ds <0.1

        ### 85% the same
        levels=content1['depth_number']
        nz = np.sum(np.abs(depth1-depth2)<0.01)
        nt = np.sum(np.abs(temp1-temp2)<0.01)
        ns = np.sum(np.abs(salinity1-salinity2)<0.01)
        ds = np.sum(np.abs(salinity1-salinity2))
        percent_z=np.float(nz)/levels
        percent_t=np.float(nt)/levels
        percent_s=np.float(ns)/levels
        if(np.all(np.isnan(ds))):
            duplicate2=percent_z>0.85 and percent_t>0.85
        else:
            duplicate2=percent_z>0.85 and percent_t>0.85 and percent_s>0.85

        duplicate_check_combined=duplicate1 or duplicate2

        if(duplicate_check_combined==True):
            ### Add meteorological wind speed and direction information for judgment--has wind speed and direction information but not same--nonduplicated
            if((content1['Wind_Direction']!='' and content2['Wind_Direction']!='') or (content1['Wind_Speed']!=999 and content2['Wind_Speed']!=999)):
                if(np.abs(content1['Wind_Speed']-content2['Wind_Speed'])>1e-3 or content1['Wind_Direction'] != content2['Wind_Direction']): 
                    duplicate_check_combined=False

        ### The same values can be observed at different times (when the water mass is relatively stable), especially in the Mediterranean Sea, Red Sea, Black Sea, Persian Gulf and other high temperature and high salt water
        if(duplicate_check_combined==True):
            # Use the datetime libraries to calculate the time difference
            date1 = datetime(content1['year'], content1['month'], content1['day'], content1['hour'], content1['minute'])
            date2 = datetime(content2['year'], content2['month'], content2['day'], content2['hour'], content2['minute'])
            duration = date1 - date2
            if duration.days < 0:
                duration = date2 - date1
            days = abs(duration.days)
            hours = abs(duration.seconds / 3600)
            if (days!=0 or hours!=0):
                ### different time-- threshold: less than 30 days
                if (days<30):
                    ### Located in the Mediterranean Sea: Latitude 30~40, longitude -5~35
                    if(content1['latitude']>30 and content1['latitude']<40 and content2['latitude']>30 and content2['latitude']<40 and content1['longitude']>-5 and content1['longitude']<35 and content2['longitude']>-5 and content2['longitude']<35 ):
                        duplicate_check_combined = False
                    ### Located on the Red Sea: latitude 12~31, longitude 32~44
                    if (content1['latitude']>12 and content1['latitude']<31 and content2['latitude']>12 and content2['latitude']<31 and content1['longitude']>32 and content1['longitude']<44 and content2['longitude']>32 and content2['longitude']<44):
                        duplicate_check_combined = False
                    ### Located on the Black Sea: latitude 41~58, longitude 28~42
                    if (content1['latitude']>41 and content1['latitude']<58 and content2['latitude']>41 and content2['latitude']<58 and content1['longitude']>28 and content1['longitude']<42 and content2['longitude']>28 and content2['longitude']<42):
                        duplicate_check_combined = False
                    ### Located in the Persian Gulf: latitude 23~30, longitude 50~60
                    if (content1['latitude']>23 and content1['latitude']<30 and content2['latitude']>23 and content2['latitude']<30 and content1['longitude']>50 and content1['longitude']<60 and content2['longitude']>50 and content2['longitude']<60):
                        duplicate_check_combined = False

        ### different ship and different location--nonduplicated
        if(content1['Platform']!=content2['Platform']):
            if(np.abs(content1['latitude'] - content2['latitude']) >= 0.01 or np.abs(content1['longitude'] - content2['longitude']) >= 0.01):
                duplicate_check_combined = False
        if(content1['Platform']!='' and content2['Platform']!='' and content1['Platform']==content2['Platform'] and '(' not in content1['Platform'] and content1['WOD_cruise_identifier']!=content2['WOD_cruise_identifier']):
            if (np.abs(content1['latitude'] - content2['latitude']) >= 0.01 or np.abs(content1['longitude'] - content2['longitude']) >= 0.01):
                duplicate_check_combined = False

        if(duplicate_check_combined==True):
            ### the type if exact duplicated
            if(np.abs(content1['latitude']-content2['latitude'])>=0.01 or np.abs(content1['longitude']-content2['longitude'])>=0.01):
                duplicate_check5=True    # wrong location
            if(content1['year']!=content2['year'] or content1['month']!=content2['month'] or content1['day']!=content2['day']):
                duplicate_check6=True    # wrong date
            if(content1['hour']!=content2['hour'] or content1['minute']!=content2['minute']):
                duplicate_check7=True   # wrong time
            if(content1['country_id']!=content2['country_id']):
                duplicate_check8=True   # wrong country code
            if(content1['dataset_id']!=content2['dataset_id']):
                duplicate_check9=True  # wrong instrument type

    #### Sixth check: Accurate repeat check
    #### All the records and meta-data are identiical(both empty strings are considered identical)
    if(duplicate_check_combined==True):
        index1=content1['access_no']==content2['access_no'] and content1['dataset_id']==content2['dataset_id']
        index2=np.abs(content1['latitude']-content2['latitude'])<1e-4 and np.abs(content1['longitude']-content2['longitude'])<1e-4
        index3=content1['year']==content2['year'] and content1['month']==content2['month'] and content1['day']==content2['day']
        index4=content1['hour']==content2['hour'] and content1['minute']==content2['minute']
        index5=np.abs(content1['sum_depth']-content2['sum_depth'])<1e-4 and np.abs(content1['sum_temp']-content2['sum_temp'])<1e-4 and np.abs(content1['sum_salinity']-content2['sum_salinity'])<1e-4
        index6=content1['probe_type']==content2['probe_type'] and content1['recorder']==content2['recorder'] and content1['maximum_depth']==content2['maximum_depth']
        index7=content1['country_id']==content2['country_id'] and content1['GMT_time']==content2['GMT_time']
        index8=content1['dbase_orig']==content2['dbase_orig'] and content1['project_name']==content2['project_name'] and content1['Platform']==content2['Platform']
        index9=content1['ocean_vehicle']==content2['ocean_vehicle'] and content1['Institute']==content2['Institute'] and content1['WOD_cruise_identifier']==content2['WOD_cruise_identifier']
        index10=np.abs(content1['std_depth']-content2['std_depth'])<1e-4 and np.abs(content1['std_temp']-content2['std_temp'])<1e-4 and np.abs(content1['std_salinity']-content2['std_salinity'])<1e-4
        ### Add meteorological wind speed and direction information for judgment
        index11=np.abs(content1['Wind_Speed']-content2['Wind_Speed'])<1e-3 and content1['Wind_Direction'] == content2['Wind_Direction']

        if(index1 and index2 and index3 and index4 and index5 and index6 and index7 and index8 and index9 and index10 and index11):
            duplicate_check10=True


    #### Seventh check: Interpolation check
    #### If the depth layers are not equal, determine whether the data has been interpolated(interpolation method: linear interpolation)
    if(content1['depth_number']!=content2['depth_number']):
        duplicate_check11_T=False
        duplicate_check11_S=False
        ###Temperature
        if(content1['depth_number'] < content2['depth_number']):
            # depth2 interpolates to the depth of depth1, temp2 compare with temp1 layer by layer
            temp2_interp=np.interp(depth1, depth2, temp2, left=np.nan, right=np.nan, period=None)
            nt = np.sum(np.abs(temp2_interp-temp1)<0.01)
            levels=content1['depth_number']
        else:
            # depth1 interpolates to the depth of depth2, temp1 compare with temp2 layer by layer
            temp1_interp=np.interp(depth2, depth1, temp1, left=np.nan, right=np.nan, period=None)
            nt = np.sum(np.abs(temp1_interp-temp2)<0.01)
            levels=content2['depth_number']
        percent_t=np.float(nt)/levels            
        duplicate_check11_T=percent_t>0.85

        ### Salinity
        if(content1['hasSalinity'] and content2['hasSalinity']):
            if(content1['depth_number'] < content2['depth_number']):
                # depth2 interpolates to the depth of depth1, salinity2 compare with salinity1 layer by layer
                salinity2_interp=np.interp(depth1, depth2, salinity2, left=np.nan, right=np.nan, period=None)
                nt = np.sum(np.abs(salinity2_interp-salinity1)<0.001) # threshold value: 0.001
                levels=content1['depth_number']
            else:
                # depth1 interpolates to the depth of depth2, salinity1 compare with salinity2 layer by layer
                salinity1_interp=np.interp(depth2, depth1, salinity1, left=np.nan, right=np.nan, period=None)
                nt = np.sum(np.abs(salinity1_interp-salinity2)<0.001) # threshold value: 0.001
                levels=content2['depth_number']
            percent_s=np.float(nt)/levels            
            duplicate_check11_S=percent_s>0.85

        duplicate_check11=duplicate_check11_T or duplicate_check11_S

        ### The GPS accuracy of CTD is very high, and if the longitude and latitude are different after three or four decimal places, it can not be determined as the same place
        if ('CTD' in content1['dataset_name'] and 'CTD' in content2['dataset_name']) or ('ctd' in content1['dataset_name'] and 'ctd' in content2['dataset_name']):
            lat_num1 = count_decimal_places(content1['latitude'])
            lon_num1 = count_decimal_places(content1['longitude'])
            lat_num2 = count_decimal_places(content2['latitude'])
            lon_num2 = count_decimal_places(content2['longitude'])
            if (lat_num1 > 2 and lon_num1 > 2 and lat_num2 > 2 and lon_num2 > 2):    # Latitude and longitude have the third decimal
                index1 = (np.abs(content1['latitude'] - content2['latitude']) > 1e-3) or (np.abs(content1['longitude'] - content2['longitude']) > 1e-3)  # latitude and longitude difference is in the third decimal
                if index1:
                    duplicate_check11 = False

    #### Eighth check: CTD double data check
    #### Two sensors of CTD work simultaneously two sets of data
    if ('CTD' in content1['dataset_name'] and 'CTD' in content2['dataset_name']) or ('ctd' in content1['dataset_name'] and 'ctd' in content2['dataset_name']):
        index1 = content1['year'] == content2['year'] and content1['month'] == content2['month'] and content1['day'] == content2['day'] and content1['hour'] == content2['hour'] and content1['minute'] == content2['minute']
        index2 = (content1['hour'] != 0 and content2['hour'] != 0 and content1['minute'] != 0 and content2['minute'] != 0) or (content1['hour'] == 0 and content2['hour'] == 0 and content1['minute'] != 0 and content2['minute'] != 0)
        index3 = (content1['ocean_vehicle'] == content2['ocean_vehicle'] and content1['ocean_vehicle'] != '') or (content1['WOD_cruise_identifier'] == content2['WOD_cruise_identifier'] and content1['WOD_cruise_identifier'] != '')
        juli = distance(content1['latitude'], content2['latitude'], content1['longitude'], content2['longitude'])
        index4 = juli <= 0.05  # distance less than 500m
        index6 = content1['access_no'] == content2['access_no']
        if (index1 and index2 and index3 and index4 and index6):
            duplicate_check12 = True
            ### The up and down data of instruments CTD are not duplicates
            if (np.abs(content1['latitude'] - content2['latitude']) < 0.1 and np.abs(content1['longitude'] - content2['longitude']) < 0.1):  # Approximate latitude and longitude
                if (content1['Cast_Direction'] != '' and content2['Cast_Direction'] != '' and content1['Cast_Direction'] != content2['Cast_Direction']):  # has Cast_Direction information; one is UP the other is DOWN
                    duplicate_check12 = False


    isDuplicate=duplicate_check1 or duplicate_check2 or duplicate_check3 or duplicate_check4 or duplicate_check5 or duplicate_check6 or duplicate_check7 or duplicate_check8 or duplicate_check9 or duplicate_check10 or duplicate_check11 or duplicate_check12

    if(isDuplicate == False):  # nonduplicated
        duplicate_multimodels=(duplicate_check1,duplicate_check2,duplicate_check3,duplicate_check4,duplicate_check5,duplicate_check6,duplicate_check7,duplicate_check8,duplicate_check9,duplicate_check10,duplicate_check11,duplicate_check12)
        return isDuplicate,duplicate_multimodels

    if(1):    
        print('\n')
        print("Spatial-temporal checks--Simultaneously but at different location: %d" %(duplicate_check1))
        print("Spatial-temporal checks--Simultaneously and co-located: %d" % (duplicate_check2))
        print("Correlation check: %d" % (duplicate_check3))
        print("Truncation check: %d" % (duplicate_check4))
        print("Layer by layer check-wrong location: %d" % (duplicate_check5))
        print("Layer by layer check-wrong date: %d" % (duplicate_check6))
        print("Layer by layer check-wrong time: %d" % (duplicate_check7))
        print("Layer by layer check-wrong country: %d" % (duplicate_check8))
        print("Layer by layer check-wrong instrument: %d" % (duplicate_check9))
        print("Accurate repeat check: %d" % (duplicate_check10))
        print("Interpolation(missing data) check: %d" % (duplicate_check11))
        print("CTD double data check: %d" % (duplicate_check12))

    if(isDuplicate==True):  #exat duplicated or near duplicated
        if(duplicate_check10==True):
            exact_duplicated=True
            isDuplicate=1  # exact duplicated
        else:
            near_duplicated=True
            isDuplicate=2  # near duplicated

    duplicate_multimodels=[duplicate_check1,duplicate_check2,duplicate_check3,duplicate_check4,duplicate_check5,duplicate_check6,duplicate_check7,duplicate_check8,duplicate_check9,duplicate_check10,duplicate_check11,duplicate_check12]

    return isDuplicate,duplicate_multimodels   ### only for N04_check_nc_duplicate_list4.py
    #return isDuplicate   ### for N04_check_nc_duplicate.py


def distance(lat1,lat2,lon1,lon2):
    r=6371 # earth radius
    N_S=r*math.pi/180   # The distance of the grid points in the north-south direction, in kilometers
    W_E=N_S*math.cos(lat1*math.pi/180)
    lon12=(lon2-lon1)*W_E
    lat12=(lat2-lat1)*N_S
    dis=math.sqrt(math.pow(lat12,2)+math.pow(lon12,2))   # In kilometers
    return dis

def truncate(value,decimal):
    beishu=math.pow(10,decimal)
    int_array=value * beishu
    int_array=int_array.astype(np.int)
    result=int_array / beishu
    return result

def check_Mode(array1,array2,thresholds=0.85):
    # division
    data_divide=array1/array2
    data_divide=data_divide[~np.isnan(data_divide)]  #Discard the case where the divisor is 0
    levels=len(data_divide)
    vals,counts=np.unique(data_divide, return_counts=True)
    counts=counts[vals!=1.]  # Not equal to 1 means the values are not equal
    isMode1=np.any((counts/levels)>thresholds)

    #subtraction
    levels=len(array1)
    vals,counts=np.unique(array1-array2, return_counts=True)
    counts=counts[vals!=0.]
    isMode2=np.any((counts/levels)>thresholds)

    isMode=isMode1 or isMode2
    return isMode

### Converts decimal to degree minutes and seconds
def dd2dms(dd):
    degree=(int)(float(dd))
    minute=(int)((float(dd)-degree)*60)
    second=(float(dd)-degree-float(minute)/60)*3600
    dms=str(degree)+' '+str(minute)+ ' '+str(second)
    return dms

def dms2dd(degree,minute,second):
    dd=degree + minute / 60 + second / 3600
    return dd

def getGMT(lat,lon):
    tf = TimezoneFinder(in_memory=True)
    local_time_zone = tf.timezone_at(lng=lon, lat=lat)
    return local_time_zone

def count_decimal_places(num):
    # Converts a floating point number to a string
    num_str = str(num)
    # If there is no decimal point, return 0
    if '.' not in num_str:
        return 0
    # Split the string into integer and decimal parts
    _, decimal_part = num_str.split('.')
    # Returns the length of the decimal part
    return len(decimal_part)




