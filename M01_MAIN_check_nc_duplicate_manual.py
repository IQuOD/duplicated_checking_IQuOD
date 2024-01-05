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

class Duplicate_check(object):

    def __init__(self):
        pass
        #self.path='/Users/zqtzt/Downloads/WODdata' ######This netCDF filepath needs to be modified to suit your case
        #self.path=filepath
        #self.path=file_path ######This netCDF filepath needs to be modified to suit your case

    def run(self,netCDF_filepath):
        while True:
            print('---------Please input two netCDF files which are potential duplicates--------')
            file1=input('The first netCDF file name is: ').rstrip().lstrip()
            file2=input('The second netCDF file name is: ').rstrip().lstrip()
            isOutput_detail = input("Output profile information or not(1: Yes; 0: No)")
			
            index_str=file1.rfind('_')
            date1=file1[index_str-14:index_str-6]
            year1=date1[0:4]
            month1=date1[4:6]
            path1=os.path.join(netCDF_filepath,year1,month1)

            index_str=file2.rfind('_')
            date2=file2[index_str-14:index_str-6]
            year2=date2[0:4]
            month2=date2[4:6]
            path2=os.path.join(netCDF_filepath,year2,month2)
            filepath1=os.path.join(path1,file1)
            filepath2=os.path.join(path2,file2)

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
                print('Duplicate result is: Near Duplicate')
            else:
                print('Duplicate result is: Not Duplicate')

    def output_detail(self,content1,content2):
        temp1=content1['temp']
        depth1=content1['depth']
        salinity1=content1['sal']

        temp2=content2['temp']
        depth2=content2['depth']
        salinity2=content2['sal']

        print('\n')
        ### Different depth number, direct output depth1, temp1, depth2, temp2
        if(content1['depth_number']!=content2['depth_number']):
            levels=np.nanmax([content1['depth_number'],content2['depth_number']])
            print('%10s %10s %10s  %10s %10s %10s' %('depth1','temp1','salinity1','depth2','temp2','salinity2'))
            for i in range(levels):
                if(i<content1['depth_number'] and i<content2['depth_number']):
                    print('%10.3f %10.4f %10.4f  %10.3f %10.4f %10.4f' %(depth1[i],temp1[i],salinity1[i],depth2[i],temp2[i],salinity2[i]))
                elif (i<content1['depth_number'] and i>=content2['depth_number']):
                    print('%10.3f %10.4f %10.4f  %10s %10s %10s' %(depth1[i],temp1[i],salinity1[i],' ',' ',' '))
                elif (i>=content1['depth_number'] and i<content2['depth_number']):
                    print('%10s %10s %10s  %10.3f %10.4f %10.4f' %(' ',' ',' ',depth2[i],temp2[i],salinity2[i]))
        elif(content1['depth_number']==content2['depth_number']):
        ### Same depth number, output deph1, depth2, depth_diff, temp1, temp2, temp_diff, sal_diff
            depth_diff=depth1-depth2
            temp_diff=temp1-temp2
            sal_diff=salinity1-salinity2
            print('%10s %10s %10s %10s %10s %10s %10s' %('depth1','depth2','depth_diff','temp1','temp2','temp_diff','sal_diff'))
            for i in range(content1['depth_number']):
                print('%10.3f %10.3f %10.3f %10.4f %10.4f %10.4f %10.4f' %(depth1[i],depth2[i],depth_diff[i],temp1[i],temp2[i],temp_diff[i],sal_diff[i]))


    def output_info_pairs(self,content1,content2):
        # output info
        print('%20s:%20s , %20s' % ('WOD_id',content1['wod_unique_id'], content2['wod_unique_id']))
        print('%20s:%20d , %20d' % ('Acess_no',content1['access_no'], content2['access_no']))
        print('%20s:%20s , %20s' % ('Dataset',content1['dataset_name'], content2['dataset_name']))
        if (content1['dataset_id']  in [2,11,8] or content2['dataset_id']  in [2,11,8]):
            print('%20s:%20s , %20s' % ('Cast_Direction', content1['Cast_Direction'], content2['Cast_Direction']))
        print('%20s:%20.4f , %20.4f' %('Lat',content1['latitude'], content2['latitude']))
        print('%20s:%20.4f , %20.4f' %('Long',content1['longitude'], content2['longitude']))
        print('%20s:%20d , %20d' %('Year',content1['year'], content2['year']))
        print('%20s:%20d , %20d' %('Month',content1['month'], content2['month']))
        print('%20s:%20s , %20s' %('Day',content1['day'], content2['day']))
        print('%20s:%20s , %20s' %('Hour',content1['hour'], content2['hour']))
        print('%20s:%20s , %20s' %('Minute',content1['minute'], content2['minute']))
        print('%20s:%20.4f , %20.4f' %('sum_depth',content1['sum_depth'], content2['sum_depth']))
        print('%20s:%20.4f , %20.4f' %('sum_temp',content1['sum_temp'], content2['sum_temp']))
        print('%20s:%20.4f , %20.4f' %('sum_salinity',content1['sum_salinity'], content2['sum_salinity']))
        print('%20s:%20s , %20s' %('Probe type',content1['probe_type'], content2['probe_type']))
        print('%20s:%20s , %20s' %('Recorder',content1['recorder'], content2['recorder']))
        print('%20s:%20d , %20d' %('Depth_number',content1['depth_number'], content2['depth_number']))
        print('%20s:%20.3f , %20.3f' %('Maximum Depth',content1['maximum_depth'], content2['maximum_depth']))
        print('%20s:%20d , %20d' %('hasTemp',content1['hasTemp'], content2['hasTemp']))
        print('%20s:%20d , %20d' %('hasSalinity',content1['hasSalinity'], content2['hasSalinity']))
        print('%20s:%20d , %20d' % ('hasOxygen', content1['hasOxygen'], content2['hasOxygen']))
        print('%20s:%20d , %20d' % ('hasChlonophyII', content1['hasChlonophyII'], content2['hasChlonophyII']))
        print('%20s:%20d , %20d' %('Country',content1['country_id'], content2['country_id']))
        print('%20s:%20.3f , %20.3f' %('GMT_time',content1['GMT_time'], content2['GMT_time']))
        if('XBT' in content1['dataset_name']):
            print('%20s:%20s , %20s' %('XBT depth fix',content1['need_z_fix'], content2['need_z_fix']))
        print('%20s:%20s , %20s' %('Database_origin',content1['dbase_orig'], content2['dbase_orig']))
        print('%20s:%20s , %20s' %('Project_name',content1['project_name'], content2['project_name']))
        print('%20s:%20s , %20s' %('Platform',content1['Platform'], content2['Platform']))
        print('%20s:%20s , %20s' %('Vehicle',content1['ocean_vehicle'], content2['ocean_vehicle']))
        print('%20s:%20s , %20s' %('Institute',content1['Institute'], content2['Institute']))
        print('%20s:%20s , %20s' %('WOD_cruise_identifier',content1['WOD_cruise_identifier'], content2['WOD_cruise_identifier']))
        print('%20s:%20s , %20s' %('Wind_Direction',content1['Wind_Direction'], content2['Wind_Direction']))
        print('%20s:%20.4f , %20.4f' %('Wind_Speed',content1['Wind_Speed'], content2['Wind_Speed']))
        print('%20s:%20.4f , %20.4f' %('Std_depth',content1['std_depth'], content2['std_depth']))
        print('%20s:%20.4f , %20.4f' %('Std_temp',content1['std_temp'], content2['std_temp']))
        print('%20s:%20.4f , %20.4f' %('Std_salinity',content1['std_salinity'], content2['std_salinity']))
        print('%20s:%20.4f , %20.4f' %('Corr_temp&depth',content1['cor_temp_depth'], content2['cor_temp_depth']))
        print('%20s:%20.4f , %20.4f' %('Corr_sal&depth',content1['cor_sal_depth'], content2['cor_sal_depth']))


    def read_nc_data(self,file):
        with nc.Dataset(file,'r') as f:
            try:
                probe_type=str(nc.chartostring(f.variables['Temperature_Instrument'][:]))
            except:
                probe_type=''
            try:
                recorder=str(nc.chartostring(f.variables['Recorder'][:]))
            except:
                recorder=''

            wod_unique_id=f.variables['wod_unique_cast'][:]
            try:
                access_no=f.variables['Access_no'][:]
            except:
                access_no=999

            time=f.variables['time']
            dtime = nc.num2date(time[-1], time.units)
            year=dtime.year
            month=dtime.month
            day=dtime.day
            hour = dtime.hour
            minute = dtime.minute

            dataset_name=str(nc.chartostring(f.variables['dataset'][:]))
            dataset_id=self.find_order_dataset(dataset_name)

            latitude=round(float(f.variables['lat'][:]),4)
            longitude=round(float(f.variables['lon'][:]),4)

            depth=f.variables['z'][:]
            depth[np.logical_or(depth>12000, depth<-10)]=np.nan

            depth_number=len(depth)
            maximum_depth=depth[-1]
            sum_depth=round(np.nansum(depth),4)
            std_depth=round(np.nanstd(depth),4)
            if(np.isnan(sum_depth)):
                sum_depth=999
            if(np.isnan(std_depth)):
                std_depth=999

            try:
                temp=f.variables['Temperature'][:]
                temp[np.logical_or(temp>40, temp<-2.5)]=np.nan
                temp2=temp[~np.isnan(temp)]
                depth2=depth[~np.isnan(temp)]
                hasTemp = 1
                sum_temp=round(np.nansum(temp),4)   
                std_temp=round(np.nanstd(temp),4)
                cor_temp_depth=np.corrcoef(temp2,depth2)[1,0]
                cor_temp_depth=round(cor_temp_depth,5)
                if(np.isnan(cor_temp_depth)):
                    cor_temp_depth=999
            except:
                temp=np.full((depth_number,1),np.nan)
                hasTemp=0
                sum_temp=0
                std_temp=999
                cor_temp_depth=999

            try:
                sal=f.variables['Salinity'][:]
                sal[np.logical_or(sal>43, sal<0)]=np.nan
                sal2=sal[~np.isnan(sal)]
                depth2=depth[~np.isnan(sal)]
                hasSalinity = 1
                sum_salinity = round(np.nansum(sal),4)
                std_salinity=round(np.nanstd(sal),4)
                cor_sal_depth=round(np.corrcoef(sal2,depth2)[1,0],5)
                if(np.isnan(cor_sal_depth)):
                    cor_sal_depth=999
            except:
                sal=np.full((depth_number,1),np.nan)
                sum_salinity=0
                hasSalinity=0
                std_salinity=999
                cor_sal_depth=999

            try:
                oxy=f.variables['Oxygen']
                hasOxygen = 1
            except:
                hasOxygen=0

            try:
                Chlonophyll=f.variables['Chlorophyll']
                hasChlonophyII = 1
            except:
                hasChlonophyII=0


            country_name=str(nc.chartostring(f.variables['country'][:]))
            country_id=self.find_id_country(country_name)

            try:
                GMT_time=round(float(f.variables['GMT_time'][:]),2)
            except:
                GMT_time=0
            try:
                WMO_id=int(f.variables['WMO_ID'][:])
            except:
                WMO_id=999

            try:
                dbase_orig = str(nc.chartostring(f.variables['dbase_orig'][:]))
            except:
                dbase_orig=''
            try:
                project_name=str(nc.chartostring(f.variables['Project'][:]))
            except:
                project_name=''
            try:
                Platform=str(nc.chartostring(f.variables['Platform'][:]))
            except:
                Platform=''

            try:
                ocean_vehicle=str(nc.chartostring(f.variables['Ocean_Vehicle'][:]))
            except:
                ocean_vehicle=''

            try:
                Institute=str(nc.chartostring(f.variables['Institute'][:]))
            except:
                Institute=''

            try:
                WOD_cruise_identifier=str(nc.chartostring(f.variables['WOD_cruise_identifier'][:]))
            except:
                WOD_cruise_identifier=''

            try:
                need_z_fix=str(nc.chartostring(f.variables['need_z_fix'][:]))
            except:
                need_z_fix=''

            try:
                Wind_Direction=str(nc.chartostring(f.variables['Wind_Direction'][:]))
            except:
                Wind_Direction=''

            try:
                Wind_Speed=f.variables['Wind_Speed'][:]
            except:
                Wind_Speed=999

            try:
                Cast_Direction=str(nc.chartostring(f.variables['Cast_Direction'][:]))
                ### Converts the Cast Direction information to uppercase
                Cast_Direction.upper()
            except:
                Cast_Direction=''

        ### Save with dictionary
        t_parameters={}
        self.add_parameters(t_parameters,Institute=Institute,need_z_fix=need_z_fix,WOD_cruise_identifier=WOD_cruise_identifier,wod_unique_id=wod_unique_id,access_no=access_no,depth=depth,temp=temp,sal=sal,dataset_id=dataset_id,dataset_name=dataset_name,latitude=latitude,longitude=longitude,probe_type = probe_type, recorder = recorder, year=year,month=month,day=day,hour=hour,minute=minute,depth_number=depth_number,maximum_depth=maximum_depth,hasTemp=hasTemp, hasSalinity=hasSalinity,hasOxygen=hasOxygen,hasChlonophyII=hasChlonophyII,country_id=country_id,GMT_time=GMT_time,WMO_id=WMO_id,dbase_orig=dbase_orig,project_name=project_name,Platform=Platform,ocean_vehicle=ocean_vehicle,sum_depth=sum_depth,sum_temp=sum_temp,sum_salinity=sum_salinity,std_depth=std_depth,std_temp=std_temp,std_salinity=std_salinity,cor_temp_depth=cor_temp_depth,cor_sal_depth=cor_sal_depth,Wind_Direction=Wind_Direction,Wind_Speed=Wind_Speed,Cast_Direction=Cast_Direction)
        return t_parameters


    def find_id_country(self,country_name):
        country_name=country_name.upper()
        try:
            country_id=t_country.c_dict[country_name]
        except:
            country_id=999

        return country_id

    def add_parameters(self,params, **kwargs):
        return params.update(kwargs)
        #https://www.zhihu.com/question/42768955/answer/94798842

    def find_order_dataset(self,dataset_name):
        dataset_name=dataset_name.lower()
        if('bod' in dataset_name or 'bottle' in dataset_name or 'ocean station' in dataset_name):
            dataset_id=1
        elif('towed' in dataset_name or 'uor' in dataset_name or 'undulating' in dataset_name):
            dataset_id=10
        elif('ctd' in dataset_name or 'xctd' in dataset_name):
            dataset_id=2
        elif('mbt' in dataset_name or 'mechanica' in dataset_name or 'mb' in dataset_name):
            dataset_id=3
        elif('xbt' in dataset_name or 'xb' in dataset_name or 'expendable' in dataset_name):
            dataset_id=4
        elif('sur' in dataset_name or 'surface' in dataset_name):
            dataset_id=5
        elif('apb' in dataset_name or 'autonomous' in dataset_name or 'animal' in dataset_name):
            dataset_id=6
        elif('mrb' in dataset_name or 'moored' in dataset_name or 'tao' in dataset_name):
            dataset_id=7
        elif('pfl' in dataset_name or 'argo' in dataset_name or 'profiling' in dataset_name):
            dataset_id=8
        elif('drb' in dataset_name or 'drifting' in dataset_name):
            dataset_id=9
        elif('gld' in dataset_name or 'glider' in dataset_name):
            dataset_id=11
        elif('dbt' in dataset_name):
            dataset_id=12
        elif('std' in dataset_name):
            dataset_id=13
        elif('microbt' in dataset_name):
            dataset_id=14
        else:
            dataset_id=999
        return dataset_id


def main(netCDF_filepath):
    dc=Duplicate_check()
    dc.run(netCDF_filepath)

if __name__ == '__main__':
    netCDF_filepath='/Users/zqtzt/Downloads/WODdata'  ### You should revise this path to fit your own purposes
    main(netCDF_filepath)