import numpy as np
import math
from timezonefinder import TimezoneFinder
from datetime import datetime


def flag_duplicate(content1,content2):
    temp1=content1['temp']
    depth1=content1['depth']
    salinity1=content1['sal']

    temp2=content2['temp']
    depth2=content2['depth']
    salinity2=content2['sal']

    flag1=0
    flag2=0
    ###Step1: 判断是否完全一样（所有信息一样）：若一样则随机选择一个
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
    # print(index1,index2,index3,index4,index5,index6,index7,index8,index9,index10)
    if(index1 and index2 and index3 and index4 and index5 and index6 and index7 and index8 and index9 and index10):
        flag1=1
        flag2=0
        print('All identical records')  #完全一样（所有信息一样）
        return flag1,flag2
    
    #判断仪器是否相同类型
    if(content1['dataset_id']!=content2['dataset_id']): 
        #如果不同，保留仪器质量更好的数据 XBT vs Baty 保留XBT； CTD vs TESAC 保留CTD；  CTDvsBottle保留CTD; CTD和Argo 保留Argo XBT v.s. CTD/Argo 都保留
        name1=content1['dataset_name'].lower()
        name2=content2['dataset_name'].lower()
        # XBT vs MBT
        if('xbt' in name1 and 'mbt' in name2):
            flag1=0;flag2=1
        elif('mbt' in name1 and 'xbt' in name2):
            flag1=1;flag2=0
        # CTD vs Bottle
        if('ctd' in name1 and 'bottle' in name2):
            flag1=0;flag2=1
        elif('bottle' in name1 and 'ctd' in name2):
            flag1=1;flag2=0     
        # CTD和Argo 保留Argo
        if('ctd' in name1 and ('pfl' in name2 or 'argo' in name2)):
            flag1=1;flag2=0
        elif(('pfl' in name1 or 'argo' in name1) and 'ctd' in name2):
            flag1=0;flag2=1
        return flag1,flag2

    #如果仪器类型相同，判断温度、盐度、溶解氧、叶绿素等变量，保留变量多的那个
    number_variable1=sum(content1['hasTemp'],content1['hasSalinity'],content1['hasOxygen'],content1['hasChlonophyll'])
    number_variable2=sum(content2['hasTemp'],content2['hasSalinity'],content2['hasOxygen'],content2['hasChlonophyll'])
    if(number_variable1>number_variable2):
        flag1=0;flag2=1
    elif(number_variable1<number_variable2):
        flag1=1;flag2=0        

    #如果仪器类型相同,判断有效元数据信息个数（非空字符串），保留元数据信息个数最多的那个
    number_valid_meta1=cal_valid_meta(content1['accession_number'],content1['probe_type'],content1['recorder'],content1['country_id'],content1['GMT_time'],content1['dbase_orig'],content1['project_name'],content1['Platform'],content1['ocean_vehicle'],content1['Institute'],content1['WOD_cruise_identifier'])
    number_valid_meta2=cal_valid_meta(content2['accession_number'],content2['probe_type'],content2['recorder'],content2['country_id'],content2['GMT_time'],content2['dbase_orig'],content2['project_name'],content2['Platform'],content2['ocean_vehicle'],content2['Institute'],content2['WOD_cruise_identifier'])


# def cal_valid_meta(*variables):
#     for variable in variables:
        

def compair(content1,content2):
    ###############以后这边是可以改进的，继续做检查！！！！ 逐个逐个检查，排除不可能的情况，依据已知的重复情况，去做排除
    temp1=content1['temp']
    depth1=content1['depth']
    salinity1=content1['sal']

    temp2=content2['temp']
    depth2=content2['depth']
    salinity2=content2['sal']

    #不重复：0  准确重复1  疑似重复2
    isDuplicate=False

    ### 2023.1.2 added flag
    near_duplicated=False
    exact_duplicated=False

    duplicate_check1=False  #时空检查:同时不同地
    duplicate_check2=False  #时空检查：同时同地
    duplicate_check3=False  #相关性检查
    duplicate_check4=False   #截断/四舍五入检查
    duplicate_check5=False  #观测层检查：错误地点 （depth-by-depth check）
    duplicate_check6=False  #观测层检查：错误日期 （depth-by-depth check）
    duplicate_check7=False  #观测层检查：错误时间 （depth-by-depth check）
    duplicate_check8=False  #观测层检查：错误国家代码 （depth-by-depth check）
    duplicate_check9=False  #观测层检查：错误仪器代码 （depth-by-depth check）
    duplicate_check10=False #准确重复检查 （depth-by-depth check）
    duplicate_check11=False #插值检查
    duplicate_check12=False #CTD同时同地同一仪器多个观测
    
    #### 第一个检查：同时不同地--走航观测同一航次同一条船在相同时间却在不同位置出现两个观测（不合理） 唯一不同的是地理位置不同 （走航观测同一时间出现在不同的位置上，且hour 不为0）
    if(content1['dataset_id'] in [1,2,3,4] and content2['dataset_id'] in [1,2,3,4]):  #XBT MBT BOTTLE CTD
        index1=content1['year'] == content2['year'] and content1['month'] == content2['month'] and content1['day'] == content2['day'] and content1['hour'] == content2['hour'] and content1['minute'] == content2['minute']
        index2=content1['hour'] != 0 and content2['hour']!=0
        # index3=content1['ocean_vehicle']==content2['ocean_vehicle'] and content1['Platform']==content2['Platform'] and content1['Institute']==content2['Institute'] and content1['WOD_cruise_identifier']==content2['WOD_cruise_identifier']  #是否需要删除wod_cruise_identifier？
        index3=content1['ocean_vehicle']==content2['ocean_vehicle'] and content1['Platform']==content2['Platform'] and '(' in content1['Platform'] and content1['Institute']==content2['Institute']
        if(index1 and index2 and index3):  #地理位置不相等，设定大于0.01
            ########如果plarform code相等，但是ocean_vehicle为空或者WOD_cruise_identifier为空，就要排除，因为一个platform code可能代表不同的船
            # if(content1['ocean_vehicle']=='' or content1['WOD_cruise_identifier']==''):
            #     duplicate_check1=False
            ###2023.3.25 修改关于判断是否是同一条船的判断
            if(content1['Platform']!=''):
                if ('(' not in content1['Platform']!=''):
                    duplicate_check1 = False
            else:  ### 否则可以继续判断是否真的是重复
                if(content1['minute'] != 0 and content2['minute']!=0):   #如果分钟有记录
                      if(np.abs(content1['latitude']-content2['latitude'])>0.01 or np.abs(content1['longitude']-content2['longitude'])>0.01):  #阈值1km
                            duplicate_check1=True  #疑似重复
                else:  #如果小时相等，但是分钟没有记录的话
                    juli=distance(content1['latitude'],content2['latitude'],content1['longitude'],content2['longitude'])
                    if(juli>30):   #以一小时为30km的最大航速算 1个小时内两个位置最大距离不太可能超过30km，如果超过了，那就是不太可能发生，就有理由相信可能是重复的
                        duplicate_check1=True  #疑似重复

    #####第二个检查：同时同地--（不考虑浮标）相同时间同一时刻（具体到分钟）、相同位置 出现多个观测
    if (content1['dataset_id'] not in [5,7,9] and content2['dataset_id'] not in [5,7,9]):
        index1=content1['year'] == content2['year'] and content1['month'] == content2['month'] and content1['day'] == content2['day'] and content1['hour'] == content2['hour'] and content1['minute'] == content2['minute']
        index2=(content1['hour'] != 0 and content2['hour']!=0 and content1['minute'] != 0 and content2['minute'] != 0) or (content1['hour'] == 0 and content2['hour']==0 and content1['minute'] != 0 and content2['minute'] != 0)
        index3=content1['ocean_vehicle']==content2['ocean_vehicle'] and content1['Platform']==content2['Platform']
        juli=distance(content1['latitude'],content2['latitude'],content1['longitude'],content2['longitude'])  #距离小于1km，容错率提高到1000米，因为可能比较老的仪器的经纬度记录没有那么高精度，可能只保留到小数点后两位（分辨率为0.01度），这个大约是1km
        index4=juli<=1  #距离相差小于1000米  ###2023.1.3 更新提高到1000米（提高容错率），而不是100米
        index5=(content1['depth_number']>1) and (content2['depth_number']>1)
        # print(index1,index2,index3,index4)
        if(index1 and index2 and index3 and index4 and index5):
            duplicate_check2=True

        #####2023.02.20 对于dataset不是bot的廓线，同一航线（有船只识别码）、同地、同时（月份相差在一个月以内）的廓线对也判定为重复
        #####2023.03.12 更改对于同一航线的判断：plantform信息存在且包含"()"内容视为存在唯一识别符，此时若plantform信息相同则代表同一航线
        #if (content1['dataset_id'] != 1 and content2['dataset_id'] != 1) and (content1['WOD_cruise_identifier'] != '' and content2['WOD_cruise_identifier'] != '') and (content1['year'] == content2['year']):
        if (content1['dataset_id'] != 1 and content2['dataset_id'] != 1) and (content1['Platform'] != '' and content2['Platform'] != '') and (content1['year'] == content2['year']):
            index1 = juli <= 0.01  # 这样就保证同地了
        #### 2023.2.23（补充）需要补充：index2判断月份是否相差在一个月之内 用time, datetime库来计算，两个廓线的时间数组相减判断是否小于1个月？
            date1 = datetime(content1['year'], content1['month'], content1['day'], content1['hour'],content1['minute'])
            date2 = datetime(content2['year'], content2['month'], content2['day'], content2['hour'],content2['minute'])
            duration=date1-date2
            if duration.days < 0:
                duration = date2 - date1
            #print(duration)
            days=abs(duration.days)
            hours=abs(duration.seconds/3600)
            #print(days)
            #print(hours)
            index2 = days ==0
            index3 = hours <= 1
            index4 = (content1['hour'] != 0 and content2['hour'] != 0 and content1['minute'] != 0 and content2['minute'] != 0)
            #####2023.03.12 更改对于同一航线的判断：plantform信息存在且包含'('')'内容视为存在唯一识别符，此时若plantform信息相同则代表同一航线
            #index5 = content1['WOD_cruise_identifier'] == content2['WOD_cruise_identifier']
            index5 = ('(' in content1['Platform'] and '(' in content2['Platform'] and content1['Platform'] == content2['Platform'])
            index6 = content1['ocean_vehicle'] == content2['ocean_vehicle'] and content1['Platform'] == content2['Platform']
            index7 = content1['depth_number'] > 1 and content2['depth_number'] > 1
            if (index1 and index2 and index4 and index3 and index5 and index6 and index7):  # 浮点数不能直接用等号去判断相等，要绝对值做差！！！！！
                duplicate_check2 = True

        ##### 2023.1.3  加上气象风速风向信息作为判断 如果气象信息如果存在（不为空）的话，如果不相同，则duplicate_check2设置为False
        if((content1['Wind_Direction']!='' and content2['Wind_Direction']!='') or (content1['Wind_Speed']!=999 and content2['Wind_Speed']!=999)):
            if(np.abs(content1['Wind_Speed']-content2['Wind_Speed'])>1e-3 or content1['Wind_Direction'] != content2['Wind_Direction']):
                duplicate_check2=False

            ##### 2023.2.14 T5的XBT数据最大观测深度一般是超过1500米的
            # if (content1['probe_type']=='XBT: T5 (SIPPICAN)' and content1['maximum_depth']<1500) :
            #     duplicate_check2 = False
            # if (content2['probe_type']=='XBT: T5 (SIPPICAN)' and content2['maximum_depth']<1500) :
            #     duplicate_check2 = False

    ###### 第三个检查：相关性检查--深度或者温度或者盐度同时加上一个整数的，或者乘上缩小放大一个倍数的的重复
    if(np.abs(content1['cor_temp_depth']-content2['cor_temp_depth'])<1e-3 or (np.abs(content1['cor_sal_depth']-content2['cor_sal_depth'])<1e-4) and content1['cor_sal_depth']!=999 and content2['cor_sal_depth']!=999 and  content1['cor_temp_depth']!=999 and content2['cor_temp_depth']!=999):
        if(content1['dataset_id'] == content2['dataset_id'] and content1['depth_number']==content2['depth_number']):
            ##### 2023.2.13 增加对位置的限制
            juli_distance = distance(content1['latitude'], content2['latitude'], content1['longitude'], content2['longitude'])
            if(juli_distance<0.5):    #距离相差小于500米
                ##### 2023.2.13 温度数据去除掉 nan 数据后依然有观测点数依然大于 2
                if(content1['depth_number']>2 and content2['depth_number']>2 and len(temp1[~np.isnan(temp1)])>2 and len(temp2[~np.isnan(temp2)])>2):
                    isMode_depth=check_Mode(depth1,depth2,0.85)
                    isMode_temp=check_Mode(temp1,temp2,0.85)
                    # print(isMode_temp)
                    # print(isMode_depth)
                    if(isMode_temp or isMode_depth):
                        duplicate_check3=True

                    ##### 2023.2.13 盐度数据去除掉 nan 数据后依然有观测点数依然大于 2
                    if(content1['cor_sal_depth']!=999 and content2['cor_sal_depth']!=999 and len(salinity1[~np.isnan(salinity1)])>2 and len(salinity2[~np.isnan(salinity2)])>2):
                        isMode_salinity=check_Mode(salinity1,salinity2,0.85)
                        if(isMode_salinity):
                            duplicate_check3=True
            

    ##### 第4个检查： 截断检查---温度/盐度记录被四舍五入/截断造成的重复
    ##### 2023.3.10 加入一个限定条件：层数不能太少（可以先设定为大于3个，至少4个）
    if(content1['depth_number']==content2['depth_number'] and np.abs(content1['sum_depth']-content2['sum_depth'])<1 and content1['depth_number']>3 and content2['depth_number']>3):
        level=content1['depth_number']
        if(np.all(np.abs(temp1-temp2)<1e-5)):  #排除数值完全相等的exact duplicates
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
    ##2023.1.4针对观测点数不同，温度/盐度记录被四舍五入/截断造成的重复
    ##### 2023.3.10 加入一个限定条件：层数不能太少（可以先设定为大于3个，至少4个）
    elif (content1['depth_number']!=content2['depth_number'] and content1['depth_number']>3 and content2['depth_number']>3):
        if (content1['depth_number'] > content2['depth_number']):
            ##depth2插值到depth1所在深度
            temp2_interp = np.interp(depth1, depth2, temp2, left=np.nan, right=np.nan, period=None)
            # print(temp1)
            # print(temp2_interp)
            level1=content2['depth_number']
            index1 = np.sum(np.abs(np.round(temp1, 2) - temp2_interp) < 1e-6) / level1 > 0.85 or np.sum(np.abs(np.round(temp2_interp, 2) - temp1) < 1e-6) / level1 > 0.85
            index2 = np.sum(np.abs(np.round(temp1, 1) - temp2_interp) < 1e-6) / level1 > 0.85 or np.sum(np.abs(np.round(temp2_interp, 1) - temp1) < 1e-6) / level1 > 0.85
            index3 = np.sum(np.abs(np.round(temp1, 0) - temp2_interp) < 1e-6) / level1 > 0.85 or np.sum(np.abs(np.round(temp2_interp, 0) - temp1) < 1e-6) / level1 > 0.85

            index4 = np.sum(np.abs(truncate(temp1, 2) - temp2_interp) < 1e-6) / level1 > 0.85 or np.sum(np.abs(truncate(temp2_interp, 2) - temp1) < 1e-6) / level1 > 0.85
            index5 = np.sum(np.abs(truncate(temp1, 1) - temp2_interp) < 1e-6) / level1 > 0.85 or np.sum(np.abs(truncate(temp2_interp, 1) - temp1) < 1e-6) / level1 > 0.85
            index6 = np.sum(np.abs(truncate(temp1, 0) - temp2_interp) < 1e-6) / level1 > 0.85 or np.sum(np.abs(truncate(temp2_interp, 0) - temp1) < 1e-6) / level1 > 0.85
        else:
            ##depth1插值到depth2所在深度
            temp1_interp=np.interp(depth2, depth1, temp1, left=np.nan, right=np.nan, period=None)
            # print(temp1_interp)
            # print(temp2)
            level2=content1['depth_number']
            index1 = np.sum(np.abs(np.round(temp1_interp, 2) - temp2) < 1e-6) / level2 > 0.85 or np.sum(np.abs(np.round(temp2, 2) - temp1_interp) < 1e-6) / level2 > 0.85
            index2 = np.sum(np.abs(np.round(temp1_interp, 1) - temp2) < 1e-6) / level2 > 0.85 or np.sum(np.abs(np.round(temp2, 1) - temp1_interp) < 1e-6) / level2 > 0.85
            index3 = np.sum(np.abs(np.round(temp1_interp, 0) - temp2) < 1e-6) / level2 > 0.85 or np.sum(np.abs(np.round(temp2, 0) - temp1_interp) < 1e-6) / level2 > 0.85

            index4 = np.sum(np.abs(truncate(temp1_interp, 2) - temp2) < 1e-6) / level2 > 0.85 or np.sum(np.abs(truncate(temp2, 2) - temp1_interp) < 1e-6) / level2 > 0.85
            index5 = np.sum(np.abs(truncate(temp1_interp, 1) - temp2) < 1e-6) / level2 > 0.85 or np.sum(np.abs(truncate(temp2, 1) - temp1_interp) < 1e-6) / level2 > 0.85
            index6 = np.sum(np.abs(truncate(temp1_interp, 0) - temp2) < 1e-6) / level2 > 0.85 or np.sum(np.abs(truncate(temp2, 0) - temp1_interp) < 1e-6) / level2 > 0.85
        if (index1 or index2 or index3 or index4 or index5 or index6):
            duplicate_check4 = 1
            ##### 2023.2.14 增加对位置的限制
            juli_distance = distance(content1['latitude'], content2['latitude'], content1['longitude'],content2['longitude'])
            if (juli_distance > 0.5):  # 距离相差小于500米
                duplicate_check4 = False

    ##第5个检查：廓线逐层检查: 如果观测点数相同的话，可以比较，deph1,depth2,depth_diff,temp1,temp2,temp_diff,sal_diff
	####temp_diff 小于阈值0.01，（看compair.py）
    duplicate_check_combined=False
    if(content1['depth_number']==content2['depth_number']):
        dz = np.sum(np.abs(depth1-depth2))
        dt = np.sum(np.abs(temp1-temp2))
        ds = np.sum(np.abs(salinity1-salinity2))
        if(np.all(np.isnan(ds))):
            duplicate1 = dz < 0.1 and dt < 0.1
        else:
            duplicate1 = dz < 0.1 and dt < 0.1 and ds <0.1

        #算相同的百分比85%
        levels=content1['depth_number']
        nz = np.sum(np.abs(depth1-depth2)<0.01)  #这个包括四舍五入的精度
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
            ###### 2023.1.2 加上气象信息作为限制条件,如果气象信息存在（不为空）且不相同，则不用下面的判断，自动判定为不重复
            if((content1['Wind_Direction']!='' and content2['Wind_Direction']!='') or (content1['Wind_Speed']!=999 and content2['Wind_Speed']!=999)):
                if(np.abs(content1['Wind_Speed']-content2['Wind_Speed'])>1e-3 or content1['Wind_Direction'] != content2['Wind_Direction']): 
                    duplicate_check_combined=False
        ####2023.3.11 增加:不同时间相邻附近位置可以观测到相同的数值，在某个时间段内水团性质导致确实可能有相同观测值的现象，尤其在地中海、红海、黑海、波斯湾等高温高盐的海水
        if(duplicate_check_combined==True):
            ### 用时间戳判定时间是否相等
            date1 = datetime(content1['year'], content1['month'], content1['day'], content1['hour'], content1['minute'])
            date2 = datetime(content2['year'], content2['month'], content2['day'], content2['hour'], content2['minute'])
            duration = date1 - date2
            if duration.days < 0:
                duration = date2 - date1
            days = abs(duration.days)
            hours = abs(duration.seconds / 3600)
            if (days!=0 or hours!=0):
                #### 2023.03.14 不同时间限定在时间差小于30天
                if (days<30):
                    ###位于地中海：纬度30~40，经度-5~35
                    if(content1['latitude']>30 and content1['latitude']<40 and content2['latitude']>30 and content2['latitude']<40 and content1['longitude']>-5 and content1['longitude']<35 and content2['longitude']>-5 and content2['longitude']<35 ):
                        duplicate_check_combined = False
                    ###位于红海：纬度12~31，经度32~44
                    if (content1['latitude']>12 and content1['latitude']<31 and content2['latitude']>12 and content2['latitude']<31 and content1['longitude']>32 and content1['longitude']<44 and content2['longitude']>32 and content2['longitude']<44):
                        duplicate_check_combined = False
                    ###位于黑海：纬度41~58，经度28~42
                    if (content1['latitude']>41 and content1['latitude']<58 and content2['latitude']>41 and content2['latitude']<58 and content1['longitude']>28 and content1['longitude']<42 and content2['longitude']>28 and content2['longitude']<42):
                        duplicate_check_combined = False
                    ###位于波斯湾：纬度23~30，经度50~60
                    if (content1['latitude']>23 and content1['latitude']<30 and content2['latitude']>23 and content2['latitude']<30 and content1['longitude']>50 and content1['longitude']<60 and content2['longitude']>50 and content2['longitude']<60):
                        duplicate_check_combined = False
        ###### 2023.03.14 加上船只信息作为限制条件:如果船只信息不相等，同时测量位置不同，判定为不重复
        if(content1['Platform']!=content2['Platform']):
        #if (content1['WOD_cruise_identifier'] != '' and content2['WOD_cruise_identifier'] != '' and content1['WOD_cruise_identifier'] != content2['WOD_cruise_identifier']):
            if(np.abs(content1['latitude'] - content2['latitude']) >= 0.01 or np.abs(content1['longitude'] - content2['longitude']) >= 0.01):
                duplicate_check_combined = False
        if(content1['Platform']!='' and content2['Platform']!='' and content1['Platform']==content2['Platform'] and '(' not in content1['Platform'] and content1['WOD_cruise_identifier']!=content2['WOD_cruise_identifier']):
            if (np.abs(content1['latitude'] - content2['latitude']) >= 0.01 or np.abs(content1['longitude'] - content2['longitude']) >= 0.01):
                duplicate_check_combined = False
        ####2023.2.14 修改
        if(duplicate_check_combined==True):
            ########判断是什么类型的准确重复
            if(np.abs(content1['latitude']-content2['latitude'])>=0.01 or np.abs(content1['longitude']-content2['longitude'])>=0.01):
                duplicate_check5=True    #### wrong location
            if(content1['year']!=content2['year'] or content1['month']!=content2['month'] or content1['day']!=content2['day']):
                duplicate_check6=True    #### wrong date
            if(content1['hour']!=content2['hour'] or content1['minute']!=content2['minute']):
                duplicate_check7=True   #wrong time  (包括时区的问题？？？)
            if(content1['country_id']!=content2['country_id']):
                duplicate_check8=True   #wrong country code
            if(content1['dataset_id']!=content2['dataset_id']):
                duplicate_check9=True  #wrong instrument type

    #### 第6个检查：准确重复检查 （判断所有元数据信息和温度、深度观测值是否完全相同）判断是否完全是相同的，所有的元数据信息和观测值都是相同的，就算是空的字符串也算是相同
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
        ##### 2023.1.3  加上气象信息作为index11 限制条件,气象信息也需要相同
        index11=np.abs(content1['Wind_Speed']-content2['Wind_Speed'])<1e-3 and content1['Wind_Direction'] == content2['Wind_Direction']

        # print(index1,index2,index3,index4,index5,index6,index7,index8,index9,index10)
        if(index1 and index2 and index3 and index4 and index5 and index6 and index7 and index8 and index9 and index10 and index11):
            duplicate_check10=True


    #### 第7个检查：插值检查 -- 如果深度层数不相等，则判断数据是否被插值过（例如缺少表面值、缺少最后一个点的值、缺少中间几个数的值） 插值使用方法：线性插值
    if(content1['depth_number']!=content2['depth_number']):
        duplicate_check11_T=False
        duplicate_check11_S=False
        ##### 温度
        if(content1['depth_number'] < content2['depth_number']):
            #depth2 插值到depth1所在深度  temp2 <-- temp2_interp, 与temp1进行一一比较
            temp2_interp=np.interp(depth1, depth2, temp2, left=np.nan, right=np.nan, period=None)
            nt = np.sum(np.abs(temp2_interp-temp1)<0.01)
            levels=content1['depth_number']
        else:
            #depth1 插值到depth2所在深度  temp1 <-- temp1_interp, 与temp2进行一一比较
            temp1_interp=np.interp(depth2, depth1, temp1, left=np.nan, right=np.nan, period=None)
            nt = np.sum(np.abs(temp1_interp-temp2)<0.01)
            levels=content2['depth_number']
        percent_t=np.float(nt)/levels            
        duplicate_check11_T=percent_t>0.85

        ###2023.1.2 增加盐度的插值检查
        if(content1['hasSalinity'] and content2['hasSalinity']):
            if(content1['depth_number'] < content2['depth_number']):
                #depth2 插值到depth1所在深度  temp2 <-- temp2_interp, 与temp1进行一一比较
                salinity2_interp=np.interp(depth1, depth2, salinity2, left=np.nan, right=np.nan, period=None)
                ##### 2023.3.9改 盐度数据从0.01改成0.001
                nt = np.sum(np.abs(salinity2_interp-salinity1)<0.001)
                levels=content1['depth_number']
            else:
                #depth1 插值到depth2所在深度  temp1 <-- temp1_interp, 与temp2进行一一比较
                salinity1_interp=np.interp(depth2, depth1, salinity1, left=np.nan, right=np.nan, period=None)
                ##### 2023.3.9改 盐度数据从0.01改成0.001
                nt = np.sum(np.abs(salinity1_interp-salinity2)<0.001)
                levels=content2['depth_number']
            percent_s=np.float(nt)/levels            
            duplicate_check11_S=percent_s>0.85

        duplicate_check11=duplicate_check11_T or duplicate_check11_S

    ###### 2022.3.9 增加：第八个检查 CTD同一仪器两个传感器同时工作得到两组相同或不同的数据的重复情况
    if ('CTD' in content1['dataset_name'] and 'CTD' in content2['dataset_name']) or ('ctd' in content1['dataset_name'] and 'ctd' in content2['dataset_name']):
        index1 = content1['year'] == content2['year'] and content1['month'] == content2['month'] and content1['day'] == content2['day'] and content1['hour'] == content2['hour'] and content1['minute'] == content2['minute']
        index2 = (content1['hour'] != 0 and content2['hour'] != 0 and content1['minute'] != 0 and content2['minute'] != 0) or (content1['hour'] == 0 and content2['hour'] == 0 and content1['minute'] != 0 and content2['minute'] != 0)
        index3 = (content1['ocean_vehicle'] == content2['ocean_vehicle'] and content1['ocean_vehicle'] != '') or (content1['WOD_cruise_identifier'] == content2['WOD_cruise_identifier'] and content1['WOD_cruise_identifier'] != '')
        index4 = np.abs(content1['latitude'] - content2['latitude']) < 1e-4 and np.abs(content1['longitude'] - content2['longitude']) < 1e-4  # 经纬度相同
        index5 = np.abs(content1['latitude'] - math.trunc(content1['latitude'])) < 1e-4 or np.abs(content2['latitude'] - math.trunc(content2['latitude'])) < 1e-4 or np.abs(content1['longitude'] - math.trunc(content1['longitude'])) < 1e-4 or np.abs(content2['longitude'] - math.trunc(content2['longitude'])) < 1e-4  # 小数点后至少有一位有效数字
        index6 = content1['access_no'] == content2['access_no']
        if (index1 and index2 and index3 and index4 and index5 and index6):
            duplicate_check12 = True


    isDuplicate=duplicate_check1 or duplicate_check2 or duplicate_check3 or duplicate_check4 or duplicate_check5 or duplicate_check6 or duplicate_check7 or duplicate_check8 or duplicate_check9 or duplicate_check10 or duplicate_check11 or duplicate_check12
    # print(duplicate_check1, duplicate_check2,duplicate_check3,duplicate_check4,duplicate_check5)

    if(isDuplicate == False):  #不是重复数据
        duplicate_multimodels=(duplicate_check1,duplicate_check2,duplicate_check3,duplicate_check4,duplicate_check5,duplicate_check6,duplicate_check7,duplicate_check8,duplicate_check9,duplicate_check10,duplicate_check11,duplicate_check12)
        return isDuplicate,duplicate_multimodels

    if(1):    
        print('\n')
        print("时空检查--同时不同地: %d" %(duplicate_check1))
        print("时空检查--同时同地: %d" % (duplicate_check2))
        print("相关性检查: %d" % (duplicate_check3))
        print("截断/四舍五入检查: %d" % (duplicate_check4))
        print("观测层检查:错误地点 %d" % (duplicate_check5))
        print("观测层检查:错误日期：%d" % (duplicate_check6)) 
        print("观测层检查:错误时间：%d" % (duplicate_check7))
        print("观测层检查:错误国家代码：%d" % (duplicate_check8))
        print("观测层检查:错误仪器类型：%d" % (duplicate_check9))
        print("准确重复（完全一致）：%d" % (duplicate_check10))
        print("插值检查（部分数据缺失）：%d" % (duplicate_check11))  #其实这个中文名应该改成  “missing data check”  数据缺失检查会好一点，因为数据缺失不一定是由插值造成的
        ###### 2023.3.9 增加
        print("CTD同时同地同一仪器多个观测: %d" % (duplicate_check12))

    if(isDuplicate==True):  #判断是准确重复还是近似重复
        if(duplicate_check10==True):
            exact_duplicated=True
            isDuplicate=1  #exact_duplicated
        else:
            near_duplicated=True
            isDuplicate=2  #near_duplicated

    duplicate_multimodels=[duplicate_check1,duplicate_check2,duplicate_check3,duplicate_check4,duplicate_check5,duplicate_check6,duplicate_check7,duplicate_check8,duplicate_check9,duplicate_check10,duplicate_check11,duplicate_check12]
    # print(duplicate_multimodels)

    #return isDuplicate,duplicate_multimodels   ### only for N04_check_nc_duplicate_list4.py
    return isDuplicate   ### for N04_check_nc_duplicate.py



def distance(lat1,lat2,lon1,lon2):
    #####开平方根
    r=6371 #地球半径
    N_S=r*math.pi/180   #这里是一个南北方向格点的距离，以km为单位
    W_E=N_S*math.cos(lat1*math.pi/180)
    lon12=(lon2-lon1)*W_E
    lat12=(lat2-lat1)*N_S
    dis=math.sqrt(math.pow(lat12,2)+math.pow(lon12,2))   #以km为单位
    return dis

def truncate(value,decimal):
    beishu=math.pow(10,decimal)
    int_array=value * beishu
    int_array=int_array.astype(np.int)
    result=int_array / beishu
    # print(result)
    return result

def check_Mode(array1,array2,thresholds=0.85):
    #除法
    data_divide=array1/array2
    data_divide=data_divide[~np.isnan(data_divide)]  #去掉除数为0的情况
    levels=len(data_divide)
    vals,counts=np.unique(data_divide, return_counts=True)
    counts=counts[vals!=1.]  #不等于1表示不是相等数值
    isMode1=np.any((counts/levels)>thresholds)

    #减法
    levels=len(array1)
    vals,counts=np.unique(array1-array2, return_counts=True)
    counts=counts[vals!=0.]
    isMode2=np.any((counts/levels)>thresholds)

    isMode=isMode1 or isMode2
    return isMode

def dd2dms(dd):   ######将十进制转换成度分秒
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




