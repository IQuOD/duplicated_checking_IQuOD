import numpy as np
import math
from timezonefinder import TimezoneFinder


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
    duplicate_check1=False
    duplicate_check2=False
    duplicate_check3=False
    duplicate_check4=False
    duplicate_check5=False
    duplicate_check12=False
    duplicate_check6,duplicate_check7,duplicate_check8,duplicate_check9,duplicate_check10,duplicate_check11=False,False,False,False,False,False

    #### 第一个检查：走航观测同一航次同一条船在相同时间却在不同位置出现两个观测（不合理） 唯一不同的是地理位置不同 （走航观测同一时间出现在不同的位置上，且hour 不为0）
    if(content1['dataset_id'] in [1,2,3,4] and content2['dataset_id'] in [1,2,3,4]):  #XBT MBT BOTTLE CTD
        index1=content1['year'] == content2['year'] and content1['month'] == content2['month'] and content1['day'] == content2['day'] and content1['hour'] == content2['hour'] and content1['minute'] == content2['minute']
        index2=content1['hour'] != 0 and content2['hour']!=0
        # index3=content1['ocean_vehicle']==content2['ocean_vehicle'] and content1['Platform']==content2['Platform'] and content1['Institute']==content2['Institute'] and content1['WOD_cruise_identifier']==content2['WOD_cruise_identifier']  #是否需要删除wod_cruise_identifier？
        index3=content1['ocean_vehicle']==content2['ocean_vehicle'] and content1['Platform']==content2['Platform'] and content1['Institute']==content2['Institute']
        ### 测试  index4=content1['access_no']!=306 and content2['access_no']!=306  ####306的这个需要谨慎排除
        if(index1 and index2 and index3):  #地理位置不相等，设定大于0.01
            ########如果plarform code相等，但是ocean_vehicle为空或者WOD_cruise_identifier为空，就要排除，因为一个platform code可能代表不同的船 
            if(content1['ocean_vehicle']=='' or content1['WOD_cruise_identifier']==''):
                duplicate_check1=False
            else:  ### 否则可以继续判断是否真的是重复
                if(content1['minute'] != 0 and content2['minute']!=0):   #如果分钟有记录
                      if(np.abs(content1['latitude']-content2['latitude'])>0.01 or np.abs(content1['longitude']-content2['longitude'])>0.01):  #阈值1km
                            duplicate_check1=2  #这里应该是疑似重复
                else:  #如果小时相等，但是分钟没有记录的话
                    juli=distance(content1['latitude'],content2['latitude'],content1['longitude'],content2['longitude'])
                    # print(juli)
                    if(juli>30):   #以一小时为30km的最大航速算 1个小时内两个位置最大距离不太可能超过30km，如果超过了，那就是不太可能发生，就有理由相信可能是重复的
                        duplicate_check1=2  #这里应该是疑似重复



    ###### 第二个检查：深度或者温度或者盐度同时加上一个整数的，或者乘上缩小放大一个倍数的的重复
    if(np.abs(content1['cor_temp_depth']-content2['cor_temp_depth'])<1e-3 or (np.abs(content1['cor_sal_depth']-content2['cor_sal_depth'])<1e-4) and content1['cor_sal_depth']!=999 and content2['cor_sal_depth']!=999 and  content1['cor_temp_depth']!=999 and content2['cor_temp_depth']!=999):
        if(content1['dataset_id'] == content2['dataset_id'] and content1['depth_number']==content2['depth_number']):
            if(content1['depth_number']>2 and content2['depth_number']>2):
                isMode_depth=check_Mode(depth1,depth2,0.9)
                isMode_temp=check_Mode(temp1,temp2,0.9)
                # print(isMode_temp)
                # print(isMode_depth)
                if(isMode_temp or isMode_depth):
                    duplicate_check3=2
                
                if(content1['cor_sal_depth']!=999 and content2['cor_sal_depth']!=999):
                    isMode_salinity=check_Mode(salinity1,salinity2,0.9)
                    if(isMode_salinity):
                        duplicate_check3=2
            
    #####第三个检查：（不考虑浮标）相同时间同一时刻（具体到分钟）、相同位置 出现多个观测
    if (content1['dataset_id'] not in [5,7,9] and content2['dataset_id'] not in [5,7,9]):
        index1=content1['year'] == content2['year'] and content1['month'] == content2['month'] and content1['day'] == content2['day'] and content1['hour'] == content2['hour'] and content1['minute'] == content2['minute']
        index2=(content1['hour'] != 0 and content2['hour']!=0 and content1['minute'] != 0 and content2['minute'] != 0) or (content1['hour'] == 0 and content2['hour']==0 and content1['minute'] != 0 and content2['minute'] != 0)
        index3=content1['ocean_vehicle']==content2['ocean_vehicle'] and content1['Platform']==content2['Platform']
        index4=np.abs(content1['latitude']-content2['latitude'])<=0.001 and np.abs(content1['longitude']-content2['longitude'])<=0.001  #100米
        index5=content1['depth_number']>1 and content2['depth_number']>1
        # print(index1,index2,index3,index4)
        if(index1 and index2 and index3 and index4 and index5):
            duplicate_check4=2

    ##### 第5个检查： 温度/盐度记录被四舍五入/截断造成的重复
    if(content1['depth_number']==content2['depth_number'] and np.abs(content1['sum_depth']-content2['sum_depth'])<1):
        level=content1['depth_number']
        if(np.all(np.abs(temp1-temp2)<1e-5)):  #排除数值完全相等的exact duplicates
            duplicate_check5=False
        else:
            index1=np.sum(np.abs(np.round(temp1,2)-temp2)<1e-6)/level>0.85 or np.sum(np.abs(np.round(temp2,2)-temp1)<1e-6)/level>0.85
            index2=np.sum(np.abs(np.round(temp1,1)-temp2)<1e-6)/level>0.85 or np.sum(np.abs(np.round(temp2,1)-temp1)<1e-6)/level>0.85
            index3=np.sum(np.abs(np.round(temp1,0)-temp2)<1e-6)/level>0.85 or np.sum(np.abs(np.round(temp2,0)-temp1)<1e-6)/level>0.85

            index4= np.sum(np.abs(truncate(temp1,2)-temp2)<1e-6)/level>0.85 or np.sum(np.abs(truncate(temp2,2)-temp1)<1e-6)/level>0.85
            index5= np.sum(np.abs(truncate(temp1,1)-temp2)<1e-6)/level>0.85 or np.sum(np.abs(truncate(temp2,1)-temp1)<1e-6)/level>0.85
            index6= np.sum(np.abs(truncate(temp1,0)-temp2)<1e-6)/level>0.85 or np.sum(np.abs(truncate(temp2,0)-temp1)<1e-6)/level>0.85
            if(index1 or index2 or index3 or index4 or index5 or index6):
                duplicate_check5=1


    # ###如果观测点数不相同的话，无法比较，直接屏幕输出 depth1,temp1,  depth2, temp2
    # if(content1['depth_number']!=content2['depth_number']):
    #     # isDuplicate=False
    #     return isDuplicate


    ##第6个检查：如果观测点数相同的话，可以比较，deph1,depth2,depth_diff,temp1,temp2,temp_diff,sal_diff
	####temp_diff 小于阈值0.01，（看compair.py）
    if(content1['depth_number']==content2['depth_number']):
        dz = np.sum(np.abs(depth1-depth2))
        dt = np.sum(np.abs(temp1-temp2))
        ds = np.sum(np.abs(salinity1-salinity2))
        if(np.all(np.isnan(ds))):
            duplicate1 = dz < 0.1 and dt < 0.1
        else:
            duplicate1 = dz < 0.1 and dt < 0.1 and ds <0.1

        #算相同的百分比75%
        levels=content1['depth_number']
        nz = np.sum(np.abs(depth1-depth2)<0.01)  #这个包括四舍五入的精度
        nt = np.sum(np.abs(temp1-temp2)<0.01)
        ns = np.sum(np.abs(salinity1-salinity2)<0.01)
        ds = np.sum(np.abs(salinity1-salinity2))
        percent_z=np.float(nz)/levels
        percent_t=np.float(nt)/levels
        percent_s=np.float(ns)/levels
        if(np.all(np.isnan(ds))):
            duplicate2=percent_z>0.75 and percent_t>0.75
        else:
            duplicate2=percent_z>0.75 and percent_t>0.75 and percent_s>0.75

        duplicate_check2=duplicate1 or duplicate2
        if(duplicate_check2==True):
            if(duplicate1==1):
                duplicate_check2=1
            else:
                duplicate_check2=2

            ########判断是什么类型的准确重复
            if(np.abs(content1['latitude']-content2['latitude'])>=0.01 or np.abs(content1['longitude']-content2['longitude'])>=0.01):
                duplicate_check6=1    #### wrong location
            elif(content1['year']!=content2['year'] or content1['month']!=content2['month'] or content1['day']!=content2['day']):
                duplicate_check7=1    #### wrong date
            elif(content1['hour']!=content2['hour'] or content1['minute']!=content2['minute']):
                duplicate_check8=1   #wrong time  (包括时区的问题？？？)
            elif(content1['country_id']!=content2['country_id']):
                duplicate_check9=1   #wrong country code
            elif(content1['dataset_id']!=content2['dataset_id']):
                duplicate_check10=1  #wrong instrument type


    #### 第7个检查：在第6个检查的基础之上，判断是否完全是相同的，所有的信息都是相同的，就算是空的字符串也算是相同
    if(duplicate_check2==True):
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
            duplicate_check11=1

    #### 第8个检查：如果深度层数不相等，则判断数据是否被插值过（例如缺少表面值、缺少最后一个点的值、缺少中间几个数的值） 插值使用方法：线性插值
    if(content1['depth_number']!=content2['depth_number']):
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
        duplicate_check12=percent_t>0.75


    isDuplicate=duplicate_check1 or duplicate_check2 or duplicate_check3 or duplicate_check4 or duplicate_check5 or duplicate_check12
    # print(isDuplicate)
    # print(duplicate_check1, duplicate_check2,duplicate_check3,duplicate_check4,duplicate_check5)
    if(isDuplicate !=1 and isDuplicate!=2):
        duplicate_multimodels=(duplicate_check1,duplicate_check3,duplicate_check4,duplicate_check5,duplicate_check6,duplicate_check7,duplicate_check8,duplicate_check9,duplicate_check10,duplicate_check11,duplicate_check12)
        # return isDuplicate,duplicate_multimodels
        return isDuplicate,duplicate_multimodels

    if(1):    
        print('\n')
        print("走航观测相同时间不同位置: %d" %(duplicate_check1))
        print("准确重复或近似重复: %d" % (duplicate_check2))
        print("深度/温度/盐度被放缩: %d" % (duplicate_check3))
        print("走航观测同一时刻同一位置（具体到分钟）多个不同观测: %d" % (duplicate_check4))
        print("温度/盐度观测记录被四舍五入/截断: %d" % (duplicate_check5))
        print("错误地理位置信息：%d" % (duplicate_check6)) 
        print("错误日期信息：%d" % (duplicate_check7))
        print("错误时间信息：%d" % (duplicate_check8))
        print("错误国家代码（不同国家重复提交）：%d" % (duplicate_check9))
        print("错误仪器类型：%d" % (duplicate_check10))
        print("所有信息及数据完全相等：%d" % (duplicate_check11))
        print("数据被插值（部分数据缺失）：%d" % (duplicate_check12))

    if(isDuplicate==True):
        if(duplicate_check1==1 or duplicate_check2==1 or duplicate_check3==1 or duplicate_check4==1 or duplicate_check5==1):
            isDuplicate=1
        else:
            isDuplicate=2

    duplicate_multimodels=[duplicate_check1,duplicate_check3,duplicate_check4,duplicate_check5,duplicate_check6,duplicate_check7,duplicate_check8,duplicate_check9,duplicate_check10,duplicate_check11,duplicate_check12]
    # print(duplicate_multimodels)
    # return isDuplicate,duplicate_multimodels   ### only for N04_check_nc_duplicate_list4.py
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

def check_Mode(array1,array2,thresholds=0.9):
    levels=len(array1)
    #除法
    # divide=array1/array2
    # divide[np.abs(divide-1)<1e-5]=np.nan
    # print(divide)
    vals,counts=np.unique(array1/array2, return_counts=True)
    counts=counts[np.logical_and(vals==1.,np.isnan(vals))]
    isMode1=np.any((counts/levels)>thresholds)

    #减法
    # differene=array1-array2
    # differene[np.abs(differene-0.)<1e-5]=np.nan
    # print(differene)
    vals,counts=np.unique(array1-array2, return_counts=True)
    counts=counts[np.logical_and(vals==0.,np.isnan(vals))]
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




