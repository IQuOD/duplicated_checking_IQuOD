%%%%%%%%%% 将excel表格的每一条廓线元数据读取到matlab中来，以mat文件形式存储
clear
clc

for nian=2008:2008
    
    file=['../WOD_',num2str(nian),'.xlsx']
    
    [DNA_series,txt,~]=xlsread(file);  %一定不要用cell数组，要不然绝对卡死
    DNA_series=single(DNA_series);
    
    %%%整理文件名 排列好
    filename_info=char(txt(:,1));
    
    
    %%%字段名称
    % variable_name=raw(1,2:end);
    variable_name={'accessin_number','dataset_id','lat','lon','year','month','day','probe_type','recorder','hour','minute','depth_number','maximum_depth','hasTemp','hasSal','hasOxygen','hasChlonophyll','country_id','GMT_time','WMO_id','dbase_orig','Project_name','plarform','vehicle','Institute','WOD_cruise_identifier','sum_temp','sum_salinity','sum_depth','std_depth','std_temp','std_salinity','corr_temp_depth','corr_sal_depth'};
    
    DNA_series(:,1)=[]; %删除WOD_unique_id
    
    txt(:,1)=[];  %删除文件名
    txt(:,2)=[]; %删除WOD_unique_id
    
    
    
    %%%%将字符串转换为ASCILL码，然后求和
    istxt=all(isnan(DNA_series));
    istxt(26)=1;  %26. WOD_cruise_identifier是字符串
    for i=1:length(istxt)
        i
        if(istxt(i))
            column=i;
            ascill_all=char(txt(:,i))+0;
            %把所有空格设置为nan
            ascill_all(ascill_all==' ')=NaN;
            %每一行字符串求和
            sum_ascill_all=sum(ascill_all,2,'omitnan');
            DNA_series(:,i)=sum_ascill_all;
        end
    end
    
    DNA_series(DNA_series==999)=NaN;  %%%%999 缺测值设置为Nan
    
    %%%%%这样就整理完了，把所有元数据信息，都用数字的形式存储到了一个大的数组中，并且每一行代表一个廓线，对应的文件名用filename_info存储，这样就可以一一对应文件名
    eval(['save DNA_summary_',num2str(nian),'.mat DNA_series variable_name filename_info'])
end


