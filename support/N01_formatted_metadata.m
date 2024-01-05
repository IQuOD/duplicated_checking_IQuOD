%%%%%%%%%% Read each profile metadata from excel table and stored in mat file form
clear
clc

for nian=1995:1995 %year
    
    file=['../Input_files/Metadata_summary_',num2str(nian),'.mat']
    load(file)

    DNA_series=single(DNA_series);
    txt(1,:)=[]; % Delete variable name
    
    %%% Organize and arrange file names
    filename_info=char(txt(:,1));
    
    
    %%% Variable name
    variable_name={'accessin_number','dataset_id','lat','lon','year','month','day','probe_type','recorder','hour','minute','depth_number','maximum_depth','hasTemp','hasSal','hasOxygen','hasChlonophyll','country_id','GMT_time','WMO_id','dbase_orig','Project_name','plarform','vehicle','Institute','WOD_cruise_identifier','sum_temp','sum_salinity','sum_depth','std_depth','std_temp','std_salinity','corr_temp_depth','corr_sal_depth'};
    
    DNA_series(:,1)=[]; % Delete WOD_unique_id
    
    txt(:,1)=[];  % Delete filename column
    txt(:,2)=[]; % Delete WOD_unique_id column
    
    
    
    %%%% Converts the string to the ASCILL code and sums
    istxt=all(isnan(DNA_series));
    istxt(26)=1;  % WOD_cruise_identifier is string
    for i=1:length(istxt)
        i
        if(istxt(i))
            column=i;
            ascill_all=char(txt(:,i))+0;
            %%% Set all NULL to nan
            ascill_all(ascill_all==' ')=NaN;
            %%% Sum per line
            sum_ascill_all=sum(ascill_all,2,'omitnan');
            DNA_series(:,i)=sum_ascill_all;
        end
    end
    
    DNA_series(DNA_series==999)=NaN;  %%%% 999 Set missing to Nan
    
    %%% All metadata information is stored in a large array in the form of numbers, each line represents a profile, and the corresponding filename is stored with filename_info to achieve a one-to-one correspondence between the filename and the data
    eval(['save ../Input_files/DNA_summary_',num2str(nian),'_20231208.mat DNA_series variable_name filename_info'])
end


