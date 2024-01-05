%%%%%%% read metadata and secondary data from WOD18 netCDF files
clear
clc

%%%%%% path with the WOD18 netCDf files 
filepath='../Input_files/WOD18_sample_1995/';  %%% This path doesn't include all the year 1995 netCDF files. It is justed a demo
%%%%%%% You need to download the netCDF files from WOD https://www.ncei.noaa.gov/access/world-ocean-database-select/dbsearch.html

filenames=dir(filepath);
filenames([filenames.isdir])=[];
n_prof=length(filenames);

DNA_series=single(NaN(n_prof,35));
txt=cell([n_prof+1,36]);

meta_name={'','wod_unique_id','accession_number','dataset_id','lat','lon','year','month','day','probe type','recorder','hour','minute','depth_number','maximum_depth','hasTemp','hasSalinity','hasOxygen','hasChlonophyll','country_name','GMT_time','WMO_ID','dbase_orig','project_name','Platform','ocean_vehicle',...
    'Institute','WOD_cruise_identifier','sum_temp','sum_salinity','sum_depth','std_depth','std_temp','std_salinity','corr_temp_depth','corr_sal_depth'};
for m=1:length(meta_name)
    txt{1,m}=meta_name{m};
end

for m=1:n_prof
    m
    file=[filepath,filenames(m).name];

    % Reading data from the NetCDF file
    try
        f = netcdf.open(file, 'NC_NOWRITE');
    catch
        continue
    end


    % Read attributes
    try
        probe_type = netcdf.getVar(f, netcdf.inqVarID(f, 'Temperature_Instrument'));
        probe_type=probe_type';
        nonEmptyid=find(probe_type==0,1)-1;
        probe_type=probe_type(1:nonEmptyid);
    catch
        probe_type = '';
    end

    try
        recorder = netcdf.getVar(f, netcdf.inqVarID(f, 'Recorder'));
        recorder=recorder';
        nonEmptyid=find(recorder==0,1)-1;
        recorder=recorder(1:nonEmptyid);
    catch
        recorder = '';
    end

    try
        need_z_fix = netcdf.getVar(f, netcdf.inqVarID(f, 'needs_z_fix'));
        need_z_fix=need_z_fix';
        nonEmptyid=find(need_z_fix==0,1)-1;
        need_z_fix=need_z_fix(1:nonEmptyid);
    catch
        need_z_fix = '';
    end

    % Read variables
    timeID = netcdf.inqVarID(f, 'time');
    time = netcdf.getVar(f, timeID);
    baseDate = datetime(1770, 1, 1,0,0,0);
    dtime = baseDate + days(time);
    xiaoshi = hour(dtime);
    fenzhong = minute(dtime);

    % Read 'date' and parse it
    date = num2str(netcdf.getVar(f, netcdf.inqVarID(f, 'date')));
    year = str2double(date(1:4));
    month = str2double(date(5:6));
    day = str2double(date(7:8));

    depthID = netcdf.inqVarID(f, 'z');
    depth = netcdf.getVar(f, depthID);
    depth(depth > 12000 | depth < -10) = NaN;

    depth_number = length(depth);
    maximum_depth = depth(end);
    sum_depth = round(nansum(depth), 4);
    std_depth = round(nanstd(depth), 4);
    if isnan(sum_depth)
        sum_depth = 999;
    end
    if isnan(std_depth)
        std_depth = 999;
    end

    % Read 'Temperature' and apply conditions
    try
        tempID = netcdf.inqVarID(f, 'Temperature');
        temp = netcdf.getVar(f, tempID);
        temp(temp > 40 | temp < -2.5) = NaN;

        temp2 = temp(~isnan(temp));
        depth2 = depth(~isnan(temp));

        hasTemp = 1;
        sum_temp = round(nansum(temp), 4);
        std_temp = round(nanstd(temp), 4);
        cor_temp_depth = round(corr(temp2, depth2), 5);

        if isnan(cor_temp_depth)
            cor_temp_depth = 999;
        end
        if isnan(sum_temp) || sum_temp == 0.0
            sum_temp = 999;
        end
        if isnan(std_temp) || std_temp == 0.0
            std_temp = 999;
        end
    catch
        hasTemp = 0;
        sum_temp = 999;
        std_temp = 999;
        cor_temp_depth = 999;
    end

    % Read 'salinity' and apply conditions
    try
        salID = netcdf.inqVarID(f, 'Salinity');
        sal = netcdf.getVar(f, salID);
        sal(sal > 43 | sal < -1) = NaN;

        sal2 = sal(~isnan(sal));
        depth2 = depth(~isnan(sal));

        hasSalinity = 1;
        sum_sal = round(nansum(temp), 4);
        std_sal = round(nanstd(temp), 4);
        cor_sal_depth = round(corr(sal2, depth2), 5);

        if isnan(cor_sal_depth)
            cor_sal_depth = 999;
        end
        if isnan(sum_sal) || sum_sal == 0.0
            sum_sal = 999;
        end
        if isnan(std_sal) || std_sal == 0.0
            std_sal = 999;
        end
    catch
        hasSalinity = 0;
        sum_sal = 999;
        std_sal = 999;
        cor_sal_depth = 999;
    end

    try
        oxygenID = netcdf.inqVarID(f, 'Oxygen');
        hasOxygen = 1;
    catch
        hasOxygen = 0;
    end
    try
        oxygenID = netcdf.inqVarID(f, 'Chlorophyll');
        hasChlonophyII = 1;
    catch
        hasChlonophyII = 0;
    end

    % Read 'wod_unique_cast'
    try
        wod_unique_id = double(netcdf.getVar(f, netcdf.inqVarID(f, 'wod_unique_cast')));
    catch
        wod_unique_id = 999;
    end

    try
        countryNameID = netcdf.inqVarID(f, 'country');
        country_name = char(netcdf.getVar(f, countryNameID));
        country_name=country_name';
        nonEmptyid=find(country_name==0,1)-1;
        country_name=country_name(1:nonEmptyid);
    catch
        country_name = '';
    end

    % Read 'GMT_time' variable
    try
        GMTTimeID = netcdf.inqVarID(f, 'GMT_time');
        GMT_time = round(double(netcdf.getVar(f, GMTTimeID)), 2);
    catch
        GMT_time = 0;
    end

    % Read 'WMO_ID' variable
    try
        WMOID = netcdf.inqVarID(f, 'WMO_ID');
        WMO_id = double(netcdf.getVar(f, WMOID));
    catch
        WMO_id = 999;
    end

    try
        dbase_orig = netcdf.getVar(f, netcdf.inqVarID(f, 'dbase_orig'));
        dbase_orig=dbase_orig';
        nonEmptyid=find(dbase_orig==0,1)-1;
        dbase_orig=dbase_orig(1:nonEmptyid);
    catch
        dbase_orig = '';
    end

    try
        project_name = netcdf.getVar(f, netcdf.inqVarID(f, 'Project'));
        project_name=project_name';
        nonEmptyid=find(project_name==0,1)-1;
        project_name=project_name(1:nonEmptyid);
    catch
        project_name = '';
    end

    try
        Platform = netcdf.getVar(f, netcdf.inqVarID(f, 'Platform'));
        Platform=Platform';
        nonEmptyid=find(Platform==0,1)-1;
        Platform=Platform(1:nonEmptyid);
    catch
        Platform = '';
    end

    try
        Ocean_Vehicle = netcdf.getVar(f, netcdf.inqVarID(f, 'Ocean_Vehicle'));
        Ocean_Vehicle=Ocean_Vehicle';
        nonEmptyid=find(Ocean_Vehicle==0,1)-1;
        Ocean_Vehicle=Ocean_Vehicle(1:nonEmptyid);
    catch
        Ocean_Vehicle = '';
    end

    try
        accession_number = netcdf.getVar(f, netcdf.inqVarID(f, 'Access_no'));
    catch
        accession_number = 999;
    end

    try
        Institute = netcdf.getVar(f, netcdf.inqVarID(f, 'Institute'));
        Institute=Institute';
        nonEmptyid=find(Institute==0,1)-1;
        Institute=Institute(1:nonEmptyid);
    catch
        Institute = '';
    end

    try
        WOD_cruise_identifier = netcdf.getVar(f, netcdf.inqVarID(f, 'WOD_cruise_identifier'));
        WOD_cruise_identifier=WOD_cruise_identifier';
        nonEmptyid=find(WOD_cruise_identifier==0,1)-1;
        WOD_cruise_identifier=WOD_cruise_identifier(1:nonEmptyid);
    catch
        WOD_cruise_identifier = '';
    end

    % Read 'dataset'
    dataset_name = char(netcdf.getVar(f, netcdf.inqVarID(f, 'dataset')));
    dataset_name=dataset_name';
    if contains(dataset_name, 'bod') || contains(dataset_name, 'bottle') || ...
            contains(dataset_name, 'ocean station') || contains(dataset_name, 'osd') || ...
            contains(dataset_name, 'low-resolution') || contains(dataset_name, 'low resolution')
        dataset_id = 1;
    elseif contains(dataset_name, 'towed') || contains(dataset_name, 'uor') || ...
            contains(dataset_name, 'undulating')
        dataset_id = 10;
    elseif contains(dataset_name, 'ctd') || contains(dataset_name, 'xctd')
        dataset_id = 2;
    elseif contains(dataset_name, 'mbt') || contains(dataset_name, 'mechanica') || ...
            contains(dataset_name, 'mb')
        dataset_id = 3;
    elseif contains(dataset_name, 'xbt') || contains(dataset_name, 'xb') || ...
            contains(dataset_name, 'expendable')
        dataset_id = 4;
    elseif contains(dataset_name, 'sur') || contains(dataset_name, 'surface')
        dataset_id = 5;
    elseif contains(dataset_name, 'apb') || contains(dataset_name, 'autonomous') || ...
            contains(dataset_name, 'animal')
        dataset_id = 6;
    elseif contains(dataset_name, 'mrb') || contains(dataset_name, 'moored') || ...
            contains(dataset_name, 'tao')
        dataset_id = 7;
    elseif contains(dataset_name, 'pfl') || contains(dataset_name, 'argo') || ...
            contains(dataset_name, 'profiling')
        dataset_id = 8;
    elseif contains(dataset_name, 'drb') || contains(dataset_name, 'drifting')
        dataset_id = 9;
    elseif contains(dataset_name, 'gld') || contains(dataset_name, 'glider')
        dataset_id = 11;
    elseif contains(dataset_name, 'dbt')
        dataset_id = 12;
    elseif contains(dataset_name, 'std')
        dataset_id = 13;
    elseif contains(dataset_name, 'microbt')
        dataset_id = 14;
    else
        dataset_id = 999;
    end


    % Read 'lat' and 'lon'
    latitude = round(double(netcdf.getVar(f, netcdf.inqVarID(f, 'lat'))), 4);
    longitude = round(double(netcdf.getVar(f, netcdf.inqVarID(f, 'lon'))), 4);

    % Get filename from the path
    [save_path, filename, ext] = fileparts(file);


    %%%%%%% put into the array and cell
    txt{m+1,1}=filenames(m).name;
    txt{m+1,10}=probe_type;
    txt{m+1,11}=recorder;
    txt{m+1,20}=country_name;
    txt{m+1,23}=dbase_orig;
    txt{m+1,24}=project_name;
    txt{m+1,25}=Platform;
    txt{m+1,26}=Ocean_Vehicle;
    txt{m+1,27}=Institute;
    txt{m+1,28}=WOD_cruise_identifier;

    DNA_series(m,1)=wod_unique_id;
    DNA_series(m,2)=accession_number;
    DNA_series(m,3)=dataset_id;
    DNA_series(m,4)=latitude;
    DNA_series(m,5)=longitude;
    DNA_series(m,6)=year;
    DNA_series(m,7)=month;
    DNA_series(m,8)=day;
    DNA_series(m,11)=xiaoshi;
    DNA_series(m,12)=fenzhong;
    DNA_series(m,13)=depth_number;
    DNA_series(m,14)=maximum_depth;
    DNA_series(m,15)=hasTemp;
    DNA_series(m,16)=hasSalinity;
    DNA_series(m,17)=hasOxygen;
    DNA_series(m,18)=hasChlonophyII;
    DNA_series(m,20)=GMT_time;
    DNA_series(m,21)=WMO_id;
    DNA_series(m,28)=sum_temp;
    DNA_series(m,29)=sum_sal;
    DNA_series(m,30)=sum_depth;
    DNA_series(m,31)=std_depth;
    DNA_series(m,32)=std_temp;
    DNA_series(m,33)=std_sal;
    DNA_series(m,34)=cor_temp_depth;
    DNA_series(m,35)=cor_sal_depth;

    netcdf.close(f);
end

save ../Input_files/Metadata_summary_1995.mat DNA_series txt
