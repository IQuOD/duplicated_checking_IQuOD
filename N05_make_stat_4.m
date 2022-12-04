
clear
clc

total_profiles=142537;  %%%1975年的所有廓线数量
duplicate_profile=433;
duplicate_precentage=duplicate_profile/total_profiles*100


filename='./potential_duplicates_output/1975/DuplicateList_potential_duplicate_ALL_1975_unique_1117.xlsx';
[num,txt,~] = xlsread(filename,'A2:BM433');




% 4. 统计国家及其和数量 第41列
% 5. 统计platform code及其和数量 第49列
% 6. 统计 WOD_cruise_identifier1 类型及数量/占比 第53列
%%%%%%%%%% 7. 统计每种重复类型的数量 及其占总重复对数的比例
% 第5列： 走航观测相同时间不同位置_1  第6列：观测被放缩_2 第7列：走航观测同一时刻同一位置多个不同观测_3
% 第8列：观测被四舍五入_4  第9列：错误地理位置信息_5 第10列：错误日期信息6 第11列：错误时间信息7 第12列：错误国家代码8
% 第13列：错误仪器类型9  第14列：所有信息及数据完全相等10



%%% 1.统计accession number 每个类别有多少个，看看是否存在有严重问题（数量多）的access_number
access_num_all=num(:,16);
[unique(access_num_all');histc(access_num_all',unique(access_num_all'));histc(access_num_all',unique(access_num_all'))./length(lat_all)*100]

access_num_copy=access_num_all;

%%%% 2. 画地理位置的图 
lat_all=num(:,18);
lon_all=num(:,20);
plot_geo_location(lon_all,lat_all,'Location of potential duplicated pairs')
saveas(gcf,'./pics/locaiton_1975.png')

% 3. 统计每种仪器类型重复数据的个数  第16列，看看是不是有一些仪器是特别有
instrument_all=txt(:,16);
instrument_type_unique=unique(instrument_all);
for i=1:length(lat_all)
    
    for ii=1:length(instrument_type_unique)
        if(instrument_all{i}==instrument_type_unique{ii})
            instrument_all_order(i)=ii;
        end
    end
end
stat_instrument=[unique(instrument_all_order);histc(instrument_all_order,unique(instrument_all_order));histc(instrument_all_order,unique(instrument_all_order))./length(lat_all)*100]
instrument_type_unique

yan'z

country_all=num(:,39)';
stat_country=[unique(country_all);histc(country_all,unique(country_all));histc(country_all,unique(country_all))./length(lat_all)*100]


platform_code_all=txt(:,49);
platform_unique=unique(platform_code_all)
% for i=1:length(lat_all)
%     
%     for ii=1:length(platform_unique)
%         if(platform_code_all{i}==platform_unique{ii})
%             platform_unique_order(i)=ii;
%         end
%     end
% end
% unique(platform_code_all)
% stat_platform=[unique(platform_unique_order);histc(platform_unique_order,unique(platform_unique_order));histc(platform_unique_order,unique(platform_unique_order))./length(lat_all)*100]
% 


WOD_cruise_identifier_all=txt(:,53);
WOD_cruise_identifier_unique=unique(WOD_cruise_identifier_all)
% for i=1:length(lat_all)
%     
%     for ii=1:length(WOD_cruise_identifier_unique)
%         if(WOD_cruise_identifier_all{i}==WOD_cruise_identifier_unique{ii})
%             WOD_cruise_identifier_order(i)=ii;
%         end
%     end
% end
% unique(WOD_cruise_identifier_all)
% stat_cruise_id=[unique(WOD_cruise_identifier_order);histc(WOD_cruise_identifier_order,unique(WOD_cruise_identifier_order));histc(WOD_cruise_identifier_order,unique(WOD_cruise_identifier_order))./length(lat_all)*100]
% 

% 第5列： 走航观测相同时间不同位置_1  第6列：观测被放缩_2 第7列：走航观测同一时刻同一位置多个不同观测_3
% 第8列：观测被四舍五入_4  第9列：错误地理位置信息_5 第10列：错误日期信息6 第11列：错误时间信息7 第12列：错误国家代码8
% 第13列：错误仪器类型9  第14列：所有信息及数据完全相等10
duplicate1=num(:,3);
duplicate2=num(:,4);
duplicate3=num(:,5);
duplicate4=num(:,6);
duplicate5=num(:,7);
duplicate6=num(:,8);
duplicate7=num(:,9);
duplicate8=num(:,10);
duplicate9=num(:,11);
duplicate10=num(:,12);

sum(duplicate1~=0),sum(duplicate1~=0)./length(duplicate1)
sum(duplicate2~=0),sum(duplicate2~=0)./length(duplicate2)
sum(duplicate3~=0),sum(duplicate3~=0)./length(duplicate3)
sum(duplicate4~=0),sum(duplicate4~=0)./length(duplicate4)
sum(duplicate5~=0),sum(duplicate5~=0)./length(duplicate5)
sum(duplicate6~=0),sum(duplicate6~=0)./length(duplicate6)
sum(duplicate7~=0),sum(duplicate7~=0)./length(duplicate7)
sum(duplicate8~=0),sum(duplicate8~=0)./length(duplicate8)
sum(duplicate9~=0),sum(duplicate9~=0)./length(duplicate9)
sum(duplicate10~=0),sum(duplicate10~=0)./length(duplicate10)



%% 存储wod_unique_id  剔除全部id，随机选取一列id来剔除 以方便测试OHC
clear
clc

filename='./potential_duplicates_output/DuplicateList_potential_duplicate_ALL_1995_unique.xlsx';
[num,txt,~] = xlsread(filename,'A2:BL312');

all_unique_id=[num(:,1),num(:,2)];
all_unique_id=all_unique_id(:);
remove_unique_id=all_unique_id;

save ./wod_unique_id_forOHC/all_duplicated_unique_ID_1995.mat remove_unique_id
clear remove_unique_id

half_unique_id=num(:,1);
remove_unique_id=half_unique_id;
save ./wod_unique_id_forOHC/half_duplicated_unique_ID_1995.mat remove_unique_id


