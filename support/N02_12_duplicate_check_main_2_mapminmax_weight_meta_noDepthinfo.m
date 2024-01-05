%%% Using mapminmax normalize each column of data, calculate the weighted average by entropy weight method, and then compare which is closer
%%% Ignore depth number and maximum depth information, the goal is to find duplicate pairs that may have been manipulated in depth number and maximum depth
clear
clc

load('DNA_summary_1975.mat')

DNA_series_meta=DNA_series(:,[1:26]);
DNA_series_meta(:,[12,13,20])=[]; % delete depth_number, maximum depth and WMO ID information

%%% Using mapminmax normalize each column of data
DNA_mapped=normalizeData(DNA_series_meta,0,1);
DNA_mapped(:,5)=0;
DNA_mapped(DNA_series_meta==0)=0;

%%% Use entropy weight method to calculate the weight
[weight]=entropy_weight(DNA_series_meta);
figure();bar(weight)

%%% Calculate the weighted average
average_DNA_single=NaN(size(DNA_mapped));
for i=1:length(weight)
    average_DNA_single(:,i)=DNA_mapped(:,i)*weight(i);
end
average_DNA=sum(average_DNA_single,2,'omitnan');

%%% Sort average_DNA in ascending order to facilitate the establishment of the later search algorithm
[average_DNA,index]=sort(average_DNA);
filename_info=filename_info(index,:);
DNA_mapped=DNA_mapped(index,:);
DNA_series=DNA_series(index,:);
DNA_series_meta=DNA_series_meta(index,:);

% figure();plot(average_DNA,'o')
%%% Cyclic search
output_variables=['filename',variable_name];
filename='./potential_duplicates_output/1975/potential_duplicate_1975_mapminmax_weight_meta_noDepthInfo.txt';
if(exist(filename))
    delete(filename)
end
fid=fopen(filename,'w+');

number_pairs=0;
number_profiles=0;
for i=1:length(average_DNA)
    i
    number1=average_DNA(i);
    difference=abs((number1-average_DNA)/number1*100);   % Calculation of percentage difference
    difference(1:i-1)=NaN;
    duplicate_number=sum(difference<0.0001);   % threshold value: 0.0001%
    if(duplicate_number>=2)
        %%% potiential duplicate
        difference(i)=NaN;
        id=[i;find(difference==nanmin(difference))];
        DNA_series_small=DNA_series(id,:);
        DNA_series_small_meta=DNA_series_meta(id,:);

        %%% Calculate how many similar fragments there are
        fragment_same_number=sum(abs(DNA_series_small_meta(1,:)-DNA_series_small_meta(2,:))<1e-5,'omitnan');
        if(DNA_series(i,2)==7 || DNA_series(i,2)==9 || DNA_series(i,2)==5)  % DRB MRB SUR
            if(fragment_same_number<23) % less than 23
                continue
            end
        else
            if(fragment_same_number<20)  % less than 20
                continue
            end
        end
       
        %%% Exclude long-term continuous observation of fixed points/nearby points(MRB¡¢Bottle¡¢SUR)
        if((DNA_series_small(1,2)==1 && DNA_series_small(2,2)==1) || (DNA_series_small(1,2)==7 && DNA_series_small(2,2)==7) || (DNA_series_small(1,2)==5 && DNA_series_small(2,2)==5))
            index1=all(DNA_series_small(1,[5,6,8,9,22,23,24])==DNA_series_small(2,[5,6,8,9,22,23,24]));  
            index2=   abs(DNA_series_small(1,27)-DNA_series_small(2,27))>0.05; % sum_temp is different
            index3= abs(DNA_series_small(1,29)-DNA_series_small(2,29))<1e-5;   % sum_depth is same
            index4=all(abs(DNA_series_small(1,[3,4])-DNA_series_small(2,[3,4]))<0.01);  % fixed point: latitude and longitude less than 0.01 degree
            if(index1 && index2 && index3 && index4)
                continue
            end
        end
        
        %%% Output filename
        for m=1:length(id)
            fprintf(fid,'%s ',filename_info(id(m),:));
        end
        fprintf(fid,'\n');
        
        number_pairs=number_pairs+1;
        number_profiles=number_profiles+duplicate_number;
        
    end
end

number_pairs
number_profiles


%%
% figure();
% plot(average_DNA2,'o');
% ylabel('Average DNA')


% for i=1:length(average_DNA)
%     if(contains(filename_info(i,:),'CASv1_T_S_19950724_00859_CTD.nc'))
%         i
%         pause
%     end
% end
