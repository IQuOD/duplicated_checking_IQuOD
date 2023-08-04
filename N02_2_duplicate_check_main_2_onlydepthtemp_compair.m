%%% Focus on searching for potential duplicate pairs that are equal in sum_depth and sum_temperature
clear
clc

for nian=1975:1975
    
    eval(['load DNA_summary_',num2str(nian),'.mat'])
    
    
    DNA_series_temp_depth=DNA_series(:,[27,29]); % sum_temp,sum_depth
    DNA_series_temp_depth(abs(DNA_series_temp_depth)>1e6)=NaN; % set missing value to Nan
    
    %%% Calculate the average of sum_depth and sum_temp
    average_DNA=nanmean(DNA_series_temp_depth,2);
 
    %%% Sort average_DNA in ascending order to facilitate the establishment of the later search algorithm
    [average_DNA,index]=sort(average_DNA);
    filename_info=filename_info(index,:);
    DNA_series=DNA_series(index,:);
    
    %%% Cyclic search
    output_variables=['filename',variable_name];
    
    filename=['./potential_duplicates_output/',num2str(nian),'/potential_duplicate_',num2str(nian),'_depth_temp.txt'];
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
            
            %%% If it's buoy data(MRB), skip
            if(DNA_series(i,2)==7) % MRB
                continue
            end
            
            if(DNA_series(i,12)<=3)  % depth number less than 3
                continue
            end
            
            
            %%% Calculate how many similar fragments there are
            fragment_same_number=sum(abs(DNA_series_small(1,:)-DNA_series_small(2,:))<1e-5,'omitnan');
            if(fragment_same_number<26)  % less than 26
                continue
            end
            
            %%% Exclude sum_depth vary greatly
            for m=2:length(id)
                if(abs(DNA_series_small(1,29)-DNA_series_small(m,29))>1)
                    id(m)=NaN;
                end
            end
            id(isnan(id))=[];
            if(length(id)<=1)
                continue
            end
            
            %%% If it is XBT CTD MBT BOT; the location difference is plus or minus 5 degrees within a month; the same probe----excludes navigation continuous observation
            %%% If type,platform, and vehicle are the same, but sum_temp,corr(temp,depth) are different, it is judged to be multiple observations on the same survey ship/platform on the same route
            if((DNA_series_small(1,2)==4 && DNA_series_small(2,2)==4) || (DNA_series_small(1,2)==2 && DNA_series_small(2,2)==2) || (DNA_series_small(1,2)==1 && DNA_series_small(2,2)==1) || (DNA_series_small(1,2)==3 && DNA_series_small(2,2)==3))
                index1=all(DNA_series_small(1,[5,6,8,23,24,26])==DNA_series_small(2,[5,6,8,23,24,26]));
                index2= abs(DNA_series_small(1,27)-DNA_series_small(2,27))>0.099; % sum_temp is different
                index3= abs(DNA_series_small(1,33)-DNA_series_small(2,33))>0.001; % cor_temp_depth is different
                index4=any(abs(DNA_series_small(1,[3,4])-DNA_series_small(2,[3,4]))<5) && any(abs(DNA_series_small(1,[3,4])-DNA_series_small(2,[3,4]))>1e-5);
                if(index1 && index2 && index3 && index4)
                    continue
                end
            end
            %%% Exclude long-term continuous observation of fixed points/nearby points(MRB¡¢Bottle¡¢SUR)
            if((DNA_series_small(1,2)==1 && DNA_series_small(2,2)==1) || (DNA_series_small(1,2)==7 && DNA_series_small(2,2)==7) || (DNA_series_small(1,2)==5 && DNA_series_small(2,2)==5))
                index1=all(DNA_series_small(1,[5,6,8,9,22,23,24])==DNA_series_small(2,[5,6,8,9,22,23,24]));  
                index2=abs(DNA_series_small(1,27)-DNA_series_small(2,27))>0.05; % sum_temp is different
                index3=abs(DNA_series_small(1,29)-DNA_series_small(2,29))<1e-5; % sum_depth is same
                index4=all(abs(DNA_series_small(1,[3,4])-DNA_series_small(2,[3,4]))<0.01);   % fixed point: latitude and longitude less than 0.01 degree
                if(index1 && index2 && index3 && index4)
                    continue
                end
            end
            
            %%% Exclude sum_salinity unequal and not 0
            if((abs(DNA_series_small(1,28)-DNA_series_small(2,28))>1e-3) && (DNA_series_small(1,28)>1e-6))
                continue
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
end
%%
% figure();
% plot(average_DNA2,'o');
% ylabel('Average DNA')