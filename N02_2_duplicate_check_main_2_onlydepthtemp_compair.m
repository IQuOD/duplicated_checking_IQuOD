%%%%  进行每一行归一化操作，然后求算术平均，接着就是搜索比较哪个更接近
%%%% 重点搜索深度和和温度和相等的
clear
clc

for nian=2008:2008
    
    eval(['load DNA_summary_',num2str(nian),'.mat'])
    
    
    DNA_series_temp_depth=DNA_series(:,[27,29]); %sum_temp,sum_depth
    DNA_series_temp_depth(abs(DNA_series_temp_depth)>1e6)=NaN; %%%空缺值missing value
    % figure()
    % plot(DNA_series_temp_depth(:,2),'o')
    
    %%%%%%%%%%  求深度和和温度和的平均
    average_DNA=nanmean(DNA_series_temp_depth,2);
    
    %对average_DNA进行升序排序，方便后面搜索算法的建立
    [average_DNA,index]=sort(average_DNA);
    filename_info=filename_info(index,:);
    DNA_series=DNA_series(index,:);
    
    %%% 循环搜索
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
        difference=abs((number1-average_DNA)/number1*100);   %差异百分比
        difference(1:i-1)=NaN;
        duplicate_number=sum(difference<0.0001);   %阈值0.001%      阈值可以在之后设定调整
        if(duplicate_number>=2)
            %%%%疑似重复
            %        pause
            %找difference相差最小和相差为0的
            difference(i)=NaN;
            id=[i;find(difference==nanmin(difference))];
            DNA_series_small=DNA_series(id,:);
            %%%%%%%%%%%%%%%%%%%如果是浮标MRB数据，先跳过 检查xxxxxxx
            if(DNA_series(i,2)==7)
                continue
            end
            
            if(DNA_series(i,12)<=3)  %%%%%%%%浮标观测忽略:深度个数小于3
                continue
            end
            
            
            %%%%%%%%%来做一下一些排除和判断 （先判断相似片段有多少个）
            fragment_same_number=sum(abs(DNA_series_small(1,:)-DNA_series_small(2,:))<1e-5,'omitnan');
            if(fragment_same_number<26)  %没有相似的片段，就要跳过了  %%%%%%%这里可以分级搜索，准确重复是严格等于32 或者31  改成27可以找出绝大部分准确重复
                continue
            end
            
            %%%%排除深度和相差很大的
            for m=2:length(id)
                if(abs(DNA_series_small(1,29)-DNA_series_small(m,29))>1)
                    id(m)=NaN;
                end
            end
            id(isnan(id))=[];
            if(length(id)<=1)
                continue
            end
            
            
            %%%%如果是XBT CTD MBT BOT，且在一个月内相差正负5度，且是同一个probe  排除走航连续观测
            %%%%type,platform,vehicle,但是sum_temp,corr(temp,depth)不一样，则判断为同一条航线 同一个调查船/平台 上的多次观测
            if((DNA_series_small(1,2)==4 && DNA_series_small(2,2)==4) || (DNA_series_small(1,2)==2 && DNA_series_small(2,2)==2) || (DNA_series_small(1,2)==1 && DNA_series_small(2,2)==1) || (DNA_series_small(1,2)==3 && DNA_series_small(2,2)==3))
                index1=all(DNA_series_small(1,[5,6,8,23,24,26])==DNA_series_small(2,[5,6,8,23,24,26]));  %都要一样
                index2= abs(DNA_series_small(1,27)-DNA_series_small(2,27))>0.099; %sum_temp不相等
                index3= abs(DNA_series_small(1,33)-DNA_series_small(2,33))>0.001;  %cor_temp_depth
                index4=any(abs(DNA_series_small(1,[3,4])-DNA_series_small(2,[3,4]))<5) && any(abs(DNA_series_small(1,[3,4])-DNA_series_small(2,[3,4]))>1e-5);
                if(index1 && index2 && index3 && index4)
                    continue
                end
            end
            %%%%%排除定点/附近点长时间连续观测 只看MRB Bottle SUR
            if((DNA_series_small(1,2)==1 && DNA_series_small(2,2)==1) || (DNA_series_small(1,2)==7 && DNA_series_small(2,2)==7) || (DNA_series_small(1,2)==5 && DNA_series_small(2,2)==5))
                index1=all(DNA_series_small(1,[5,6,8,9,22,23,24])==DNA_series_small(2,[5,6,8,9,22,23,24]));  %都要一样
                index2=abs(DNA_series_small(1,27)-DNA_series_small(2,27))>0.05; %sum_temp不相等
                index3=abs(DNA_series_small(1,29)-DNA_series_small(2,29))<1e-5;  %sum_depth相等
                index4=all(abs(DNA_series_small(1,[3,4])-DNA_series_small(2,[3,4]))<0.01);  %定点：经纬度小于0.01°
                if(index1 && index2 && index3 && index4)
                    continue
                end
            end
            
            %%%%%%%%排除盐度和不相等的，且不为0
            if((abs(DNA_series_small(1,28)-DNA_series_small(2,28))>1e-3) && (DNA_series_small(1,28)>1e-6))
                continue
            end
            
            
            %         %%%%%输出原始数据文件
            %         fprintf(fid,'%s\n','【Potential Duplicates pairs】:');
            %         for m=1:length(id)
            %             fprintf(fid,'%s ',filename_info(id(m),:));
            %             fprintf(fid,'%3d %.4f %.4f %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %.4f %.4f %3d %.4f %.4f %.4f %.4f %.4f\n',DNA_series(id(m),:));
            %         end
            %         fprintf(fid,'\n');
            
            
            %%%%%输出原始数据文件
            for m=1:length(id)
                fprintf(fid,'%s ',filename_info(id(m),:));
            end
            fprintf(fid,'\n');
            
            %                pause
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