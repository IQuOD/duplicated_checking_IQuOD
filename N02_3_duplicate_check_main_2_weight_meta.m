%%%%  进行每一行归一化操作，然后求算术平均，接着就是搜索比较哪个更接近
%%%% 这个有问题！太多bottle了，没有其他的仪器 只有一个CTD  很奇怪！！！！
%%检查程序！！！！！
%%我觉得这个可以把MRB加进来！！！

clear
clc

for nian=2008:2008
    eval(['load DNA_summary_',num2str(nian),'.mat'])
    
    DNA_series_meta=DNA_series(:,[1:26]);
    DNA_series_meta(:,20)=[]; %第20列 WMO_ID信息量很少
    %标准化操作mapstd  将数据标准化为均值为0，方差为1的数据    mapminmax归一化
    %归一化和标准化，都可以探讨
    % DNA_mapped=mapstd(DNA_series,0,1);
    DNA_mapped=mapminmax(DNA_series_meta,0,1);
    DNA_mapped(DNA_series_meta==0)=0;
    
    %%%%%%%%%%求权重
    [weight]=entropy_weight(DNA_series_meta);
    % weight=[weight(1:19),nan,weight(20:end)];
    
    %每一行求算术平均  -->后面可以加权平均   后面可以经度纬度权重加大
    % average_DNA=nanmean(DNA_mapped,2);
    %%%%%%%%%%  加权平均
    average_DNA_single=NaN(size(DNA_mapped));
    for i=1:length(weight)
        average_DNA_single(:,i)=DNA_mapped(:,i)*weight(i);
    end
    average_DNA=sum(average_DNA_single,2,'omitnan');
    
    
    %对average_DNA进行升序排序，方便后面搜索算法的建立
    [average_DNA,index]=sort(average_DNA);
    filename_info=filename_info(index,:);
    DNA_mapped=DNA_mapped(index,:);
    DNA_series=DNA_series(index,:);
    DNA_series_meta=DNA_series_meta(index,:);
    
    %%% 循环搜索
    output_variables=['filename',variable_name];
    
    filename=['./potential_duplicates_output/',num2str(nian),'/potential_duplicate_',num2str(nian),'_weight_meta.txt'];
    if(exist(filename))
        delete(filename)
    end
    fid=fopen(filename,'w+');
    % for i=1:length(output_variables)
    %     fprintf(fid,'%s ',output_variables{i})
    % end
    % output_filename='potential_duplicates.xlsx';
    
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
            DNA_series_small_meta=DNA_series_meta(id,:);
            %%%%%%%%%%%%%%%%%%%如果是浮标MRB数据，先跳过 检查xxxxxxx
            %         if(DNA_series(i,1)==7)
            %             continue
            %         end
            
            
            %%%%%%%%%来做一下一些排除和判断 （先判断相似片段有多少个）
            %没有相似的片段，就要跳过了
            fragment_same_number=sum(abs(DNA_series_small_meta(1,:)-DNA_series_small_meta(2,:))<1e-5,'omitnan');
            if(DNA_series(i,2)==7 || DNA_series(i,2)==9 || DNA_series(i,2)==5)  %浮标限制窄一些 DRB MRB SUR
                if(fragment_same_number<23)
                    continue
                end
            else
                if(fragment_same_number<22)
                    continue
                end
            end
            
            %%%%如果是XBT CTD MBT BOT，且在一个月内相差正负5度，且是同一个probe  排除走航连续观测
            %%%%type,platform,vehicle,但是sum_temp,corr(temp,depth)不一样，则判断为同一条航线 同一个调查船/平台 上的多次观测
            if((DNA_series_small(1,2)==4 && DNA_series_small(2,2)==4) || (DNA_series_small(1,2)==2 && DNA_series_small(2,2)==2) || (DNA_series_small(1,2)==1 && DNA_series_small(2,2)==1) || (DNA_series_small(1,2)==3 && DNA_series_small(2,2)==3))
                index1=all(DNA_series_small(1,[5,6,8,23,24,26])==DNA_series_small(2,[5,6,8,23,24,26]));  %都要一样
                index2=abs(DNA_series_small(1,27)-DNA_series_small(2,27))>0.099;
                index3= abs(DNA_series_small(1,33)-DNA_series_small(2,33))>0.001;  %cor_temp_depth
                index4=any(abs(DNA_series_small(1,[3,4])-DNA_series_small(2,[3,4]))<5) && any(abs(DNA_series_small(1,[3,4])-DNA_series_small(2,[3,4]))>1e-5);
                if(index1 && index2 && index3 && index4)
                    continue
                end
            end
            %%%%%排除定点/附近点长时间连续观测 只看MRB Bottle SUR
            if((DNA_series_small(1,2)==1 && DNA_series_small(2,2)==1) || (DNA_series_small(1,2)==7 && DNA_series_small(2,2)==7) || (DNA_series_small(1,2)==5 && DNA_series_small(2,2)==5))
                index1=all(DNA_series_small(1,[5,6,8,9,22,23,24])==DNA_series_small(2,[5,6,8,9,22,23,24]));  %都要一样
                index2=   abs(DNA_series_small(1,27)-DNA_series_small(2,27))>0.05; %sum_temp不相等
                index3= abs(DNA_series_small(1,29)-DNA_series_small(2,29))<1e-5;  %sun_depth相等
                index4=all(abs(DNA_series_small(1,[3,4])-DNA_series_small(2,[3,4]))<0.01);  %定点：经纬度小于0.01°
                if(index1 && index2 && index3 && index4)
                    continue
                end
            end
            
            
            
            %%%%%输出原始数据文件
            for m=1:length(id)
                fprintf(fid,'%s ',filename_info(id(m),:));
            end
            fprintf(fid,'\n');
            
            %%%%%输出原始数据文件
            %         fprintf(fid,'%s\n','【Potential Duplicates pairs】:');
            %         for m=1:length(id)
            %             fprintf(fid,'%s ',filename_info(id(m),:));
            %             fprintf(fid,'%3d %.4f %.4f %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %.4f %.4f %3d %.4f %.4f %.4f %.4f %.4f\n',DNA_series(id(m),:));
            %         end
            %         fprintf(fid,'\n');
            %                pause
            number_pairs=number_pairs+1;
            number_profiles=number_profiles+duplicate_number;
            
        end
    end
    
    number_pairs
    number_profiles
end


%%
% figure();
% plot(average_DNA2,'o');
% ylabel('Average DNA')
