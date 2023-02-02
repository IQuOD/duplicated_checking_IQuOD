%%%利用主成分分析得到每条廓线的DNA
%%%主成分贡献率阈值设置为95%

clear
clc

for nian=1995:1995
    eval(['load DNA_summary_',num2str(nian),'.mat'])
    
    DNA_series_copy=DNA_series;
    DNA_series_copy(:,20)=[];   %第20列 WMO_ID信息量很少，去掉
    DNA_series_copy(isnan(DNA_series_copy))=0;  %缺省值设置为0
    
    %%%利用主成分分析计算每条廓线的DNA
    x=DNA_series_copy;
    z=zscore(x);    %数据按列标准化
    R=cov(z);       %协方差矩阵
    [V,D]=eig(R);   %计算协方差矩阵的特征向量和特征根
    d=diag(D);      
    eig1=sort(d,'descend');  %将主成分按从大到小排列
    v=fliplr(V);    
    s=0;
    i=0;
    while s/sum(eig1)<0.95   %主成分贡献率大于95%
        i=i+1;
        s=s+eig1(i);
    end
    DNA_new=z*v(:,1:i);
%     w=100*eig1/sum(eig1);
%     figure(1)
%     pareto(w)                 %绘制贡献率直方图
    

    [weight]=entropy_weight(DNA_new);
%     figure(); bar(weight)

    %%%  加权平均
    average_DNA_single=NaN(size(DNA_new));
    for i=1:length(weight)
        average_DNA_single(:,i)=DNA_new(:,i)*weight(i);
    end
    average_DNA=sum(average_DNA_single,2,'omitnan');
    
    %对average_DNA进行升序排序，方便后面搜索算法的建立
    [average_DNA,index]=sort(average_DNA);
    filename_info=filename_info(index,:);
    DNA_new=DNA_new(index,:);
    DNA_series=DNA_series(index,:);
    
    %%% 循环搜索
    output_variables=['filename',variable_name];
    
    filename=['./potential_duplicates_output/',num2str(nian),'/potential_duplicat_',num2str(nian),'_PCA_95_allinfo.txt']
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
        duplicate_number=sum(difference<0.001);   %阈值0.01%      
        if(duplicate_number>=2)
            difference(i)=NaN;
            id=[i;find(difference==nanmin(difference))];
            DNA_series_small=DNA_series(id,:);
            %%%如果是浮标SUR/MRB/DRB数据，先跳过 检查xxxxxxx
            if(DNA_series(i,2)==5||DNA_series(i,2)==7||DNA_series(i,2)==9) %分别代表SUR MRB DRB仪器
                continue
            end
            
            %%%针对深度或者温度或者盐度同时加上几个数的，或者乘上缩小放大一个倍数的，直接输出 相关系数肯定是相等的 %%%%其实这个可以去到另外一个模块
            %保留
            if(any(abs(DNA_series_small(1,[33,34])-DNA_series_small(2,[33,34]))<1e-4))
                %%%输出原始数据文件
                for m=1:length(id)
                    fprintf(fid,'%s ',filename_info(id(m),:));
                end
                fprintf(fid,'\n');
                
                number_pairs=number_pairs+1;
                number_profiles=number_profiles+duplicate_number;
                continue
            end
            
            %%%来做一下一些排除和判断 （先判断相似片段有多少个）
            fragment_same_number=sum(abs(DNA_series_small(1,:)-DNA_series_small(2,:))<1e-5,'omitnan');
            if(fragment_same_number<27)  %没有相似的片段，就要跳过了  %%%这里可以分级搜索，准确重复是严格等于33 或者32  改成27可以找出绝大部分准确重复
                continue
            end 
            
            %%%如果是XBT CTD MBT BOT，且在一个月内相差正负5度，且是同一个probe  排除走航连续观测
            %%%type,platform,vehicle,但是sum_temp,corr(temp,depth)不一样，则判断为同一条航线 同一个调查船/平台 上的多次观测
            if((DNA_series_small(1,2)==4 && DNA_series_small(2,2)==4) || (DNA_series_small(1,2)==2 && DNA_series_small(2,2)==2) || (DNA_series_small(1,2)==1 && DNA_series_small(2,2)==1) || (DNA_series_small(1,2)==3 && DNA_series_small(2,2)==3))
                index1=all(DNA_series_small(1,[5,6,8,23,24,26])==DNA_series_small(2,[5,6,8,23,24,26]));  %都要一样
                index2= abs(DNA_series_small(1,27)-DNA_series_small(2,27))>0.099; %sum_temp不相等
                index3= abs(DNA_series_small(1,33)-DNA_series_small(2,33))>0.001;  %cor_temp_depth
                index4=any(abs(DNA_series_small(1,[3,4])-DNA_series_small(2,[3,4]))<5) && any(abs(DNA_series_small(1,[3,4])-DNA_series_small(2,[3,4]))>1e-5);
                if(index1 && index2 && index3 && index4)
                    continue
                end
            end
            %%%排除定点/附近点长时间连续观测 只看MRB Bottle SUR
            if((DNA_series_small(1,2)==1 && DNA_series_small(2,2)==1) || (DNA_series_small(1,2)==7 && DNA_series_small(2,2)==7) || (DNA_series_small(1,2)==5 && DNA_series_small(2,2)==5))
                index1=all(DNA_series_small(1,[5,6,8,9,22,23,24])==DNA_series_small(2,[5,6,8,9,22,23,24]));  %都要一样
                index2=abs(DNA_series_small(1,27)-DNA_series_small(2,27))>0.05; %sum_temp不相等
                index3=abs(DNA_series_small(1,29)-DNA_series_small(2,29))<1e-5;  %sum_depth相等
                index4=all(abs(DNA_series_small(1,[3,4])-DNA_series_small(2,[3,4]))<0.01);  %定点：经纬度小于0.01°
                if(index1 && index2 && index3 && index4)
                    continue
                end
            end
            
            %%%输出原始数据文件
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
end