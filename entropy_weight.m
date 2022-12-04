
function [weight]=entropy_weight(DNA_series)

%标准化操作mapstd  将数据标准化为均值为0，方差为1的数据    mapminmax归一化
%归一化和标准化，都可以探讨
% DNA_mapped=mapstd(DNA_series,0,1);
DNA_series(1,4)=1994;
DNA_mapped=mapminmax(DNA_series',0,1);
DNA_mapped=DNA_mapped';
DNA_mapped(DNA_series==0)=0;


B=DNA_mapped;
%B矩阵中的每一列均已经归一化至0-1
B(B==0)=0.00001;
B(B==1)=0.99999;
 
[n,m]=size(B); % n个样本, m个指标
%%计算第j个指标下，第i个样本占该指标的比重p(i,j)
p=NaN(n,m);
dd=sum(B,1,'omitnan');
for j=1:m
    p(:,j)=B(:,j)./dd(j);
end
% for i=1:n
%     i
%     for j=1:m
%         p(i,j)=B(i,j)/sum(B(:,j),'omitnan');
%     end
% end
%%计算第j个指标的熵值e(j)
k=1/log(n);
% e=NaN();
for j=1:m
%     j
    e(j)=-k*sum(p(:,j).*log(p(:,j)),'omitnan');
end
d=ones(1,m)-e; %计算信息熵冗余度
weight=d./sum(d,'omitnan'); %求权值w 主要需要的结果是这个
s=100*weight*B'; %求综合得分

end