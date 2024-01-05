
function [weight]=entropy_weight(DNA_series)

DNA_mapped=normalizeData(DNA_series,0,1);
DNA_mapped(DNA_series==0)=0;


B=DNA_mapped;
%%% Each column in the B matrix has been normalized to 0-1
B(B==0)=0.00001;
B(B==1)=0.99999;
 
[n,m]=size(B); % n samples, m variables
%%% Calculate the proportion p(i,j) of the ith sample in the jth variable.
p=NaN(n,m);
dd=sum(B,1,'omitnan');
for j=1:m
    p(:,j)=B(:,j)./dd(j);
end
%%% Calculate the entropy of the jth variable e(j)
k=1/log(n);
for j=1:m
    e(j)=-k*sum(p(:,j).*log(p(:,j)),'omitnan');
end
d=ones(1,m)-e; % calculate information entropy redundancy
weight=d./sum(d,'omitnan'); % calculate the weight w
s=100*weight*B'; % calculate composite score

end