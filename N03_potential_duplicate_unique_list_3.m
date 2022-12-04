%%%%%%%% duplicate list 去重取唯一
%%%%%%%% N02步骤中尝试搜索了很多个潜在的重复对文件，现在需要把这些文件里面的文件对进行合并（多个txt文件合并成一个大的txt文件）
%%%% 以便在下一步duplicate 检查中进行处理；但是合并的时候可能会有重复的文件对，这个时候需要“去重”，保证文件名的唯一性
clear
clc

files_txt=dir('./potential_duplicates_output/2008/potential_dup*.txt');

filenames_column1={};
filenames_column2={};
filenames_combines={};
m=1;
for i=1:length(files_txt)
    
    file1=[files_txt(i).folder,'/',files_txt(i).name];
    fid=fopen(file1,'r');
    
    while ~feof(fid)
        str=fgetl(fid);
        str=strtrim(str);
        s=regexp(str,'\s+','split');
        for k=2:length(s)
            filename1=s{1};
            filename2=s{k};
            
            filename_combine=[filename1,' ',filename2];
            
            filenames_column1{m}=filename1;
            filenames_column2{m}=filename2;
            filenames_combines{m}=filename_combine;
            m=m+1;
        end
    end
end

%%%%去重
filename_unique_pairs=unique(filenames_combines);


%%%%%% 重新输出新的potential_list
fid=fopen('./potential_duplicate_ALL_2008_unique_1119.txt','w+');
%%%%%输出原始数据文件
for m=1:length(filename_unique_pairs)
    s=regexp(filename_unique_pairs{m},'\s+','split');
    output_filename1=s{1};
    output_filename2=s{2};
    
    fprintf(fid,'%s   %s',output_filename1,output_filename2);
    fprintf(fid,'\n');
    clear s
end

length(filename_unique_pairs)
