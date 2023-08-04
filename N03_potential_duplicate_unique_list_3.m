%%% In the N02 step, we  search for many potential duplicate files, now we need to merge the file pairs
%%% Merge multiple txt files into one large txt file to ensure the uniqueness of the file name
clear
clc

files_txt=dir('./potential_duplicates_output/1975/potential_dup*.txt');

filenames_column1={};
filenames_column2={};
filenames_combines={};
m=1;
for i=1:length(files_txt)
    
    file1=[files_txt(i).folder,'/',files_txt(i).name];
    fid=fopen(file1,'r');
    
    while ~feof(fid)
        str=fgetl(fid); % read a line in a file
        str=strtrim(str); % crop the Spaces at the beginning and end of the string
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

%%% Ensure the uniqueness of the file name
filename_unique_pairs=unique(filenames_combines);


%%% Reoutput potential_list
fid=fopen('./potential_duplicate_ALL_1992_unique_1119.txt','w+');
%%% Output filename
for m=1:length(filename_unique_pairs)
    s=regexp(filename_unique_pairs{m},'\s+','split');
    output_filename1=s{1};
    output_filename2=s{2};
    
    fprintf(fid,'%s   %s',output_filename1,output_filename2);
    fprintf(fid,'\n');
    clear s
end

length(filename_unique_pairs)
