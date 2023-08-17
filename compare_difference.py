###比较按行归一化和按列归一化结果的不同
# f1=open('Filename_potential_duplicate_ALL_1995_unique_1119.txt','r')
# f2=open('Filename_potential_duplicate_ALL_1995_unique_1119_new.txt','r')

f1=open('Filename_potential_duplicate_ALL_2011_unique_1119.txt','r')
f2=open('Filename_potential_duplicate_ALL_2011_unique_1119 _new.txt','r')

txt1=f1.read()
txt2=f2.read()

f1.close()
f2.close()

line1=txt1.split()
line2=txt2.split()

outfile=open('different_f1f2.txt','w')

for i in line1:
    #查看1中文件是否在2中存在
    if i not in line2:
        outfile.write(i)
        outfile.write('\n')
outfile.write('Above content in file1,but not in file2')
outfile.write('\n')
outfile.write('\n')
for j in line2:
    #查看1中文件是否在2中存在
    if j not in line1:
        outfile.write(j)
        outfile.write('\n')
outfile.write('Above content in file2,but not in file1')
outfile.write('\n')

print('核对结束！')



