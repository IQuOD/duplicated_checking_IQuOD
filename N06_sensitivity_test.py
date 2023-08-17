### This program is used for sensitivity test, and is used to compare the difference before and after the modified program
f1=open('Filename_potential_duplicate_ALL_2011_unique_1119_new.txt','r')
f2=open('Filename_potential_duplicate_ALL_2011_unique_1119_new2.txt','r')

txt1=f1.read()
txt2=f2.read()

f1.close()
f2.close()

line1=txt1.split()
line2=txt2.split()

#outfile=open('sensitivity_test1975.txt','w')
#outfile=open('sensitivity_test1995.txt','w')
outfile=open('sensitivity_test2011.txt','w')


for i in line1:
    # Check whether the file in file1 exists in file2
    if i not in line2:
        outfile.write(i)
        outfile.write('\n')
outfile.write('Above content in file1,but not in file2')
outfile.write('\n')
outfile.write('\n')
for j in line2:
    # Check whether the file in file2 exists in file1
    if j not in line1:
        outfile.write(j)
        outfile.write('\n')
outfile.write('Above content in file2,but not in file1')
outfile.write('\n')

print('Finish!')



