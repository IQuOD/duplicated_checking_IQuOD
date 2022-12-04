import netCDF4 as nc
import numpy as np
import math
import os
import country_table as t_country
warnings.filterwarnings('ignore')
from check_functions import compair_main

class duplicate_flag(object):
    def __init__(self):
	    self.path='/Users/zqtzt/Downloads/WOD_rawdata/'
	    # print(self.path)

	def run(self):
		while True:
			print('---------请输入两个重复的netCDF文件名称--------')
			file1=input('第一个nc文件名称是：').rstrip().lstrip()
            file2=input('第二个nc文件名称是：').rstrip().lstrip()

            index_str=file1.rfind('_')
            date1=file1[index_str-14:index_str-6]
            year1=date1[0:4]
            month1=date1[4:6]
            path1=self.path+'/'+year1+'/'+month1
            index_str=file2.rfind('_')

            date2=file2[index_str-14:index_str-6]
            year2=date2[0:4]
            month2=date2[4:6]
            path2=self.path+'/'+year2+'/'+month2
            filepath1=os.path.join(path1,file1)
            filepath2=os.path.join(path2,file2)

            ###读取第一个nc文件数据
            content1=self.read_nc_data(filepath1)   #content1是一个字典
            ###读取第二个文件数据
            content2=self.read_nc_data(filepath2)

            ######判断保留哪个
            [flag1,flag2]=compair_main.flag_duplicate(content1,content2)

            #####写进原本netCDF文件
            self.write_flag_netCDF(filepath1,flag1)
            self.write_flag_netCDF(filepath2,flag2)

def main():
	#create object
	df=duplicate_flag()
	df.run()



if __name__ == '__main__':
	main()