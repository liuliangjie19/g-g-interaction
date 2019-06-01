#!/usr/bin/python
import sys
import os
import numpy as np

#print "GGMU"
data_dir = sys.argv[1]
snp_info = sys.argv[2]
sample_info = sys.argv[3]
final_output = sys.argv[4]

if (data_dir[-1]!='/'):
    data_dir = data_dir+"/"
dir = os.popen("ls "+data_dir+"*.txt")
dir_list = dir.readlines()
dir.close()
#print data_dir, snp_info, final_output, len(dir_list)
snps = open(snp_info, "r")
snp_list = snps.readlines()
snps.close()
sams = open(sample_info, "r")
sam_list = sams.readlines()
sams.close()
final_data = np.empty([len(sam_list)-1,len(snp_list)+1])
rowname = []
rowname.append("SampleID")
rowname.append("Stat")
#print len(snp_list),len(sam_list),final_data[0,1]
for item in snp_list:
    item = item[0:len(item)-1]
    linelist = item.split("\t")
    rowname.append([linelist[5],linelist[6]])
print rowname
for item in sam_list:
    item = item[0:len(item)-1]
    linelist = item.split("\t")
#    if (linelist[0]=='ID'):
#        continue
