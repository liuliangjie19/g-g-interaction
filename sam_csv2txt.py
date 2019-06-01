#!/usr/bin/python
import os
import sys
#print "Hala Madrid!"
def sam_csv2txt(file_name, outfilename):
    fo1 = open(file_name, "r")
    fo2 = open(outfilename, "w")
    sample_line = fo1.readline()
    i=0
    j=0
    while(sample_line):
        sample_line = sample_line[0:-1]
        line_list = sample_line.split(",")
        sample_line = fo1.readline()
        if (line_list[0]=='ID'):
            fo2.write("NO.\tPlate\tWell_NO.\tSample\tStat\n")
            continue
        i = int(line_list[0])
        j = int(line_list[1])
        if (line_list[10]=='SZ'):
            fo2.write(str(i)+"\t"+str(j)+"\t"+str(i-(j-1)*384)+"\t"+line_list[9]+"\t"+line_list[10]+"\n")
        else:
            fo2.write(str(i)+"\t"+str(j)+"\t"+str(i-(j-1)*384)+"\t"+line_list[9]+"\t"+"Healthy\n")
        #break
    fo1.close()
    fo2.close()

def info_check(sam, snp, plate):
    file1 = open(sam, "r")
    file2 = open(snp, "r")
    sam_line = file1.readline()
    snp_line = file2.readline()
    sam_dic = {}
    snp_dic = {}
    while(sam_line):
        sam_list = sam_line.split("\t")
        #print len(sam_list)
        if (sam_list[0]!='NO.' and sam_list[1]==plate):
            sam_dic[int(sam_list[2])] = sam_list[3]
        sam_line = file1.readline()
    while(snp_line):
        snp_list = snp_line.split("\t")
        #print len(snp_list)
        if (len(snp_list)==10 and snp_list[0]!='Well'):
            #print snp_line
            snp_dic[int(snp_list[0])] = snp_list[1]
        snp_line = file2.readline()
    for k in sam_dic.keys():
        if (sam_dic[k]!=snp_dic[k]):
            print k,sam_dic[k],snp_dic[k]
    file1.close()
    file2.close()

sample_file_name = sys.argv[1]
out_file_name = sys.argv[2]
#print sample_file_name,out_file_name
#FILE = open(sample_file_name, "r")
sam_csv2txt(sample_file_name, out_file_name)
#info_check(sample_file_name, out_file_name, 4)
