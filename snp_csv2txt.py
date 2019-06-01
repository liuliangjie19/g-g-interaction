#!/usr/bin/python
#print "Hala Madrid!"
import sys
import os

def get_vaild_snp(csv_file, txt_file):
    csv = open(csv_file, "r")
    txt = open(txt_file, "w")
    csv_line = csv.readline()
    while(csv_line):
        linelist = csv_line.split(",")
        if (linelist[0]=='ID'):
            csv_line = csv_line.replace(",", "\t")
            txt.write(csv_line)
            csv_line = csv.readline()
            continue
        if (linelist[11]=='OK'):
            csv_line = csv_line.replace(",", "\t")
            txt.write(csv_line)
        csv_line = csv.readline()
    csv.close()
    txt.close()

snp_csv_file = sys.argv[1]
output_file = sys.argv[2]
#snp_data_dir = sys.argv[3]
get_vaild_snp(snp_csv_file, output_file)
#print snp_csv_file, output_file
