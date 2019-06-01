import urllib2
import sys
import ssl
import re

context = ssl._create_unverified_context()

#print "GGMU!"
snpfile = sys.argv[1]
#print snpfile
SNP = open(snpfile, "r")
snplist = SNP.readlines()
SNP.close()
snp_pattern = re.compile(r'<span>chr[0-9XY]+:\d+\s*</span>')
gene_pattern = re.compile(r'/gene/\d+')

Smap = open("snp_map.txt", "w")
Gmap = open("gene_map.txt", "w")

website = "http://www.ncbi.nlm.nih.gov"
ua_headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.99 Safari/537.36"}
#print snplist
for snp in snplist:
    snp = snp[0:-1]
    if snp[0:2]!='rs':
        continue
    http = website+"/snp/"+snp
    #print http
    request = urllib2.Request(http,headers=ua_headers)
    try:
        response = urllib2.urlopen(request, context = context)
    except:
        print "wrong snp ID!"
        continue
    html = response.read()
    if (snp_pattern.search(html)):
        #print html[re.search(r'<span>chr[0-9XY]+:\d+\s*</span>', html).span()[0]:re.search(r'<span>chr[0-9XY]+:\d+\s*</span>',html).span()[1]]
        tmp = html[re.search(r'<span>chr[0-9XY]+:\d+\s*</span>', html).span()[0]:re.search(r'<span>chr[0-9XY]+:\d+\s*</span>',html).span()[1]]
        print snp+"\t"+tmp[re.search(r'chr[0-9XY]+:\d+\s', tmp).span()[0]:re.search(r'chr[0-9XY]+:\d+\s', tmp).span()[1]]
        Smap.write(snp+"\t"+tmp[re.search(r'chr[0-9XY]+:\d+\s', tmp).span()[0]:re.search(r'chr[0-9XY]+:\d+\s', tmp).span()[1]]+"\n")
    else:
        print "no snp"
    if (gene_pattern.search(html)):
        geneid = html[re.search(r'/gene/\d+',html).span()[0]:re.search(r'/gene/\d+',html).span()[1]]
        #print geneid
        try:
            request = urllib2.Request(website+geneid, headers = ua_headers)
            response = urllib2.urlopen(request, context = context)
            html = response.read()
            #print html[1:10]
            gene_symbol_pattern = re.compile(r'noline\W{2}[A-Z0-9]+<span')
            chr_pattern = re.compile(r'<td>[0-9XY]+</td>')
            position_pattern = re.compile(r'<td>.+\(\d+\.\.\d+.*\)\s*</td>')
            if (chr_pattern.search(html) and position_pattern.search(html) and gene_symbol_pattern.search(html)):
                #print html[re.search(r'<td>[0-9XY]+</td>', html).span()[0]:re.search(r'<td>[0-9XY]+</td>', html).span()[1]]
                tmp = html[re.search(r'<td>[0-9XY]+</td>', html).span()[0]:re.search(r'<td>[0-9XY]+</td>', html).span()[1]]
                chr = tmp[re.search(r'[0-9XY]+', tmp).span()[0]:re.search(r'[0-9XY]+', tmp).span()[1]]
                #print html[re.search(r'<td>.+\(\d+\.\.\d+.*\)\s*</td>', html).span()[0]:re.search(r'<td>.+\(\d+\.\.\d+.*\)\s*</td>', html).span()[1]]
                tmp = html[re.search(r'<td>.+\(\d+\.\.\d+.*\)\s*</td>', html).span()[0]:re.search(r'<td>.+\(\d+\.\.\d+.*\)\s*</td>', html).span()[1]]
                pos = tmp[re.search(r'\d+\.\.\d+', tmp).span()[0]:re.search(r'\d+\.\.\d+', tmp).span()[1]]
                #print html[re.search(r'noline\W{2}[A-Z0-9]+<span', html).span()[0]:re.search(r'noline\W{2}[A-Z0-9]+<span', html).span()[1]]
                tmp = html[re.search(r'noline\W{2}[A-Z0-9]+<span', html).span()[0]:re.search(r'noline\W{2}[A-Z0-9]+<span', html).span()[1]]
                symbol = tmp[re.search(r'[A-Z0-9]+', tmp).span()[0]:re.search(r'[A-Z0-9]+', tmp).span()[1]]
                print symbol+"\t"+chr+"\t"+pos
                Gmap.write(symbol+"\t"+chr+"\t"+pos+"\n")
            else:
                print "no gene position infomation!"
        except:
            print "wrong gene ID"
            continue
    else :
        print "no gene"
    #print html

Smap.close()
Gmap.close()
