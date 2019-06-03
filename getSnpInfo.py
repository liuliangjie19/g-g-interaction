import sys
import urllib2
import ssl
import re
from bs4 import BeautifulSoup

context = ssl._create_unverified_context()

#print "GGMU!"
snpfile = sys.argv[1]
#print snpfile
SNP = open(snpfile, "r")
result = open("snpGeneinfo.csv", "w")
snplist = SNP.readlines()
SNP.close()
snplist.pop(0)
snp_pattern = re.compile(r'<span>chr[0-9XY]+:\d+\s*</span>')
gene_pattern = re.compile(r'/gene/\d+')
website = "http://www.ncbi.nlm.nih.gov"
ua_headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.99 Safari/537.36"}

for snp in snplist:
    snp = snp[0:-1]
    print(snp+" start")
    #print(snp)
    http = website+"/snp/"+snp
    request = urllib2.Request(http,headers=ua_headers)
    try:
        response = urllib2.urlopen(request, context = context)
    except:
        print "wrong snp ID!"
        continue
    html = response.read()
    info = BeautifulSoup(html,features="lxml")
    #print info.prettify()
    sum = info.select('div[class="summary-box usa-grid-full"]')[0]
    sum1 = sum.select('dl[class="usa-width-one-half"]')[0]
    sum2 = sum.select('dl[class="usa-width-one-half"]')[1]
    spanlist = sum1.select('span')
    #pos_pattern = re.compile(r'chr[0-9XY]+:\d+')
    #all_pattern = re.compile(r'[ATCG]>[ATGC]')
    ddlist = sum1.select('dd')
    all = ddlist[2].string
    frelist = ddlist[-1].select('div')
    re_1000G_pattern = re.compile(r'1000G\)')
    for fre in frelist:
        try:
            if (re_1000G_pattern.search(fre.string)):
                #print(snp+fre.string)
                frequency = fre.string
        except:
            continue
    frelist = ddlist[-1].select('span')
    if (re_1000G_pattern.search(frelist[0].string)):
        #print(snp+frelist[0].string)
        frequency = frelist[0].string
    pos = re.sub(r'chr', "", spanlist[0].string)
    pos = re.sub(r':', "\t", pos)
    #print(snp+"\t"+pos+"\t"+re.sub(r'\s', "", all[re.search(r'[ATCG][ATCG\s/>]*[ATCG]',all).span()[0]:re.search(r'[ATCG][>ATCG\s/]*[ATCG]',all).span()[1]])+"\t"+re.sub(r'\s',"",frequency))
    divlist = sum2.select('dd')
    try:
        gene = divlist[1].select('div')[0].string
    except:
        gene = divlist[1].select('span')[0].string
    gene = re.sub(r'\s', "", gene)
    gene = re.sub(r':.+',"", gene)
    genelinklist = info.find_all("a", string = gene)
    #genelink = genelinklist[0]['href']
    try:
        genelink = genelinklist[0]['href']
        request = urllib2.Request(website+genelink, headers = ua_headers)
        response = urllib2.urlopen(request, context = context)
        html = response.read()
        info = BeautifulSoup(html,features="lxml")
        table = info.select('table[class="jig-ncbigrid"]')[0]
        genechr = table.select('td')[3].string
        generegion = table.select('td')[4].string
        generegion = generegion[re.search(r'\d+\.\.\d+', generegion).span()[0]:re.search(r'\d+\.\.\d+', generegion).span()[1]]
        generegion = re.sub(r'\.\.', "\t", generegion)
        #print genechr+generegion
        #print(snp+"\t"+pos+"\t"+re.sub(r'\s', "", all[re.search(r'[ATCG][ATCG\s/>]*[ATCG]',all).span()[0]:re.search(r'[ATCG][>ATCG\s/]*[ATCG]',all).span()[1]])+"\t"+re.sub(r'\s',"",frequency)+"\t"+gene+"\t"+genechr+"\t"+generegion)
        result.write(snp+"\t"+pos+"\t"+re.sub(r'\s', "", all[re.search(r'[ATCG][ATCG\s/>]*[ATCG]',all).span()[0]:re.search(r'[ATCG][>ATCG\s/]*[ATCG]',all).span()[1]])+"\t"+re.sub(r'\s',"",frequency)+"\t"+gene+"\t"+genechr+"\t"+generegion+"\n")
    except:
        #print(snp+"\t"+pos+"\t"+re.sub(r'\s', "", all[re.search(r'[ATCG][ATCG\s/>]*[ATCG]',all).span()[0]:re.search(r'[ATCG][>ATCG\s/]*[ATCG]',all).span()[1]])+"\t"+re.sub(r'\s',"",frequency)+"\t"+"no gene")
        result.write(snp+"\t"+pos+"\t"+re.sub(r'\s', "", all[re.search(r'[ATCG][ATCG\s/>]*[ATCG]',all).span()[0]:re.search(r'[ATCG][>ATCG\s/]*[ATCG]',all).span()[1]])+"\t"+re.sub(r'\s',"",frequency)+"\t"+"no gene"+"\n")
    print(snp+" End")
    #print(snp+"\t"+pos+"\t"+re.sub(r'\s', "", all[re.search(r'[ATCG][ATCG\s/>]*[ATCG]',all).span()[0]:re.search(r'[ATCG][>ATCG\s/]*[ATCG]',all).span()[1]])+"\t"+re.sub(r'\s',"",frequency)+"\t"+gene+"\t"+genechr+"\t"+generegion)
    #print(spanlist)
    #if (snp_pattern.search(html)):
        #print html[re.search(r'<span>chr[0-9XY]+:\d+\s*</span>', html).span()[0]:re.search(r'<span>chr[0-9XY]+:\d+\s*</span>',html).span()[1]]
    #    tmp = html[re.search(r'<span>chr[0-9XY]+:\d+\s*</span>', html).span()[0]:re.search(r'<span>chr[0-9XY]+:\d+\s*</span>',html).span()[1]]
    #    print snp+"\t"+tmp[re.search(r'chr[0-9XY]+:\d+\s', tmp).span()[0]:re.search(r'chr[0-9XY]+:\d+\s', tmp).span()[1]]
    #    tmp = html[re.search(r'<div\sclass=\"summary-box', html).span()[0]:re.search(r'<div\sclass=\"summary-box', html).span()[1]]
    #    print tmp
        #Smap.write(snp+"\t"+tmp[re.search(r'chr[0-9XY]+:\d+\s', tmp).span()[0]:re.search(r'chr[0-9XY]+:\d+\s', tmp).span()[1]]+"\n")
    #else:
    #    print "no snp"
