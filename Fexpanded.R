#g1 = read.csv("clear_SLC6A4.txt", sep = '\t', header = F)
#g2 = read.csv("clear_TH.txt", sep = '\t', header = F)
#pheo = g1$V2
#g1d = g1[4:ncol(g1)-1]
#g2d = g2[4:ncol(g2)-1]
#Aggregator(Y=pheo, g1d, g2d)$p.value

geno_info=read.table("~/Desktop/学位论文/FRGEpistasis/inst/extdata/simGeno-chr2.raw", header = T)
gLst = read.csv("~/Desktop/学位论文/FRGEpistasis/inst/extdata/gene.list.csv")
library(snpStats)
snpinfo = read.table("~/Desktop/学位论文/newFourier/data/inter.ped", header = F)
newped = snpStats::read.pedfile(file="~/Desktop/学位论文/newFourier/data/inter.ped", snps = "~/Desktop/学位论文/newFourier/data/inter.info")
names(snpinfo)[1:6]<-names(geno_info)[1:6]
newposi = read.table("~/Desktop/学位论文/newFourier/data/inter.txt", header = T)
newmap<-newposi
newmap[,1]<-newposi[,3]
newmap[,3]<-rep(0,dim(newmap)[1])
names(newmap)<-c('chr','snp_id','distant','posi')
newmap$snp_id<-as.character(newmap$snp_id)
for(i in 1:dim(newmap)[1]){
  newmap[i,2]<-gsub('rs','r',newmap[i,2])
}
casesnp = read.csv("~/Desktop/学位论文/newFourier/data/case.csv")
casesnp$X.3<-rep(1,dim(casesnp)[1])
ctrlsnp = read.csv("~/Desktop/学位论文/newFourier/data/control.csv")
ctrlsnp$X.3<-rep(0,dim(ctrlsnp)[1])
names(ctrlsnp)[6:143]<-gsub('X','rs',names(ctrlsnp)[6:143])
names(casesnp)[6:143]<-gsub('r','rs',names(casesnp)[6:143])
snpmatrix<-rbind(casesnp,ctrlsnp)
##combine the table of case and control
newcol<-as.data.frame(matrix(rep(0,dim(snpmatrix)[1]),ncol = 1))
rownames(newcol)<-rownames(snpmatrix)
newcol<-as.data.frame(matrix(rep(0,dim(snpmatrix)[1]),ncol = 1))
rownames(newcol)<-rownames(snpmatrix)
#names(newcol)<-names(geno_info)[5]
fisnp<-cbind(snpmatrix[,1:4],newcol,snpmatrix[,-c(1:4)])##add a new col to the table of case+control
names(fisnp)[1:6]<-names(geno_info)[1:6]

#recode the snp with minor allele frequence
#load the prepared function lociencode
source('c:/Users/jiaxin/Desktop/??????ϵ????/newcasuality/recode_minor_Freq.R')


snpm<-lociencode(fisnp[,-c(1:6)])##convert the allele to the (0,1,2)num
snpm<-cbind(fisnp[,1:6],snpm)
snpm$IID<-snpm$FID
newpheno<-snpm[,c(2,6)]
genelist<-as.data.frame(table(newposi$Genenames))
genebed<-read.csv('c:/Users/jiaxin/Desktop/??????ϵ????/newfrge/mapfile.csv')

newglst<-cbind(genelist,genebed)[,-2]
names(newglst)<-names(gLst)
##define the extension scope of gene region
rng=0
fdr=0.05
## output data structure
out_epi <- data.frame( )
gene_info<-read.csv('c:/Users/jiaxin/Desktop/??????ϵ????/sorteddata/gene_info.csv')
formercase<-read.csv('c:/Users/jiaxin/Desktop/??????ϵ????/sorteddata/finalsnpm.csv')
#formercase<-cbind(pheno$PHENOTYPE,formercase);names(formercase)[1]<-'pheno'
missing<-names(formercase)%in%gene_info$SNPnames
miss_snp<-names(formercase)[!missing]
totalsnp<-c(as.vector(gene_info$SNPnames),miss_snp)
#totalsnp<-c(as.vector(gene_info$SNPnames),miss_snp)
searchedname<-'MAOA MAOA SLC1A2 NR3C1 CYP2D6 LOC100288866 GRIN2A AKT3 BDNF FAM MAPK8 GRM7 GAD1 GRIN2B ABCB6 GRIK4 ABCB1 NTRK2 ABCB1 NULL NULL MC2R ABCB1 CRHBP CYP2C9'
searchedname<-strsplit(searchedname,split = ' ')
totalgene<-c(as.character(gene_info$Genenames),searchedname[[1]])
#gene_info$SNPnames<-totalsnp
rawvec<-totalgene
names(rawvec)<-totalsnp
newgeno<-snpm[,-c(1:6)]
glist<-as.character(rawvec[names(newgeno)])
#newposi<-subset(newposi,SNPnames%in%names(newgeno))

#define a posi vecector instead use the newposi
posivec<-as.vector(newposi$Position)
names(posivec)<-as.vector(newposi$SNPnames)
geno_expn=matrix(0,dim(snpm)[1],1)
geno_expn=geno_expn[,-1]
newglst<-read.csv('c:/Users/jiaxin/Desktop/??????ϵ????/sorteddata/newglst.csv')
#newgeno<-read.csv('c:/Users/jiaxin/Desktop/??????ϵ????/newfrge/snpm.csv')
newgeno<-snpm[,-c(1:6)]
for(gene in unique(newglst$Gene_Symbol)){
  #for(gene in unique(newglst$Gene_Symbol)){
  #gene<-as.numeric(gene)
  #print(gene)
  gene<-as.character(gene)
  pergene <- fourierExpansion(gene,newgeno,newglst,newmap,rng)
  colnames(pergene)<-paste('gene',gene,seq(1:dim(pergene)[2]),sep = '_')
  geno_expn<-cbind(geno_expn,pergene)
  print(paste('gene',gene,'done!',sep = ' '))
  
}

#write.csv(geno_expn,'c:/Users/jiaxin/Desktop/??????ϵ????/sorteddata/expanded_genotype.csv',
#          row.names = F)

sam_info = read.csv("sam_info.txt", sep = '\t', header = T)
snp_info = read.csv("snp_info.txt", sep = '\t', header = T)
my_gene = unique(snp_info$Gene)




