library(BiocInstaller)
source("https://bioconductor.org/biocLite.R")
biocLite('FRGEpistasis')
library(FRGEpistasis)
library(fda)
library(MASS)
library(stats)
library(factoextra)
#pergene <- fourierExpansion(gene,newgeno,newglst,newmap,rng)

snp2gene = read.csv("snp_info.txt", sep = '\t', header = T)
gene_list = unique(snp2gene$Gene)
gene_map = read.csv("gene_map.txt", sep = '\t', header = F)
gene_map = unique(gene_map[gene_map$V1%in%gene_list,])
colnames(gene_map) = c("Gene_Symbol", "Chromosome", "Start", "End")
snp_map = read.csv("newmap.csv", sep = "\t", header = F)
colnames(snp_map) = c("chr", "snp_id", "distant", "posi")
rng = 0
geno_num = read.csv("geno_num.2.txt", sep = '\t', header = T)
geno_num = geno_num[,-c(57,95)]
snp_list = colnames(geno_num)[-c(1,2)]
snp.id = snp_map$snp_id
snp.pos = as.vector(snp_map$posi)
names(snp.pos) = as.vector(snp.id)
for(gene in gene_list){
  gene = "DRD1"
  snp.in.gene = snp2gene[snp2gene$Gene%in%gene,]$Assays.name
  print(gene)
  geno.num.gene = cbind(geno_num[,c(1,2)], geno_num[,colnames(geno_num)%in%snp.in.gene])
  geno.num.gene[geno.num.gene<0]=NA
  geno.num.gene = na.omit(geno.num.gene)
  geno.num.gene.stat = geno.num.gene[,c(1,2)]
  geno.num.gene = geno.num.gene[,-c(1,2)]
  print(dim(geno.num.gene))
  nsnps = dim(geno.num.gene)[2]
  gil = list()
  pos = as.numeric(snp.pos[names(snp.pos)%in%colnames(geno.num.gene)])
  if(nsnps>1) {
    if ( is.null(pos) ) {
      pos <-(0:(nsnps-1))/(nsnps-1)
    }
    else {
      idx = order(pos)
      geno.num.gene = geno.num.gene[,idx]
      pos = pos[idx]
      pos<-(pos-pos[1])/(pos[nsnps]-pos[1])
      print(pos)
    }
    # use PCA to determine the number of basis
    eigenval<-prcomp(geno.num.gene)$sd^2
    print(gene)
    sum_eigen=sum(eigenval)
    tmp=0
    n_of_basis=0
    for(i in 1:length(eigenval))
    {
      tmp=eigenval[i]+tmp
      n_of_basis=i;
      if(tmp>=0.8*sum_eigen) {
        break
      }
    }
    n_of_basis=floor(n_of_basis/2)*2+1 
    print(n_of_basis)
    #make n_of_basis_A the odd number
    #end of setting the basis number
    frange <-c(pos[1], pos[length(pos)])
    fbasis<-create.fourier.basis(frange,nbasis=n_of_basis)
    phi=eval.basis(pos,fbasis);
    gil$zeta = t(ginv(t(phi)%*%phi)%*%t(phi)%*%t(geno.num.gene))
    gil$zeta = cbind(geno.num.gene.stat, gil$zeta)
  }
  else {
    gil$zeta  <-cbind(geno.num.gene.stat, geno.num.gene)
  }
  filename = paste(gene, "fex.csv", sep = "_")
  write.csv(gil$zeta, file = filename)
}

print(gene_list)
gene_list1 = gene_list[-11]
pval.mat = matrix(data = NA, nrow = length(gene_list1), ncol = length(gene_list1))
rownames(pval.mat) = colnames(pval.mat) = gene_list1
for(gene1 in 1:length(gene_list1)) {
  gene.fex1 = read.csv(paste("Fex/", gene_list1[gene1], "_fex.csv", sep = ""))
  n1 = dim(gene.fex1)[2]-3
  #print(n1)
  for(gene2 in 1:length(gene_list1)) {
    if (gene2 == gene1){
      pval.mat[gene1, gene2] = 1
      next
    }
    gene.fex2 = read.csv(paste("Fex/", gene_list1[gene2], "_fex.csv", sep = ""))
    n2 = dim(gene.fex2)[2]-3
    gene.fex = merge(gene.fex1, gene.fex2, by=c("Sample.ID", "Stat"), all=FALSE)
    pval.vec = rep(NA, times = n1*n2)
    n = 1
    for(i in 4:(4+n1-1)) {
      for(j in (4+n1+1):(4+n1+n2)) {
        #print(anova(glm(gene.fex[,2]~gene.fex[,i]*gene.fex[,j], family="binomial"), test="Chisq"))
        pval.vec[n] = anova(glm(gene.fex[,2]~gene.fex[,i]*gene.fex[,j], family="binomial"), test="Chisq")[4,5]
        n = n+1
      }
    }
    print(gene1)
    print(gene2)
    print(min(pval.vec))
    pval.mat[gene1, gene2] = pval.mat[gene2, gene1] = min(pval.vec)
  }
}

write.csv(pval.mat, file = "pvalmatrix.csv")
pvalmat1 = read.csv("result_p_aggregator.csv", header = T)
pvalmat1 = rbind(pvalmat1, rep(NA, times = 17))
colnames(pvalmat1)[6] = "G72/G30.txt"
rownames(pvalmat1) = colnames(pvalmat1)
colnames(pvalmat1)
rownames(pvalmat1)
for(i in 1:17){
  for(j in i:17){
    if (i==j){
      pvalmat1[j,i]=1
    }
    pvalmat1[j,i] = pvalmat1[i,j]
  }
}
#pval.mat[rownames(pval.mat)%in%gene_list1[1],colnames(pval.mat)%in%gene_list1[2]]
#pvalmat1[rownames(pvalmat1)%in%paste(gene_list1[1], "txt", sep="."), colnames(pvalmat1)%in%paste(gene_list1[1], "txt", sep=".")]
write.csv(pvalmat1, file = "pvalmatrixbefore.csv")
pcomp = matrix(data = NA, nrow = 2, ncol = length(gene_list1)*length(gene_list1))
comp = 1
p1 = 0
p2 = 0
for(gene1 in 1:length(gene_list1)) {
  for(gene2 in 1:length(gene_list1)){
    pval1 = pval.mat[rownames(pval.mat)%in%gene_list1[gene1],colnames(pval.mat)%in%gene_list1[gene2]]
    pval2 = pvalmat1[rownames(pvalmat1)%in%paste(gene_list1[gene1], "txt", sep="."), colnames(pvalmat1)%in%paste(gene_list1[gene2], "txt", sep=".")]
    #print(gene_list1[gene1])
    #print(gene_list1[gene2])
    #print(pval1)
    #print(pval2)
    if (pval1<0.05){
      print(paste("1", gene_list1[gene1], gene_list1[gene2], pval1, pval2, sep = "__"))
      p1 = p1 + 1
    }
    if (pval2<0.05){
      print(paste("2", gene_list1[gene1], gene_list1[gene2], pval1, pval2, sep = "__"))
      p2 = p2 + 2 
    }
    pcomp[1,comp] = pval1
    pcomp[2,comp] = pval2
    comp = comp+1
  }
  #break
}
plot(pcomp[1,],pcomp[2,])
x = c(0:1)
y = x
lines(y~x)
pcomp1=pcomp[,which(pcomp[1,]<0.05 | pcomp[2,]<0.05)]
plot(pcomp1[1,], pcomp1[2,])
lines(y~x)
plot(pcomp2[1,], pcomp2[2,])
lines(y~x)
plot(density(COMT.fex[which(COMT.fex$Stat==1),]$X1),xlim=c(-1,3), ylim=c(0,2), col=1)
lines(density(COMT.fex[which(COMT.fex$Stat==0),]$X1), col=2)
lines(density(COMT.fex[which(COMT.fex$Stat==1),]$X2), col=3)
lines(density(COMT.fex[which(COMT.fex$Stat==0),]$X2), col=4)
lines(density(COMT.fex[which(COMT.fex$Stat==1),]$X3), col=5)
lines(density(COMT.fex[which(COMT.fex$Stat==0),]$X3), col=6)
lines(density(COMT.fex[which(COMT.fex$Stat==1),]$X4), col=7)
lines(density(COMT.fex[which(COMT.fex$Stat==0),]$X4), col=8)
lines(density(COMT.fex[which(COMT.fex$Stat==1),]$X5), col=9)
lines(density(COMT.fex[which(COMT.fex$Stat==0),]$X5), col=10)
lines(density(COMT.fex[which(COMT.fex$Stat==1),]$X4))
lines(density(COMT.fex[which(COMT.fex$Stat==1),]$X5))
