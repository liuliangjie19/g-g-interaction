setwd("~/Desktop/data_from_HG/my_own_package/GGIGCI/")

source("R/from.assay.to.ped.R")
source("R/read.assay.file.R")
source("R/igci.R")
source("R/read.ped.to.snpMatrix.R")
source("R/impute.res.R")
source("R/GBIGM.without.phen.R")
source("R/GBIGM.R")
source("R/gene.gene.interaction.R")
source("R/gene.gene.causal.R")
require(openssl)

data.total.1536 <- read.assay.file(dir = "data_first_20190116/", pattern = "AZ")
selfhead(data.total.1536)
snp.info <- read.xls("~/Desktop/data_from_HG/第一批练习.xls", sheet = 3)
sample.info <- read.xls("~/Desktop/data_from_HG/第一批练习.xls", sheet = 2)

snp.info <- snp.info[which(snp.info$Status=="OK"), ]
sample.info <- sample.info[which(sample.info$Sex!=""), ]
for(i in 1:dim(snp.info)[1]) {
  x <- snp.info$ABI.ID[i]
  snp.info$ABI.ID[i] <- gsub("_+", "_", x)
  #print(x)
}

data.ped <- from.assay.to.ped(data.total.1536, snp.info)
sample.info <- data.frame(Fam.id = sample.info$Unique.Identifier, Sample.Name = sample.info$Unique.Identifier, Paternal.ID = rep(0, dim(sample.info)[1]),
                          Maternal.ID = rep(0, dim(sample.info)[1]), Sex = sample.info$Sex, Phen = sample.info$Classify.)
sample.info[which(sample.info$Sex=="M"), ]$Sex <- 1
sample.info[which(sample.info$Sex=="F"), ]$Sex <- 2
sample.info[which(sample.info$Phen=="NZ"), ]$Phen <- 1
sample.info[which(sample.info$Phen=="SZ"), ]$Phen <- 2

data.ped <- from.assay.to.ped(data.total.1536, snp.info = snp.info, sample.info = sample.info, 
                              is.sample.info = T)
write.table(data.ped, file = "data/data.ped", quote = F, sep = "\t", col.names = F, row.names = F)
write.table(data.ped, file = "data/data.with.header.ped", quote = F, sep = "\t", col.names = T, row.names = F)
#res <- read.ped.to.snpMatrix(ped = "data/data.ped", info = snp2gene, 
                             snp.list = snp.list, ped.header = F, ped.sep = "\t")
res <- read.ped.to.snpMatrix(ped = "data/data.with.header.ped", ped.header = T, info = snp2gene)

temp.im <- impute.res(res, removed = "sample")
#removed 45 sample in the raw data

#gene2:COMT
#gene1:DLG4

temp.snp.list1 <- temp.im$gene.info[which(temp.im$gene.info$Genenames==gene1), ]$SNPnames
temp.snp.list2 <- temp.im$gene.info[which(temp.im$gene.info$Genenames==gene2), ]$SNPnames
G1 <- temp.im$SnpMatrix[, colnames(temp.im$SnpMatrix)%in%temp.snp.list1]
G2 <- temp.im$SnpMatrix[, colnames(temp.im$SnpMatrix)%in%temp.snp.list2]
phen <- temp.im$sample.info$affected
#GBIGM.without.phen(G1, G2, n.times = 100)
GBIGM(phen, G1, G2, n.times = 100)
temp.ggiud <- gene.gene.interaction(temp.im, n.times = 1000, selected.gene.list = c("DRD1", "DRD2", "DRD3"))
gene.gene.causal(temp.ggiud, n.times = 100, sample.size = 700, fex = F)


