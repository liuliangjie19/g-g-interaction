setwd("~/Desktop/data_from_HG/my_own_package/GGIGCI/")

source("R/from.assay.to.ped.R")
source("R/read.assay.file.R")
source("R/igci.R")
source("R/from.ped.to.snpMatrix.R")
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
write.table(data.ped, file = "data/data.ped", quote = F, sep = "\t", col.names = F)
read.ped.to.snpMatrix("data_first_20190116/data.ped", "data_first_20190116/data.txt")
