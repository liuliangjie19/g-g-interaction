#!/usr/bin/env Rscipt

setwd("~/Desktop/data_from_HG/my_own_package/GGIGCI/")

convert.sds.to.allele <- function(rs, allele1 = NA, allele2 = NA, is.num = F) {
  if (is.num) {
    vic <- "1 1"
    fam <- "2 2"
    both <- "1 2"
  } else {
    if (is.na(allele1) | is.na(allele2)){
      stop("allele1 and allele2 are required when is.num is FALSE.\n")
    }
    vic <- paste0(allele1, " ", allele1)
    fam <- paste0(allele2, " ", allele2)
    both <- paste0(allele1, " ", allele2)
  }
  miss <- "0 0"
  out <- rep(NA, length(rs))
  for(i in 1:length(rs)) {
    n <- rs[i]
    if (n=="vic") {
      out[i] <- vic
    } 
    else if (n == "fam") {
      out[i] <- fam
    }
    else if (n=="Both") {
      out[i] <- both
    }
    else {
      out[i] <- miss
    }
  }
  return(out)
}

from.assay.to.ped <- function(assay, snp.info, sample.info, is.sample.info = F, 
                              col.names.convert = F, fam.id = F, is.num = T, write.file = F){
  
  if (!is.data.frame(assay) | !is.data.frame(snp.info)){
    stop("assay or snp.info must be a data.frame!\n")
  }
  if (col.names.convert) {
    if (!"ABI.ID"%in%colnames(snp.info)){
      stop("The information of ABI.ID is required when col.names.convert is TRUE.\n")
    }
  }
  if (is.sample.info) {
    if (!is.data.frame(sample.info)){
      stop("The sample.info is required when is.sample.info is TRUE.\n")
    }
  }
  if (fam.id) {
    if (is.data.frame(sample.info)){
      stop("The sample.info is required when fam.id is TRUE.\n")
    }
  }
  
  if(col.names.convert) {
    if (dim(snp.info[which(is.na(snp.info$ABI.ID)& is.na(snp.info$Assays.name)), ])[1]!=0){
      stop("Required information of ABI.ID or Assays.name missing.\n")
    }
    for(i in 1:dim(snp.info)[1]) {
      if (is.na(snp.info[i, ]$Assays.name)) {
        snp.info[i, ]$Assays.name <- snp.info[i, ]$ABI.ID
      }
      if (is.na(snp.info[i, ]$ABI.ID)) {
        snp.info[i, ]$ABI.ID <- snp.info[i, ]$Assays.name 
      }
    }
    old.cnames <- snp.info$ABI.ID
    new.cnames <- snp.info$Assays.name
    temp.cname <- new.cnames
    
    names(temp.cname) <- old.cnames
    #print(temp.cname)
    for(i in 1:dim(assay)[2]) {
      n <- colnames(assay)[i]
      if (n%in%names(temp.cname)){
        #print(temp.cname[n])
        colnames(assay)[i] <- temp.cname[n]
      }
    }
    #print(colnames(assay))
  }
  
  #out <- list()
  #ped <- NA
  c <- dim(assay)[1]
  r <- dim(assay)[2]
  ped <- data.frame(Fam.id = assay$Sample.Name, Sample.Name = assay$Sample.Name, 
                    Paternal.ID = rep(0, c), Maternal.ID = rep(0, c), Sex = rep(0, c),
                    Phen = rep(-9, c))
  for (i in 3:r) {
    temp.snp <- colnames(assay)[i]
    allele1 <- snp.info[which(snp.info$Assays.name==temp.snp), ]$allele1
    allele2 <- snp.info[which(snp.info$Assays.name==temp.snp), ]$allele2
    temp.snp.vector <- convert.sds.to.allele(assay[, i], allele1 = allele1, allele2 = allele2, is.num = is.num)
    ped[, i+4] <- temp.snp.vector
    colnames(ped)[i+4] <- temp.snp
  }
  
  if (is.sample.info) {
    ped <- ped[, -c(1,3,4,5,6)]
    ped <- merge(sample.info, ped, by = c("Sample.Name"), all.x = T)
  }
  
  return(ped)
}

