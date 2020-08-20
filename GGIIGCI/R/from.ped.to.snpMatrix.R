#!/usr/bin/env Rscript

setwd("~/Desktop/data_from_HG/my_own_package/GGIGCI/")

read.ped.to.snpMatrix <- function(ped, info) {
  #selfhead(ped, nrow = 7)
  res = list()
  #require(snpStats)
  if (class(ped)!="character") {
    stop("Please, enter a valid path file for genotype data.")
  }
  file.type <- substr(ped, nchar(ped)-3, nchar(ped))
  if (file.type!=".ped") {
    stop(paste0("The file type of ", ped, " is not .ped.\n"))
  }
  if (dim(info)[1]!=4) {
    stop("The cols of info should be Chr, Gene, Snp and Position.\n")
  }
  if (dim(info)[2]!=dim(ped)[1]-6) {
    stop("The number of SNPs in info isn't equal to the number of SNPs in the ped data.\n")
  }
  
  
  return(res)
}
