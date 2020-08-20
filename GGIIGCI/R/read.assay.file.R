#!/usr/bin/env Rscript

#setwd("~/Desktop/macdesktop/data_from_HG/my_own_package/GGIGCI/")

read.assay.file <- function(dir="./", pattern, plate.num = 4, skip = 11, out.plate = F){
  
  #print(dir)
  plate.list <- paste0("plate", c(1:plate.num))
  files.plate <- dir(dir, pattern = pattern)
  output <- NA
  for(plate in plate.list) {
    #file.plate <- dir(dir, pattern = pattern)
    temp.output <- NA
    for(file in files.plate) {
      temp.strs <- strsplit(file, split = "-")[[1]]
      temp.snp.name <- temp.strs[3]
      temp.plate <- temp.strs[4]
      temp.gene.name <- temp.strs[2]
      
      if (temp.plate==plate){
        message(paste0("merging the site of ", temp.gene.name, ":", temp.snp.name))
        temp.file.data <- read.table(file = paste0(dir, file), skip = skip, sep = "\t", header = T)
        #message(dim(temp.file.data))
        temp.data <- temp.file.data[, c(1,2,6)]
        plate.vector <- rep(temp.plate, dim(temp.data)[1])
        #temp.data <- cbind(plate.vector, temp.data)
        colnames(temp.data)[3] <- temp.snp.name
        if (!is.data.frame(temp.output)){
          temp.output <- temp.data
          if (out.plate) {
            temp.output <- cbind(plate.vector, temp.output)
          }
        }
        else {
          temp.output <- merge(temp.output, temp.data, by = c("Well", "Sample.Name"), all.x = T)
        }
      }
    }
    if (!is.data.frame(output)){
      output <- temp.output
    }
    else {
      output <- rbind(output, temp.output)
    }
  }
  return(output)
}
