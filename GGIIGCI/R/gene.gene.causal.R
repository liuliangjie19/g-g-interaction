#!/usr/bin/env Rscript

setwd("~/Desktop/data_from_HG/my_own_package/GGIGCI/")

gene.gene.causal <- function(data.ggiud, n.times = 1000, sample.size = 100, 
                             ref = c("Gaussian", "uniform"), est = c("entropy", "intergral"),
                             interaction.p = 0.05, fex = T) {
  #if (!inherits(data.ggigci, "GGIGCI.data")){
  #  stop("data must be one of the class 'GGIGCI.data'(the output of the function impute.res).")
  #}
  if (!inherits(data.ggiud, "GGIUD")){
    stop("data must be one of the class 'GGIUD'(the output of the function gene.gene.interaction).")
  }
  
  ref <- match.arg(ref)
  est <- match.arg(est)
  gene.table <- data.ggiud$GENE.table
  snp.info <- data.ggiud$SnpMatrix
  gene.info <- data.ggiud$gene.info
  sample.info <- data.ggiud$sample.info
  phen <- as(sample.info[,6], "numeric")
  selected.gene.pairs <- data.ggiud$GGI.ud[which(data.ggiud$GGI.ud$P.val<=interaction.p), ]
  gene.f.table <- data.frame(gene1 = selected.gene.pairs[, 1], gene2 = selected.gene.pairs[, 2], 
                             f1 = rep(NA, length(selected.gene.pairs[,1])), f2 = rep(NA, length(selected.gene.pairs[, 2])),
                             P.value = rep(0, length(selected.gene.pairs[, 1])), causal.direction = rep(0, length(selected.gene.pairs[, 1])))
  gene.causal.table <- data.frame(from = selected.gene.pairs[, 1], to = selected.gene.pairs[, 2], 
                                  F.value = rep(0, length(selected.gene.pairs[, 1])))
  n.pairs <- dim(selected.gene.pairs)[1]
  for(i in seq_len(n.pairs)){
    #print(selected.gene.pairs[i,])
    gene1 <- selected.gene.pairs[i, 1]
    gene2 <- selected.gene.pairs[i, 2]
    snp.gene1 <- gene.info[which(gene.info[,2]==gene1), 3]
    snp.gene2 <- gene.info[which(gene.info[,2]==gene2), 3]
    snp.info.gene1 <- as(snp.info[, colnames(snp.info)%in%snp.gene1], "numeric")
    snp.info.gene2 <- as(snp.info[, colnames(snp.info)%in%snp.gene2], "numeric")
    
    if (fex){
      x.gene1 <- fex(cbind(phen, snp.info.gene1))
      x.gene2 <- fex(cbind(phen, snp.info.gene2))
    }
    else {
      x.gene1 <- apply(cbind(phen, 2*base3to10(snp.info.gene1)), 1, sum)
      x.gene2 <- apply(cbind(phen, 2*base3to10(snp.info.gene2)), 1, sum)
    }
    pairs.data <- data.frame(gene1 = x.gene1, gene2 = x.gene2)
    temp.result <- IGCIfor2value(pairs.data, times = n.times, sample.size = sample.size, refMeasure = ref, estimator = est)
    
    gene.f.table[i, 3:6] <- c(temp.result$f.mean[1], temp.result$f.mean[2], temp.result$p.value, temp.result$causal.direction)
    if (gene.f.table[i, 6]=="gene1->gene2"){
      gene.causal.table[i, ] <- c(gene.f.table[i, 1], gene.f.table[i, 2], gene.f.table[i, 3])
    }
    else if (gene.f.table[i, 6]=="gene2->gene1") {
      gene.causal.table[i, ] <- c(gene.f.table[i, 2], gene.f.table[i, 1], gene.f.table[i, 4])
    }
  }
  out <- list(SnpMatrix = snp.info, gene.info = gene.info, sample.info = sample.info, gg.causal = gene.causal.table, gene.table = gene.table)
  class(out) <- c("GGC.data", "list")
  return(out)
}

fex <- function(gene){
  return(apply(gene, 1, sum))
}

base3to10 <- function(gene){
  n.SNP<-ncol(gene)
  long<-nrow(gene)
  tab<-matrix(NA,ncol=n.SNP,nrow=long)
  for (i in seq_len(n.SNP)){
    tab[,i]<-3^(i-1)*gene[,i]
  }
  b10<-apply(tab,1,sum)
  return(b10)
}
