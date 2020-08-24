#!/usr/bin/env Rscript

setwd("~/Desktop/data_from_HG/my_own_package/GGIGCI/")

gene.gene.interaction <- function(data, n.times = 1000, selected.gene.list=NULL){
  
  if (!inherits(data, "GGIGCI.data")) {
    stop("data must be one of the class 'GGIGCI.data'(the output of the function impute.res).")
  }
  snp.info <- data$SnpMatrix
  gene.info <- data$gene.info
  sample.info <- data$sample.info
  #phen <- sample.info[,6]
  if (is.null(selected.gene.list)) {
    selected.gene.list <- unique(gene.info[, 2])
  }
  else if (any(!selected.gene.list%in%gene.info[,2])) {
    stop("invaild gene selected.")
  }
  #else {
  gene.info <- gene.info[which(gene.info[, 2]%in%selected.gene.list), ]
  #}
  if (any(sample.info[,2]!=rownames(snp.info))){
    warning("miss match between sample info matrix and SnpMatrix.")
  }
  phen <- sample.info[,6]
  p.val.matrix <- diag(0, nrow = length(selected.gene.list))
  colnames(p.val.matrix) <- selected.gene.list
  rownames(p.val.matrix) <- selected.gene.list
  s.val.matrix <- diag(0, nrow = length(selected.gene.list))
  colnames(s.val.matrix) <- selected.gene.list
  rownames(s.val.matrix) <- selected.gene.list
  
  interaction.pairs <- combn(selected.gene.list, 2)
  gene.tabe <- data.frame(Genenames = selected.gene.list, Snps = table(gene.info[,2])[selected.gene.list])[,-2]
  ggi.ud <- data.frame(from = interaction.pairs[1, ], to = interaction.pairs[2, ], 
                       P.val = rep(0, length(interaction.pairs[1, ])), S.val = rep(0, length(interaction.pairs[1, ])))
  nc.i <- ncol(interaction.pairs)
  
  for(n in 1:nc.i){
    gene1 <- interaction.pairs[1, n]
    gene2 <- interaction.pairs[2, n]
    snp.gene1 <- gene.info[which(gene.info[, 2]==gene1), 3]
    snp.gene2 <- gene.info[which(gene.info[, 2]==gene2), 3]
    snp.info.gene1 <- snp.info[, colnames(snp.info)%in%snp.gene1]
    snp.info.gene2 <- snp.info[, colnames(snp.info)%in%snp.gene2]
    temp.result <- GBIGM(phen = phen, gene1 = snp.info.gene1, gene2 = snp.info.gene2, n.times = n.times)
    p.val.matrix[gene1, gene2] <- p.val.matrix[gene2, gene1] <- temp.result[1]
    s.val.matrix[gene1, gene2] <- s.val.matrix[gene2, gene1] <- temp.result[2]
    ggi.ud[n, 3:4] <- temp.result
  }
  
  snp.matrix <- snp.info[, colnames(snp.info)%in%gene.info[,3]]
  out <- list(SnpMatrix = snp.matrix, gene.info = gene.info, sample.info = sample.info, 
              GGI.ud = ggi.ud, GENE.table = gene.tabe, P.matrix = p.val.matrix, S.matrix = s.val.matrix)
  class(out) <- c("GGIUD", "list")
  return(out)
}