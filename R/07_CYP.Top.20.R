library(readr)
library(tmod)
library(magrittr)
detach("package:dplyr", unload=TRUE)
library(dplyr)
library(foreach)
library(doMC)
registerDoMC(cores=6)

#CYP genes
matrix <- read_tsv("./data/processed/allen_HBA_brainarea_vs_genes_exp_qcNames.tsv", guess_max=22000)
length(unique(matrix$gene_symbol))
length(unique(adult$gene_symbol))

matrix %<>% mutate(isCYP = grepl("^CYP[0-9]" , gene_symbol))
sort((matrix %>% filter(isCYP))$gene_symbol)
length((matrix %>% filter(isCYP))$gene_symbol)
cypGenes <- sort((matrix %>% filter(isCYP))$gene_symbol)

regionNames <- setdiff(colnames(matrix), c("gene_symbol", "isCYP"))

allCYPResults<-foreach(region=regionNames,.combine=rbind) %dopar% {
  matrix %<>% arrange_(paste0("desc(`",region,"`)"))
  sortedGenes <- matrix$gene_symbol
  top20 <- head(sortedGenes, n=20)
  intersectSize <- length(intersect(top20, cypGenes))
  data.frame(region, intersectSize)
}

filter( allCYPResults, intersectSize > 2)
filter( allCYPResults, intersectSize > 1)
