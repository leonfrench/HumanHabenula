library(AnnotationDbi)
library(annotate)
library(org.Hs.eg.db)
library(readxl)
detach("package:dplyr", unload=TRUE)
library(dplyr)
library(magrittr)
library(limma)
library(ggplot2)
library(readr)

#deprecated

sourceExpression <- "allen_HBA"
#probe_strategy <- "reannotator"
probe_strategy <- "symbolsQC"
allsampleAnnot = NULL
allExpression = NULL



for (donorFolder in list.files(paste0("./data/raw/",sourceExpression,"/"), pattern = "normalized_microarray.*")) {
  sampleAnnot <- read_csv(paste0("./data/raw/",sourceExpression,"/",donorFolder,"/SampleAnnot.csv"))
  
  #check if donor has habenula
  if (nrow(filter(sampleAnnot, grepl("lateral habenula", structure_name))) > 0) {
    expressionMatrix <- read_csv(paste0("./data/raw/",sourceExpression,"/",donorFolder, "/MicroarrayExpression.csv"), col_names=F) 
    probeInfo <- read_csv(paste0("./data/raw/",sourceExpression,"/",donorFolder,"/Probes.csv"))
    
    expressionMatrix %<>% rename(probe_id = X1)
    dim(expressionMatrix)
    
    strip_left_right <- function(structure_name) {
      tokens <- trimws(unlist(strsplit(structure_name, ",")))
      tokens <- tokens[tokens != "left"]
      tokens <- tokens[tokens != "right"]
      cleaned_name <- paste(tokens, collapse = ", ")
      cleaned_name
    }
    sampleAnnot %<>% rowwise() %>% mutate(structure_name_left_right_stripped = strip_left_right(structure_name))
    sampleAnnot %<>% mutate(donorID = donorFolder)
    sampleAnnot %<>% mutate(uniqueID = paste("ID", structure_id, slab_num, well_id, polygon_id, donorID, sep=".")) %>% select(uniqueID, everything())
    #sampleAnnot$order <- paste0("X", 2:(nrow(sampleAnnot)+1))
    colnames(expressionMatrix) <- c("probe_id", sampleAnnot$uniqueID)
    
    expressionMatrix <- inner_join(probeInfo %>% select(probe_id, probe_name), expressionMatrix) %>% select(-probe_id)
    
    #bind cols of expression matrix
    allExpression <- bind_cols(allExpression, expressionMatrix)
    
    #bind rows of sample annot
    allsampleAnnot <- bind_rows(allsampleAnnot, sampleAnnot)
  }
}
#set in same order
sampleAnnot <- allsampleAnnot
expressionMatrix <- allExpression[, c("probe_name", sampleAnnot$uniqueID)]

length(unique(sampleAnnot$structure_name_left_right_stripped))

#setup design file
sampleAnnot %<>% mutate(isHabenula = grepl("lateral habenula", structure_name))
sampleAnnot %>% filter(isHabenula) %>% select(donorID, everything())

sampleAnnot %>% select(uniqueID, isHabenula)
designMatrix <- model.matrix(~ sampleAnnot$isHabenula + sampleAnnot$donorID)

expressionMatrix <- as.data.frame(expressionMatrix)
rownames(expressionMatrix) <- expressionMatrix$probe_name
expressionMatrix$probe_name <- NULL

fit <- lmFit(expressionMatrix,designMatrix)
fit <- eBayes(fit)

topTable(fit,coef="sampleAnnot$isHabenulaTRUE")

limmaResults <- inner_join(as_tibble(fit$p.value[,"sampleAnnot$isHabenulaTRUE", drop=F], rownames="probe_name"), as_tibble(fit$t[,"sampleAnnot$isHabenulaTRUE", drop=F], rownames="probe_name"), by="probe_name")
limmaResults %<>% rename( p.value=`sampleAnnot$isHabenulaTRUE.x`, t = `sampleAnnot$isHabenulaTRUE.y`)

if(probe_strategy =="reannotator") {
  reannotator <- read_tsv("./data/raw/gene_symbol_annotations/AllenInstitute_custom_Agilent_Array.txt")
  reannotator %<>% select(probe_name=`#PROBE_ID`, Gene_symbol)
  #uses semicolon to separate genes
  reannotator %<>% mutate(Gene_symbol = gsub(";", "|", Gene_symbol))
  
  length(intersect(reannotator$probe_name, probeInfo$probe_name))
  length(setdiff(reannotator$probe_name, probeInfo$probe_name))
  length(setdiff(probeInfo$probe_name, reannotator$probe_name))
  
  #join
  limmaResults <- left_join(limmaResults, reannotator) %>% select(probe_name, Gene_symbol, p.value, t)
  limmaResults <- left_join(limmaResults, probeInfo %>% dplyr::select(probe_name, allenEntrez = entrez_id))
} else if (probe_strategy =="passQC") {
  #don't use - this losses probes that are habenula specific
  qc_table <- read_xlsx("./data/raw/gene_symbol_annotations/Miller et al. doi.org_10.1186_1471-2164-15-154 12864_2013_7016_MOESM8_ESM.xlsx",skip=1)
  qc_table %<>% select(probe_name = X__1, gene_symbol = X__2, m, b, is_qc_pass = `Pass?`)
  qc_table %>% group_by(is_qc_pass) %>% summarize(n=n())
  #17389 entrez IDs in the mri transcriptomics
  
  print(paste("Probes before filtering", nrow(qc_table)))
  qc_table %<>% filter(is_qc_pass != "FALSE") #count NA values as passing
  print(paste("Probes after filtering for QC", nrow(qc_table)))
  length(unique(qc_table$gene_symbol))
  
  
  #use gene symbol to get NCBI ID, then get updated symbol
  symbolToID <- probeInfo %>% select(gene_symbol, entrez_id) %>% distinct()
  symbolToID %<>% filter(!is.na(entrez_id)) 
  symbolToID %<>% mutate(new_symbol = getSYMBOL(as.character(entrez_id), data='org.Hs.eg')) 
  
  qc_table <- left_join(qc_table, symbolToID)
  qc_table %<>% mutate(legacySymbol = gene_symbol) #save old symbol
  qc_table %<>% mutate(gene_symbol = if_else(is.na(new_symbol), gene_symbol, new_symbol)) %>% dplyr::select(gene_symbol, entrez_id,  everything(), -new_symbol)
  
  qc_table <- qc_table %>% filter(!grepl("A_", gene_symbol)) %>% filter(!grepl("CUST_", gene_symbol)) 
  
  #31452
  print(paste("Probes after filtering for geneSymbol", nrow(qc_table)))
  print(length(unique(qc_table$gene_symbol)))
  
  limmaResults <- inner_join(limmaResults, qc_table)
} else if(probe_strategy =="symbolsQC") {
  qc_table <- read_xlsx("./data/raw/gene_symbol_annotations/Miller et al. doi.org_10.1186_1471-2164-15-154 12864_2013_7016_MOESM8_ESM.xlsx",skip=1)
  qc_table %<>% select(probe_name = X__1, gene_symbol = X__2, is_qc_pass = `Pass?`)
  
  #use gene symbol to get NCBI ID, then get updated symbol
  symbolToID <- probeInfo %>% select(gene_symbol, entrez_id) %>% distinct()
  symbolToID %<>% filter(!is.na(entrez_id)) 
  symbolToID %<>% mutate(new_symbol = getSYMBOL(as.character(entrez_id), data='org.Hs.eg')) 
  
  qc_table <- left_join(qc_table, symbolToID)
  qc_table %<>% mutate(legacySymbol = gene_symbol) #save old symbol
  qc_table %<>% mutate(gene_symbol = if_else(is.na(new_symbol), gene_symbol, new_symbol)) %>% dplyr::select(gene_symbol, entrez_id,  everything(), -new_symbol)
  
  qc_table <- qc_table %>% filter(!grepl("A_", gene_symbol)) %>% filter(!grepl("CUST_", gene_symbol)) 
  #remove numeric gene symbols - probably excel date conversions
  qc_table %<>% filter(is.na(as.numeric(gene_symbol))) 
  
  print(paste("Probes after filtering for _ and numeric geneSymbols", nrow(qc_table)))
  print(paste("Gene count",length(unique(qc_table$gene_symbol))))
  
  limmaResults <- inner_join(limmaResults, qc_table)
}

limmaResults %<>% mutate(p.adj = p.adjust(p.value, method="fdr"))
length(unique(limmaResults %>% filter(p.adj < 0.05) %>% .$gene_symbol ))
limmaResults %>% filter(p.adj < 0.05, t > 0) %>% arrange(p.adj) %>% head()
limmaResults %>% filter(p.adj < 0.05, t > 0) %>% arrange(p.adj) %>% head()
limmaResults %>% filter(grepl("GPR151", gene_symbol))

#summarize to genes, use lowest p-value
gene_summary <- limmaResults %>% group_by(gene_symbol) %>% arrange(p.value) %>% 
  summarize(p.value = first(p.value), direction=sign(first(t)))
#convert to ranks
gene_summary %<>% mutate(pValueWithDirection = direction * (nrow(gene_summary) - rank(p.value)))
#create a region by gene matrix
gene_summary %>% select(gene_symbol, targetRegion = pValueWithDirection)

write_csv(limmaResults, paste0("./results/limma/Habenula.3donors.lateral.",sourceExpression,".csv"))
