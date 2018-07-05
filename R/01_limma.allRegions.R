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

#uncomment line below to set adult or fetal expression source
sourceExpression <- "allen_HBA"
#sourceExpression <- "allen_human_fetal_brain"

if (sourceExpression == "allen_HBA") {
  folderPattern <- "normalized_microarray.*"
  sampleFilename <- "SampleAnnot.csv"
  probeFilename <- "Probes.csv"
  expressionFilename <- "MicroarrayExpression.csv"
} else if (sourceExpression == "allen_human_fetal_brain") {
  folderPattern <- "lmd_matrix_.*"
  sampleFilename <- "columns_metadata.csv"
  probeFilename <- "rows_metadata.csv"
  expressionFilename <- "expression_matrix.csv"
}
allsampleAnnot = NULL
allExpression = NULL

for (donorFolder in list.files(paste0("./data/raw/",sourceExpression,"/"), pattern = folderPattern)) {
  sampleAnnot <- read_csv(paste0("./data/raw/",sourceExpression,"/",donorFolder,"/",sampleFilename))
  
  #check if donor has habenula
  if (nrow(filter(sampleAnnot, grepl("lateral habenula", structure_name))) > 0) {
    print(donorFolder)
  }
  if (nrow(filter(sampleAnnot, grepl("medial habenula", structure_name))) > 0) {
    print(paste("medial", donorFolder))
  }
}

#fetal brain has more/better gene symbol mappings
probeInfo <- read_csv(paste0("./data/raw/allen_human_fetal_brain/lmd_matrix_12566/rows_metadata.csv"))
probeInfo %<>% rename(probe_id = probeset_id, probe_name = probeset_name)

for (donorFolder in list.files(paste0("./data/raw/",sourceExpression,"/"), pattern = folderPattern)) {
  sampleAnnot <- read_csv(paste0("./data/raw/",sourceExpression,"/",donorFolder,"/",sampleFilename))
  
  #check if donor has habenula
  if (nrow(filter(sampleAnnot, grepl("lateral habenula", structure_name))) > 0) {
    expressionMatrix <- read_csv(paste0("./data/raw/",sourceExpression,"/",donorFolder, "/",expressionFilename), col_names=F) 
    
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
    if (sourceExpression == "allen_HBA") {
      sampleAnnot %<>% mutate(uniqueID = paste("ID", structure_id, slab_num, well_id, polygon_id, donorID, sep=".")) %>% select(uniqueID, everything())
    } else if (sourceExpression == "allen_human_fetal_brain") {
      sampleAnnot %<>% mutate(uniqueID = paste("ID", structure_id, well_id, donorID, sep=".")) %>% select(uniqueID, everything())
    }

    colnames(expressionMatrix) <- c("probe_id", sampleAnnot$uniqueID)
    
    expressionMatrix <- inner_join(probeInfo %>% select(probe_id, probe_name), expressionMatrix) %>% select(-probe_id)
    
    #write out for estimating cell-type proportions
    #needs gene symbol
    
    #bind cols of expression matrix
    allExpression <- bind_cols(allExpression, expressionMatrix)
    
    #bind rows of sample annot
    allsampleAnnot <- bind_rows(allsampleAnnot, sampleAnnot)
  }
}

#set in same order
sampleAnnot <- allsampleAnnot
expressionMatrix <- allExpression[, c("probe_name", sampleAnnot$uniqueID)]

print(paste("Number of unique regions:", length(unique(sampleAnnot$structure_name_left_right_stripped))))
print(paste("Number of samples:", nrow(sampleAnnot)))

expressionMatrix <- as.data.frame(expressionMatrix)
rownames(expressionMatrix) <- expressionMatrix$probe_name
expressionMatrix$probe_name <- NULL

#create probe mapping file here
qc_table <- read_xlsx("./data/raw/gene_symbol_annotations/Miller et al. doi.org_10.1186_1471-2164-15-154 12864_2013_7016_MOESM8_ESM.xlsx",skip=1)
qc_table %<>% select(probe_name = X__1, gene_symbol = X__2, is_qc_pass = `Pass?`)

#fix numeric gene symbols from excel date conversions
qc_table <- inner_join( qc_table, probeInfo %>% select(probe_name, for_excel_problem = gene_symbol))
qc_table %<>% mutate(gene_symbol = if_else(!is.na(as.numeric(gene_symbol)), for_excel_problem, gene_symbol)) %>% select(-for_excel_problem)

#use gene symbol to get NCBI ID, then get updated symbol
symbolToID <- probeInfo %>% select(gene_symbol, entrez_id) %>% distinct()
symbolToID %<>% filter(!is.na(entrez_id)) 
symbolToID %<>% mutate(new_symbol = getSYMBOL(as.character(entrez_id), data='org.Hs.eg')) 

qc_table <- left_join(qc_table, symbolToID)
qc_table %<>% mutate(legacySymbol = gene_symbol) #save old symbol
qc_table %<>% mutate(gene_symbol = if_else(is.na(new_symbol), gene_symbol, new_symbol)) %>% dplyr::select(gene_symbol, entrez_id,  everything(), -new_symbol)

#remove those not mapping
qc_table <- qc_table %>% filter(!grepl("A_", gene_symbol)) %>% filter(!grepl("CUST_", gene_symbol)) 
fetalProbes %<>% filter(gene_symbol != "na")

print(paste("Probes after filtering for _ and numeric geneSymbols", nrow(qc_table)))
print(paste("Gene count",length(unique(qc_table$gene_symbol))))

#handle rownames
write_tsv(expressionMatrix %>% mutate(probe_name = rownames(expressionMatrix)) %>% select(probe_name, everything()), 
          paste0("./data/processed/",sourceExpression, "_sampleMatrix_qcNames.tsv"))
write_tsv(sampleAnnot, paste0("./data/processed/",sourceExpression, "_sampleMatrix_colAnnotations_qcNames.tsv"))


regionByGene <- NULL
for (targetRegion in unique(sampleAnnot$structure_name_left_right_stripped)) {
  print(targetRegion)
  sampleAnnot %<>% mutate(isTarget = (targetRegion == structure_name_left_right_stripped))
  sampleAnnot %>% filter(isTarget) %>% select(donorID, everything())
  
  designMatrix <- model.matrix(~ sampleAnnot$isTarget + sampleAnnot$donorID)
  
  fit <- lmFit(expressionMatrix,designMatrix)
  fit <- eBayes(fit)
  
  limmaResults <- inner_join(as_tibble(fit$p.value[,"sampleAnnot$isTargetTRUE", drop=F], rownames="probe_name"), as_tibble(fit$t[,"sampleAnnot$isTargetTRUE", drop=F], rownames="probe_name"), by="probe_name")
  limmaResults %<>% rename( p.value=`sampleAnnot$isTargetTRUE.x`, t = `sampleAnnot$isTargetTRUE.y`)
  
  limmaResults <- inner_join(limmaResults, qc_table, by= "probe_name")
  #may have slight bias for longer genes with more probes
  gene_summary <- limmaResults %>% group_by(gene_symbol) %>% arrange(p.value) %>% 
    summarize(p.value = first(p.value), direction=sign(first(t)))
  #convert to ranks
  gene_summary %<>% mutate(pValueWithDirection = direction * (nrow(gene_summary) - rank(p.value)))

  #write out for target region
  if(grepl("habenula", targetRegion)) {
    print("Writing out")
    write_csv(limmaResults, paste0("./results/limma/", targetRegion,".",sourceExpression,".csv"))
    write_csv(gene_summary, paste0("./results/limma/", targetRegion,".",sourceExpression,".geneSummary.csv"))
  }
  
  gene_summary %<>% select(gene_symbol, !! targetRegion := pValueWithDirection)
  
  
  #join to a region by gene matrix
  if (is.null(regionByGene)) { 
    regionByGene <- gene_summary 
  } else { 
    regionByGene <- inner_join(regionByGene, gene_summary, by= "gene_symbol") 
  }
}

write_tsv(regionByGene, paste0("./data/processed/",sourceExpression, "_brainarea_vs_genes_exp_qcNames.tsv"))
