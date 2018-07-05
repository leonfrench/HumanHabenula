library(reshape2)
library(readr)
library(dplyr)
library(magrittr)
library(homologene)
library(markerGeneProfile)
sourceExpression <- "allen_HBA" #just use the adult data

#template code copied from Ogan Mancarci and the markerGeneProfile readme
probe2gene <- read_csv("./results/limma/lateral habenular nucleus.allen_HBA.csv") %>% select(probe_name, gene_symbol)
probeLevelExpression <- read_tsv("./data/processed/allen_HBA_sampleMatrix_qcNames.tsv")
probeLevelColAnnotations <- read_tsv("./data/processed/allen_HBA_sampleMatrix_colAnnotations_qcNames.tsv")

byRegionAllDonors = NULL #initialize to blank - does each donor individually
for (donor in unique(probeLevelColAnnotations$donorID)) {
  print(donor)
  samplesForDonor <- probeLevelColAnnotations %>% filter(donorID == donor)

  expressionMatrix <- probeLevelExpression %>% select(one_of(c("probe_name",samplesForDonor$uniqueID)))
  
  dim(expressionMatrix)
  
  #join
  expressionMatrix <- inner_join(probe2gene, expressionMatrix)
  expressionMatrix %<>% rename(Gene.Symbol = gene_symbol, Probe = probe_name)
  
  #setup gene lists
  otherGeneListsFolder <- "./data/other gene lists/"
  DarmanisLists <- list()
  for(geneListFilename in list.files(otherGeneListsFolder, pattern = "Darmanis.*txt", full.names = T)) {
    print(geneListFilename)
    genesOfInterest <- read.csv(geneListFilename,header=F,stringsAsFactors = F)
    shortName <- gsub(".txt","",gsub(paste0(".*/"),"", geneListFilename))
    DarmanisLists[shortName] <- list(genesOfInterest$V1)
  }
  
  ExpressoLists <- list()
  for(geneListFilename in list.files(otherGeneListsFolder, pattern = "NeuroExpresso.All.*txt", full.names = T)) {
    print(geneListFilename)
    genesOfInterest <- read.csv(geneListFilename,header=F,stringsAsFactors = F)$V1
    shortName <- gsub(".txt","",gsub(paste0(".*/"),"", geneListFilename))
    ExpressoLists[shortName] <- list(mouse2human(genesOfInterest)$humanGene)
  }
  
  medExp = expressionMatrix %>% 
    sepExpr() %>% {.[[2]]} %>%
    unlist %>% median
  
  expressionMatrix = mostVariable(expressionMatrix, threshold = medExp, threshFun=median)
  #this is code from Ogan 
  estimationsDarm =  mgpEstimate(exprData=expressionMatrix,
                                 genes=DarmanisLists,
                                 geneColName='Gene.Symbol',
                                 outlierSampleRemove=F, # should outlier samples removed. This is done using boxplot stats.
                                 geneTransform = function(x){ x }, # this is the default option for geneTransform
                                 groups=NULL, #if there are experimental groups provide them here. if not desired set to NULL
                                 seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
                                 removeMinority = TRUE) # removes genes if they are the minority in terms of rotation sign from estimation process
  
  estimationsExpresso =  mgpEstimate(exprData=expressionMatrix,
                                     genes=ExpressoLists,
                                     geneColName='Gene.Symbol',
                                     outlierSampleRemove=F, # should outlier samples removed. This is done using boxplot stats.
                                     geneTransform = function(x){ x }, # this is the default option for geneTransform
                                     groups=NULL, #if there are experimental groups provide them here. if not desired set to NULL
                                     seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
                                     removeMinority = TRUE) # removes genes if they are the minority in terms of rotation sign from estimation process
  
  estimationsDarm <- as.data.frame(estimationsDarm$simpleScaledEstimation)
  estimationsDarm$uniqueID <- rownames(estimationsDarm)

  estimationsExpresso <- estimationsExpresso$simpleScaledEstimation
  estimationsExpresso[names(estimationsExpresso[lapply(estimationsExpresso, length) == 0])] <- NULL
  estimationsExpresso <- as.data.frame(estimationsExpresso)
  estimationsExpresso$uniqueID <- rownames(estimationsExpresso)

  
  samplesForDonor <- inner_join(samplesForDonor, estimationsDarm)
  samplesForDonor <- inner_join(samplesForDonor, estimationsExpresso)
  
  #setup folders
  
  dir.create(paste0("./data/processed/", sourceExpression), showWarnings = FALSE)
  dir.create(paste0("./data/processed/", sourceExpression, "/", donor), showWarnings = FALSE)

  write_csv(samplesForDonor, paste0("./data/processed/", sourceExpression,"/",donor,"/SampleAnnot.cellEstimates.csv"))
  
  cellColumns <- intersect(colnames(samplesForDonor), union(names(DarmanisLists), names(ExpressoLists)))
  byRegion <- samplesForDonor %>% ungroup() %>% group_by(structure_name_left_right_stripped) %>% summarize_at(cellColumns, mean)
  byRegion %<>% mutate(donorID = donor)
  #bind them across brains
  byRegionAllDonors <- bind_rows(byRegionAllDonors, byRegion)
}

byRegionAllDonors


#from C8H10N4O2
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}
#NA's
byRegionAllDonors %>% group_by(donorID) %>% summarise_all(funs(sum(is.na(.)))) %>% melt() %>% filter(value > 0)
goodRegions <- as.character(byRegionAllDonors %>% summarise_all(funs(sum(is.na(.)))) %>% melt() %>% filter(value < 1) %>% .$variable)

#keep regions with no NA's
byRegionAllDonors %<>% select(one_of(goodRegions))

#normalize a column for each donor
byRegionAllDonorsScaled <- byRegionAllDonors %>% group_by(donorID) %>% mutate_if(is.numeric, scale_this) 
#average for each donor
byRegionAllDonorsScaled %<>% group_by(structure_name_left_right_stripped) %>% select(-donorID) %>% summarize_all(funs(mean(., na.rm = TRUE)))

#check for NA's again
byRegionAllDonorsScaled %>% summarise_all(funs(sum(is.na(.)))) %>% melt() %>% filter(value != 0)

write_csv(byRegionAllDonorsScaled, paste0("./results/Adult.Regional.cellEstimates.csv"))

