sourceExpression <- "adult"
#sourceExpression <- "fetal"

#geneSetName <- "CYP" 
#geneSetName <- "GO" 
#geneSetName <- "PhenoCarta" 
#geneSetName <- "DisGeNet"
#geneSetName <- "Custom" 
geneSetName <- "Extra" 

targetGroup <- "positive regulation of hemopoiesis"
targetGroup <- "Drug-induced depressive state"
targetGroup <- "glycosaminoglycan binding"
targetGroup <- "regulation of odontogenesis of dentin-containing tooth"
targetGroup <- "icosanoid secretion"
targetGroup <- "leukocyte chemotaxis"
targetGroup <- "acetylcholine binding"
targetGroup <- "immunoglobulin binding"
targetGroup <- "neurotransmitter transporter activity"
targetGroup <- "positive regulation of hemopoiesis"
targetGroup <- "Uremia"
targetGroup <-"positive regulation of alcohol biosynthetic process"
targetGroup <-"FTO Obesity cluster"





if (sourceExpression == "adult") {
  medial <- read_csv("./results/limma/medial habenular nucleus.allen_HBA.geneSummary.csv.addedStats.csv", guess_max = 130000)
  lateral <- read_csv("./results/limma/lateral habenular nucleus.allen_HBA.geneSummary.csv.addedStats.csv", guess_max = 130000)
} else if (sourceExpression == "fetal") {
  medial <- read_csv("./results/limma/medial habenular nucleus.allen_human_fetal_brain.geneSummary.csv.addedStats.csv", guess_max = 130000)
  lateral <- read_csv("./results/limma/lateral habenular nucleus.allen_human_fetal_brain.geneSummary.csv.addedStats.csv", guess_max = 130000)
}

if (geneSetName == "PhenoCarta") {
  geneSetsToUse <- loadPhenocarta("human", allGenes)
} else if (geneSetName == "CYP") {
  geneSetsToUse <- loadCypGenes(allGenes)
} else if (geneSetName == "Custom") {
  geneSetsToUse <- loadFileSets(geneSetName)
} else if (geneSetName == "Extra") {
  geneSetsToUse <- loadFileSets(geneSetName)
} else if (geneSetName == "DisGeNet") {
  geneSetsToUse <- loadDisGeNetSets(allGenes)
} else if (geneSetName == "GO") {
  if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { #assume it's already loaded - needs a fix to see if the variable is declared
  } else {
    geneSetsGO <- loadGOSets(allGenes)
  }
  geneSetsToUse <- geneSetsGO
}

baseFilename <- paste0("./results/limma/filtered/", geneSetName)

targetGroupID <- dplyr::filter(tbl_df(geneSetsToUse$MODULES), Title == targetGroup)$ID
filter(medial, gene_symbol %in% unlist(geneSetsToUse$MODULES2GENES[targetGroupID])) %>% 
  dplyr::arrange(desc(pValueWithDirection)) %>% 
  write_tsv(paste0(baseFilename, ".medial.", targetGroup, ".",sourceExpression,".addedStats.tsv"))

filter(lateral, gene_symbol %in% unlist(geneSetsToUse$MODULES2GENES[targetGroupID])) %>% 
  dplyr::arrange(desc(pValueWithDirection)) %>% 
  write_tsv(paste0(baseFilename, ".lateral.", targetGroup,".",sourceExpression, ".addedStats.tsv"))

