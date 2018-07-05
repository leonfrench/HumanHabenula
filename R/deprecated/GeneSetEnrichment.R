#Leon is using GOSOURCEDATE: 2017-Mar29 (type GO.db to find out)

Sys.info()["nodename"]

filename <- "./results/limma/Habenula.3donors.lateral.allen_HBA.csv"
cat(paste("Working directory:", getwd()))
cat(paste("Using input file:",filename))
baseFilename <- gsub(".csv", "", filename)

otherGeneListsFolder <- "./data/other gene lists/"

library(ggsignif)
library(cowplot)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(magrittr)
library(homologene) #install via install_github('oganm/homologene')
library(org.Hs.eg.db)
library(AnnotationDbi)
library(annotate)
library(GO.db)
library(tmod)
library(metap)
library(reshape2)
library(xlsx)

geneStatistics <- read_csv(filename) 

#update gene symbols
symbol_to_entrez <- geneStatistics %>% dplyr::select(probe_name, allenEntrez) %>% filter(!is.nan(allenEntrez)) %>% filter(!is.na(allenEntrez)) %>% distinct()
symbol_to_entrez %<>% mutate(new_symbol = getSYMBOL(as.character(allenEntrez), data='org.Hs.eg')) 

geneStatistics <- inner_join(geneStatistics, symbol_to_entrez)

geneStatistics %<>% mutate(pValueWithDirection = if_else(t > 0, nrow(geneStatistics) - rank(p.value), -1* nrow(geneStatistics) + rank(p.value)))


paste("Genes with negative correlations:", dplyr::filter(geneStatistics, p.adj < 0.05, t < 0) %>% summarize(n = n()))
paste("Genes with positive correlations:", dplyr::filter(geneStatistics, p.adj < 0.05, t > 0) %>% summarize(n = n()))

geneStatistics %<>% dplyr::select(-Gene_symbol) %>% dplyr::rename(geneSymbol = new_symbol) 
geneStatistics %<>% dplyr::select(probe_name, geneSymbol, everything())
#add gene names
geneStatistics %<>% mutate(name = unlist(lookUp(as.character(allenEntrez), "org.Hs.eg", "GENENAME"))) %>% dplyr::select(geneSymbol, name, everything())

#sort 
geneStatistics <- arrange(geneStatistics, desc(pValueWithDirection))

write_csv(geneStatistics, paste0(baseFilename, ".addedStats.csv"))

sortedGenes <- geneStatistics$geneSymbol

######################################################
# AUC via tmod
######################################################

if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { #assume it's already loaded - needs a fix to see if the variable is declared
} else {
  go_object <- as.list(org.Hs.egGO2ALLEGS)
  
  symbolsInGO <- getSYMBOL(unique(unlist(go_object)), data='org.Hs.eg')
  
  #build GO sets for tmod -slow
  tmodNames <- data.frame()
  modules2genes <- list()
  goGroupName <- names(go_object)[1]
  showMethods(Term)
  
  goCount <- length(go_object)
  count <- 1
  for(goGroupName in names(go_object)) {
    if (count %% 1000 == 0) print(paste(count, "of", goCount))
    count <- count + 1
    
    goGroup <- go_object[goGroupName]
    geneIDs <- unique(unlist(goGroup, use.names=F))  #discard evidence codes
    genesymbols <- unique(getSYMBOL(geneIDs, data='org.Hs.eg'))
    
    genesymbols <- intersect(genesymbols, sortedGenes)
    if (!(length(genesymbols) >= 10 & length(genesymbols) <= 200)) next();
    
    modules2genes[goGroupName] <- list(genesymbols)
    
    tmodNames <- rbind(tmodNames, data.frame(ID=goGroupName, Title = Term(goGroupName)))
  }
  geneSetsGO <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
}


result <- tbl_df(tmodUtest(c(sortedGenes), mset=geneSetsGO, qval = 1, filter = T))
result %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% dplyr::select(-adj.P.Val) 
result %<>% rowwise() %>% mutate(aspect = Ontology(ID))
(result %<>% group_by(U, N1, AUC, P.Value) %>% summarize(MainTitle = first(Title),  ID=paste(ID, collapse=","), aspect= first(aspect), allNames = if_else(n() > 1, paste(Title[2:length(Title)], collapse=","), "")))
result %<>% ungroup() %>% mutate(adj.P.Value=p.adjust(P.Value, method = "fdr")) %>% arrange(P.Value) #adjust again
result$rank <- 1:nrow(result)
result %<>% dplyr::select(MainTitle, geneCount = N1, AUC, P.Value, adj.P.Value, everything(), -U) 

sum(result$adj.P.Value < 0.05)

head(filter(result, AUC > 0.5),n=15)

result$adj.P.Value <- signif(result$adj.P.Value, digits=3)
result$AUC <- signif(result$AUC, digits=3)
result$P.Value <- signif(result$P.Value, digits=3)

# Output tables for top ten positively and negatively enriched GO groups for paper
(result.up <- head(filter(result, AUC > 0.5) %>% dplyr::select(-ID), n=10))

#check one p-value to test TMOD
write_csv( result,  paste(baseFilename,".GO.results.csv",sep=""))
write_csv( dplyr::select(result.up, Name = MainTitle,`Gene Count` = geneCount, AUROC = AUC,  `p` = P.Value, `pFDR` = adj.P.Value, aspect),  paste0(baseFilename,".GO.up10.csv"))


##### Look at what proportion of the results match certain categories

cat(paste("Number of significant GO groups",as.character(lengths(result %>% filter (adj.P.Value < 0.05))[1])))

cat(paste("Number of significant positively enriched GO groups",as.character(lengths(result %>% filter (adj.P.Value < 0.05, AUC > 0.5))[1])))

cat(paste("Number of significant negatively enriched GO groups",as.character(lengths(result %>% filter (adj.P.Value < 0.05, AUC < 0.5))[1])))

source("./R Code/ROCPlots.R")

dplyr::filter(tbl_df(geneSetsGO$MODULES), Title == "detection of chemical stimulus involved in sensory perception")$ID
dplyr::filter(tbl_df(geneSetsGO$MODULES), Title == "transcriptional repressor activity, RNA polymerase II transcription regulatory region sequence-specific binding")$ID

#plots <- createPlots(sortedGenes, c("GO:0098798", "GO:1905368", "GO:0005882","GO:0031490","GO:0005814","GO:0033038","GO:0032452", "GO:0001227", "GO:0045178", "GO:0001047"), geneSetsGO)

plots <- createPlots(sortedGenes, c("GO:0050907", "GO:0098798", "GO:1905368", "GO:0033038","GO:0006521"), geneSetsGO)
(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.6),scale = 0.95,labels = c("A", "B"))) #add labels = c("A", "B"), for manuscript
#save as 11x11 PDF

#######################
# hypergeometric testing
########################
result <- tmodHGtest(fg= head(sortedGenes,n=50), bg=sortedGenes, mset=geneSetsGO,qval = 1.1, filter=T)

result %<>% rowwise() %>% mutate(aspect = Ontology(ID))
(result %<>% group_by(U, N1, AUC, P.Value) %>% summarize(MainTitle = first(Title),  ID=paste(ID, collapse=","), aspect= first(aspect), allNames = if_else(n() > 1, paste(Title[2:length(Title)], collapse=","), "")))
result %<>% ungroup() %>% mutate(adj.P.Value=p.adjust(P.Value, method = "fdr")) %>% arrange(P.Value) #adjust again
result$rank <- 1:nrow(result)
result %<>% dplyr::select(MainTitle, geneCount = N1, AUC, P.Value, adj.P.Value, everything(), -U) 



#################################################################
#################################################################
tmodNames <- data.frame()
modules2genes <- list()


for(geneListFilename in list.files(otherGeneListsFolder, pattern = ".*txt", full.names = T)) {
  print(geneListFilename)
  genesOfInterest <- read.csv(geneListFilename,header=F,stringsAsFactors = F)
  shortName <- gsub(".txt","",gsub(paste0(".*/"),"", geneListFilename))
  
  genesOfInterest$term <- shortName
  
  #already a human gene list
  if (grepl(pattern = "Darmanis.", geneListFilename  ) | grepl(pattern = "Mistry", geneListFilename) | grepl(pattern = "HouseKeeping", geneListFilename  ) | grepl(pattern = "human", geneListFilename  )) {
    modules2genes[shortName] <- list(genesOfInterest$V1)
  } else { #needs conversion from mouse
    print(" converting from mouse to human")
    modules2genes[shortName] <- list(mouse2human(genesOfInterest$V1)$humanGene)
  }
  
  tmodNames <- rbind(tmodNames, data.frame(ID=shortName, Title = shortName))
}
geneSets <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
geneSetsCellType <- geneSets #for later reuse

result <- tmodUtest(sortedGenes, mset=geneSets, qval = 1, filter = F)
result <- tbl_df(result) %>% dplyr::select(Title, geneCount =N1,AUC,  P.Value, adj.P.Val, -ID)
result %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr")) #tmod runs one-sided tests

result$adj.P.Val <- signif(result$adj.P.Val, digits=3)
result$AUC <- signif(result$AUC, digits=3)

writeTableWrapper <- function(prefixFilter, result) {
  subsetResult <- filter(result, grepl(prefixFilter, Title))
  subsetResult$oldTitle <- subsetResult$Title
  subsetResult$Title <- gsub(paste0(prefixFilter, "."), "", subsetResult$Title)
  subsetResult$Title <- gsub("[.]", " ", subsetResult$Title)
  subsetResult$Title <- gsub("[_]", " ", subsetResult$Title)
  subsetResult$Title <- gsub("Neuron interneuron", "Interneuron", subsetResult$Title)
  subsetResult$Title <- gsub("ligo", "ligodendrocyte", subsetResult$Title)
  subsetResult$adj.P.Val <- p.adjust(subsetResult$P.Value, method="fdr")
  subsetResult$adj.P.Val <- signif(subsetResult$adj.P.Val, digits=3)
  write_csv(dplyr::select(subsetResult, `Cell-type or class` = Title,`Gene Count` = geneCount, AUROC = AUC,  `pFDR` = adj.P.Val), paste0(baseFilename,".", prefixFilter,".csv")) 
  subsetResult
}

(zeiselResult <- writeTableWrapper("Zeisel", result))
(darmResult <- writeTableWrapper("Darmanis", result))
writeTableWrapper("NeuroExpresso.Cortex", result)

writeTableWrapper("Thakurela", result)
writeTableWrapper("Custom", result)


plots <- createPlots(sortedGenes, as.character(zeiselResult$oldTitle), geneSets, customNames=zeiselResult$Title, filter=F)
(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.9),scale = 0.95,labels = c("A", "B"))) #add labels = c("A", "B"), for manuscript
#save as 11x11


plots <- createPlots(sortedGenes, as.character(darmResult$oldTitle), geneSets, customNames=darmResult$Title, filter=F)
(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.8), labels = c("A", "B"))) #add labels = c("A", "B"), for manuscript
plots$rasterPlot #save as 10x4 pdf
#filter(geneStatistics, geneSymbol %in% geneSets["Darmanis.Oligo"]$GENES$ID)





