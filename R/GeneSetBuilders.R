otherGeneListsFolder <- "./data/other gene lists/"
phenocarta_folder <- "./data/other gene lists/phenocarta/"

loadPhenocarta <- function(taxon, geneBackground) {
  #guess column types from the whole dataset, basically
  #phenocarta <- read_tsv(paste0(phenocarta_folder, "AllPhenocartaAnnotations.downloadedOct28.2016.tsv"), skip = 4, guess_max = 130000)
  phenocarta <- read_tsv(paste0(phenocarta_folder, "AllPhenocartaAnnotations.downloadedMay14.2018.tsv"), skip = 4, guess_max = 130000)
  phenocarta$ID <- gsub("http://purl.obolibrary.org/obo/", "", phenocarta$`Phenotype URIs`)
  phenocarta <- dplyr::filter(phenocarta, Taxon == taxon) %>% dplyr::select(symbol = `Gene Symbol`, name = `Phenotype Names`, ID) %>% filter(symbol %in% geneBackground) %>% distinct()
  geneLists <- phenocarta %>% group_by(ID) %>% summarize(name = paste(unique(name), collapse = ","), genes = unique(list(symbol)), size = dplyr::n()) %>% filter(size > 5 & size < 200) 
  #distinct(geneLists)
  namedLists <- geneLists$genes
  names(namedLists) <- geneLists$ID
  idToName <- data.frame(ID = geneLists$ID, Title = geneLists$name)
  geneSets <- makeTmod(modules = idToName, modules2genes = namedLists)
  geneSets
}

#cyp Genes
loadCypGenes <- function(geneBackground) {
  cypGenes <- geneBackground[grepl("^CYP[0-9]", geneBackground)]
  tmodNames <- data.frame()
  modules2genes <- list()
  modules2genes["CYP genes"] <- list(cypGenes)
  tmodNames <- data.frame(ID="CYP genes", Title = "CYP genes")
  
  cyp3Genes <- geneBackground[grepl("^CYP3A", geneBackground)]
  modules2genes["CYP3A genes"] <- list(cyp3Genes)
  tmodNames <- rbind(tmodNames, data.frame(ID="CYP3A genes", Title = "CYP3A genes"))

  modules2genes["CYP genes minus 3A"] <- list(setdiff(cypGenes, cyp3Genes))
  tmodNames <- rbind(tmodNames, data.frame(ID="CYP genes minus 3A", Title = "CYP genes minus 3A"))

  geneSetsCYP <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
  geneSetsCYP
}


loadGOSets <- function(geneBackground) {
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
    
    genesymbols <- intersect(genesymbols, geneBackground)
    if (!(length(genesymbols) >= 10 & length(genesymbols) <= 200)) next();
    
    modules2genes[goGroupName] <- list(genesymbols)
    
    tmodNames <- rbind(tmodNames, data.frame(ID=goGroupName, Title = Term(goGroupName)))
  }
  geneSetsGO <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
}

loadDisGeNetSets <- function(geneBackground) {
  disgenet <- read_tsv(paste0(otherGeneListsFolder, "/DisGeNET/curated_gene_disease_associations.tsv.gz"))
  disgenet %<>% dplyr::select(symbol = geneSymbol, name = diseaseName, ID = diseaseId)
  disgenet %<>% filter(symbol %in% geneBackground)
  
  geneLists <- group_by(disgenet, ID) %>% 
    dplyr::summarise(name = paste(unique(name), collapse = ","), genes = unique(list(symbol)), size = dplyr::n()) %>% 
    filter(size >= 10 & size <= 200) 
  namedLists <- geneLists$genes
  names(namedLists) <- geneLists$ID
  idToName <- data.frame(ID = geneLists$ID, Title = geneLists$name)
  geneSets <- makeTmod(modules = idToName, modules2genes = namedLists)
  geneSets
}

loadFileSets <- function(prefix) {
  tmodNames <- data.frame()
  modules2genes <- list()

  for(geneListFilename in list.files(otherGeneListsFolder, pattern = paste0(prefix,".*txt"), full.names = T)) {
    print(geneListFilename)
    genesOfInterest <- read.csv(geneListFilename,header=F,stringsAsFactors = F)
    shortName <- gsub(".txt","",gsub(paste0(".*/"),"", geneListFilename))
    shortName <- gsub(paste0(prefix,"."),"",shortName)
    
    genesOfInterest$term <- shortName
    
    modules2genes[shortName] <- list(genesOfInterest$V1)
    
    tmodNames <- rbind(tmodNames, data.frame(ID=shortName, Title = shortName))
  }
  geneSets <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
  geneSets
}
