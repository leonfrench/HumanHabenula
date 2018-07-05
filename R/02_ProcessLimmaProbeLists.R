library(org.Hs.eg.db)
library(annotate)
library(readr)
library(magrittr)
detach("package:dplyr", unload=TRUE)
library(dplyr)

#target <- "lateral"
target <- "medial"

adultFilename <- paste0("./results/limma/",target," habenular nucleus.allen_HBA.geneSummary.csv")
fetalFilename <- paste0("./results/limma/",target," habenular nucleus.allen_human_fetal_brain.geneSummary.csv")
#load fetal and adult data
adult <- read_csv(adultFilename)
fetal <- read_csv(fetalFilename)

adultProbes <- read_csv(paste0("./results/limma/",target," habenular nucleus.allen_HBA.csv"))
fetalProbes <- read_csv(paste0("./results/limma/",target," habenular nucleus.allen_human_fetal_brain.csv"))

#check counts
length(unique(adultProbes$gene_symbol))
length(unique(fetalProbes$gene_symbol))
length(intersect(adultProbes$gene_symbol, fetalProbes$gene_symbol))
length(unique(adultProbes$probe_name))
length(unique(fetalProbes$probe_name))
length(intersect(adultProbes$probe_name, fetalProbes$probe_name))

adultProbes %<>% mutate(adj.p.value = p.adjust(p.value, method="fdr")) %>% arrange(p.value) %>% 
  mutate(sigUpRegulated = adj.p.value < 0.05 & t > 0) %>% print()
fetalProbes %<>% mutate(adj.p.value = p.adjust(p.value, method="fdr")) %>% arrange(p.value) %>% 
  mutate(sigUpRegulated = adj.p.value < 0.05 & t > 0) %>% print()

supplementTable <- inner_join(
  adultProbes %>% select(probe_name, gene_symbol, p.value, t, adj.p.value), 
  fetalProbes %>% select(probe_name, p.value, t, adj.p.value), suffix = c(".adult", ".fetal")
  , by="probe_name")

paste("Significantly enriched probes for adult", target, ":", nrow(adultProbes %>% filter(sigUpRegulated)))
paste("Significantly enriched probes for fetal", target, ":", nrow(fetalProbes %>% filter(sigUpRegulated)))

sigAdultGenes <- unique(adultProbes %>% filter(sigUpRegulated) %>% .$gene_symbol)
sigFetalGenes <- unique(fetalProbes %>% filter(sigUpRegulated) %>% .$gene_symbol)

paste("Significantly enriched genes for adult", target, ":",  length(sigAdultGenes))
paste("Significantly enriched genes for fetal", target, ":",  length(sigFetalGenes))
paste("Intersecting genes", target, ":",  length(intersect(sigAdultGenes, sigFetalGenes)))


#get the probe level corrected max p-value
threshold <- max((adultProbes %>% filter(sigUpRegulated))$p.value)
adult %<>% arrange(-pValueWithDirection) %>% mutate(is_probe_corrected_significant = p.value <= threshold & direction > 0) %>% print()
threshold <- max((fetalProbes %>% filter(sigUpRegulated))$p.value)
fetal %<>% arrange(-pValueWithDirection) %>% mutate(is_probe_corrected_significant = p.value <= threshold & direction > 0) %>% print()

adult %<>% mutate(rank = rank(-pValueWithDirection))
fetal %<>% mutate(rank = rank(-pValueWithDirection))

geneToName <- adultProbes %>% select(gene_symbol, entrez_id) %>% distinct() %>% filter(!is.na(entrez_id))

#add gene names
geneToName %<>% mutate(name = unlist(lookUp(as.character(entrez_id), "org.Hs.eg", "GENENAME"))) %>% distinct() %>% print()

adult <- left_join(adult, geneToName) %>% select(gene_symbol, name, entrez_id, everything())
fetal <- left_join(fetal, geneToName) %>% select(gene_symbol, name, entrez_id, everything())
supplementTable <- left_join(supplementTable, geneToName) %>% select(gene_symbol, name, entrez_id, everything())

#adult %>% group_by(gene_symbol) %>% summarize(n=n()) %>% filter(n >1)

#top 20
#symbol, name, number of significant probes, p-value, fetal rank
sigProbeCount <- adultProbes %>% filter(sigUpRegulated) %>% group_by(gene_symbol) %>% summarize(sigProbes = dplyr::n())
adult20 <- adult %>% arrange(rank) %>% head(20)
adult20 %<>% inner_join(sigProbeCount)
adult20 %<>% mutate(p.value = signif(p.value, digits=3)) %>% 
  select(`Gene Symbol` = gene_symbol, `Name` = name, `Significant probes` = sigProbes, `p-value` = p.value) %>% 
  print()

write_csv(adult, paste0(adultFilename, ".addedStats.csv") )
write_csv(fetal, paste0(fetalFilename, ".addedStats.csv") )
write_csv(supplementTable, paste0(adultFilename, ".supplementTable.csv") )
write_csv(adult20, paste0(adultFilename, ".top20.csv") )

#show overlapping
#compare medial and lateral top 20
adult20 <- adult %>% arrange(rank) %>% head(20)
inner_join(adult20, fetal %>% select(gene_symbol, fetal_rank = rank)) %>% filter(fetal_rank < 21)
inner_join(adult20, fetal %>% select(gene_symbol, fetal_rank = rank)) %>% filter(fetal_rank < 21) %>% .$gene_symbol

