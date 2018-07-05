library(cowplot)
library(reshape2)
detach("package:dplyr", unload=TRUE)
library(dplyr)
library(magrittr)
library(readr)

allGOResults <- NULL
for (file in list.files(path = "./results/", pattern = "-GO.tsv")) {
  result <- read_tsv(paste0("./results/", file)) %>% mutate(source = file)
  allGOResults <- bind_rows(result, allGOResults)
}
allGOResults %<>% tidyr::separate(source, sep="[-]", into=c("age", "region")) #this warning is fine

#check for non matches
allGOResults %>% group_by(ID) %>% summarize(n=dplyr::n()) %>% filter(n!=4)

#unique groups:
length(unique(allGOResults$ID))

#how many significant?
allGOResults %>% filter(age == "adult") %>% group_by(region) %>% filter(adj.P.Value < 0.05) %>% summarize(n=dplyr::n())

#how many significant and best AUC
allGOResults %>% filter(age == "adult") %>% group_by(region) %>% filter(betterAUCCount == 0) %>% summarize(n=dplyr::n())

nicotineGroups <- allGOResults %>% filter(grepl( "nicotine", MainTitle)) %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Val = p.adjust(P.Value, method="fdr"))
nicotineGroups %<>% filter(age=="adult")
nicotineGroups %>% ungroup() %>% select(MainTitle) %>% distinct()
nicotineGroups %>% filter(adj.P.Val < 0.05)

alcoholGroups <- allGOResults %>% filter(grepl( "alcohol", MainTitle)) %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Val = p.adjust(P.Value, method="fdr"))
alcoholGroups %<>% filter(age=="adult")
alcoholGroups %>% ungroup() %>% select(MainTitle) %>% distinct()
alcoholGroups %>% filter(adj.P.Val < 0.05) 
alcoholGroups %>% filter(MainTitle == "response to alcohol")
alcoholGroups %>% filter(adj.P.Val < 0.05) %>% select(MainTitle, ID) #look at which genes

allGOResults %>% filter(grepl( "cocaine", MainTitle)) %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Val = p.adjust(P.Value, method="fdr"))
allGOResults %>% filter(grepl( "morphine", allNames)) %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Val = p.adjust(P.Value, method="fdr"))
allGOResults %>% filter(grepl( "amphet", MainTitle)) %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Val = p.adjust(P.Value, method="fdr"))

allGOResults %>% filter(MainTitle == "response to xenobiotic stimulus") %>% arrange(P.Value) %>% filter(age=="adult")
allGOResults %>% filter(MainTitle == "xenobiotic metabolic process") %>% arrange(P.Value) %>% filter(age=="adult")

lateralTop15 <- allGOResults %>% filter(age=="adult", region == "lateral habenular nucleus") %>% arrange(P.Value) %>% head(15) %>% select(-age, -region)
medialTop15 <- allGOResults %>% filter(age=="adult", region == "medial habenular nucleus") %>% arrange(P.Value) %>% head(15) %>% select(-age, -region)

intersect(lateralTop15$MainTitle, medialTop15$MainTitle)

medialTop15 %>% mutate(AUC = signif(AUC, digits=3), adj.P.Value = signif(adj.P.Value, digits=3)) %>%
  select(Name = MainTitle, Genes = geneCount, AUROC=AUC, `Specificity Rank` = betterAUCCount, `p-valueFDR` = adj.P.Value) %>%
  write_tsv("./results/adult-medial habenular nucleus-GO.top15.tsv")

lateralTop15 %>% mutate(AUC = signif(AUC, digits=3), adj.P.Value = signif(adj.P.Value, digits=3)) %>%
  select(Name = MainTitle, Genes = geneCount, AUROC=AUC, `Specificity Rank` = betterAUCCount, `p-valueFDR` = adj.P.Value) %>%
  write_tsv("./results/adult-lateral habenular nucleus-GO.top15.tsv")

#filter for specific and significant in fetal
sigAndSpec <- allGOResults %>% group_by(ID, region) %>% filter(age=="adult", adj.P.Value < 0.05, betterAUCCount <= 1) %>% arrange(rank)
sigAndSpec %>% group_by(region) %>% summarize(n=dplyr::n()) 
  
fetalEnriched <- allGOResults %>% filter(age=="fetal", P.Value < 0.05) %>% select(ID, region, fetalrank = rank, fetalAUC= AUC, fetalP = P.Value)
inner_join(sigAndSpec, fetalEnriched, by = c("ID", "region")) %>% group_by(region) %>% summarize(n=dplyr::n()) 

medialSigAndSpec <- inner_join(sigAndSpec, fetalEnriched, by = c("ID", "region")) %>% filter(region=="medial habenular nucleus") 

#cut into shorter list for lateral (AUC count == 0)
lateralSigAndSpec <- inner_join(sigAndSpec, fetalEnriched, by = c("ID", "region")) %>% filter(region=="lateral habenular nucleus", betterAUCCount==0)


lateralSigAndSpec %>% ungroup() %>% mutate(AUC = signif(AUC, digits=3), adj.P.Value = signif(adj.P.Value, digits=3), fetalAUC = signif(fetalAUC, digits=3), fetalP = signif(fetalP, digits=3)) %>%
  select(Name = MainTitle, Genes = geneCount, AUROC=AUC, `Specificity Rank` = betterAUCCount, `p-valueFDR` = adj.P.Value, `Fetal AUROC` = fetalAUC,  `Fetal p-value` = fetalP) %>%
  write_tsv("./results/adult-lateral habenular nucleus-GO.SigAndSpec.tsv")
#lateral results
binom.test(nrow(lateralSigAndSpec), nrow(sigAndSpec %>% filter(region=="lateral habenular nucleus", betterAUCCount==0)), p = 0.05, alternative = c("greater"), conf.level = 0.95)

medialSigAndSpec %>% ungroup() %>% mutate(AUC = signif(AUC, digits=3), adj.P.Value = signif(adj.P.Value, digits=3), fetalAUC = signif(fetalAUC, digits=3), fetalP = signif(fetalP, digits=3)) %>%
  select(Name = MainTitle, Genes = geneCount, AUROC=AUC, `Specificity Rank` = betterAUCCount, `p-valueFDR` = adj.P.Value, `Fetal AUROC` = fetalAUC,  `Fetal p-value` = fetalP) %>%
  write_tsv("./results/adult-medial habenular nucleus-GO.SigAndSpec.tsv")
binom.test(nrow(medialSigAndSpec), nrow(sigAndSpec %>% filter(region=="medial habenular nucleus")), p = 0.05, alternative = c("greater"), conf.level = 0.95)


#filter on better hits only
allGOResults %>% group_by(MainTitle, region) %>% filter(betterAUCCount < 3) %>% summarize(n = dplyr::n()) %>% filter(n>1)

allGOResults %>% group_by(MainTitle, region) %>% filter(betterAUCCount < 4) %>% summarize(n = dplyr::n()) %>% filter(n>1)



