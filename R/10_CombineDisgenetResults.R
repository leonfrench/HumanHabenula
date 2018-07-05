library(cowplot)
library(reshape2)
detach("package:dplyr", unload=TRUE)
library(dplyr)
library(magrittr)
library(readr)

allDisResults <- NULL
for (file in list.files(path = "./results/", pattern = "-DisGeNet.tsv")) {
  result <- read_tsv(paste0("./results/", file)) %>% mutate(source = file)
  allDisResults <- bind_rows(result, allDisResults)
}
allDisResults %<>% tidyr::separate(source, sep="[-]", into=c("age", "region")) #this warning is fine

#check for non matches
allDisResults %>% group_by(ID) %>% summarize(n=dplyr::n()) %>% filter(n!=4)

#unique groups:
length(unique(allDisResults$ID))

#how many significant?
allDisResults %>% filter(age == "adult") %>% group_by(region) %>% filter(adj.P.Value < 0.05) %>% summarize(n=dplyr::n())

allDisResults %>% filter(age == "adult") %>% group_by(region) %>% filter(betterAUCCount == 0) %>% summarize(n=dplyr::n())

alcoholGroups <- allDisResults %>% filter(grepl( "[Aa]lcohol", MainTitle)) %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Val = p.adjust(P.Value, method="fdr"))
alcoholGroups %<>% filter(age=="adult") %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr"))
alcoholGroups %>% ungroup() %>% select(MainTitle) %>% distinct()
alcoholGroups %>% arrange(P.Value) %>% select(MainTitle, region, everything(), -otherNames, -otherNames, -aspect)

allDisResults %>% filter(grepl( "cocaine", MainTitle, ignore.case = T)) %>% filter(age=="adult") %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr"))
allDisResults %>% filter(grepl( "amphet", MainTitle, ignore.case = T)) %>% filter(age=="adult") %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr"))
allDisResults %>% filter(grepl( "Marijuana", MainTitle, ignore.case = T)) %>% filter(age=="adult") %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr"))

depressionGroups <- allDisResults %>% filter(MainTitle %in% c("Depressive Symptoms", "Depressed mood", "Drug-induced depressive state", "Depression, Bipolar")) %>% arrange(P.Value) %>% group_by(age, region) %>% mutate(adj.P.Val = p.adjust(P.Value, method="fdr"))
depressionGroups %<>% filter(age=="adult") %>% group_by(age, region) %>% mutate(adj.P.Value = p.adjust(P.Value, method="fdr"))
depressionGroups %>% ungroup() %>% select(MainTitle) %>% distinct()
depressionGroups %>% arrange(P.Value) %>% select(MainTitle, region, everything(), -aspect, -otherNames)


lateralTop10 <- allDisResults %>% filter(age=="adult", region == "lateral habenular nucleus") %>% arrange(P.Value) %>% head(10) %>% select(-age, -region)
medialTop10 <- allDisResults %>% filter(age=="adult", region == "medial habenular nucleus") %>% arrange(P.Value) %>% head(10) %>% select(-age, -region)

intersect(lateralTop10$MainTitle, medialTop10$MainTitle)

#dont write out medial
medialTop10 %>% mutate(AUC = signif(AUC, digits=3), adj.P.Value = signif(adj.P.Value, digits=3)) %>%
  select(Name = MainTitle, Genes = geneCount, AUROC=AUC, `Specificity Rank` = betterAUCCount, `p-valueFDR` = adj.P.Value)

lateralTop10 %>% mutate(AUC = signif(AUC, digits=3), adj.P.Value = signif(adj.P.Value, digits=3)) %>%
  select(Name = MainTitle, Genes = geneCount, AUROC=AUC, `Specificity Rank` = betterAUCCount, `p-valueFDR` = adj.P.Value) %>%
  write_tsv("./results/adult-lateral habenular nucleus-DisGeNet.top10.tsv")
