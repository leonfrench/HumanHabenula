library(cowplot)
library(reshape2)
detach("package:dplyr", unload=TRUE)
library(dplyr)
library(magrittr)
library(readr)

allColeResults <- NULL
for (file in list.files(path = "./results/", pattern = "-Extra.tsv")) {
  result <- read_tsv(paste0("./results/", file)) %>% mutate(source = file)
  allColeResults <- bind_rows(result, allColeResults)
}
allColeResults %<>% tidyr::separate(source, sep="[-]", into=c("age", "region")) #this warning is fine

coleLists <- c("CTRA.down","eudonic.up","hedonic.down","hedonic.up","eudonic.down","CTRA.up")
allColeResults %<>% filter(MainTitle %in% coleLists) 
allColeResults %>% arrange(betterAUCCount)

allColeResults %<>% mutate(AUCSig = if_else(P.Value < 0.05, " *", "")) %>% mutate(AUCSig = paste0(signif(AUC, digits=3), AUCSig))

(summaryTable <- allColeResults %>% select(MainTitle, AUCSig, age, region) %>% mutate(location = paste(age, region)) %>% select(-age, -region) %>% spread(location, AUCSig) %>% arrange(desc(MainTitle)))
write_csv(summaryTable, path="./results/SummaryForCole.csv")
