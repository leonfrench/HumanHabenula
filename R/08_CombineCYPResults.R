
allResults <- NULL
for (file in list.files(path = "./results/", pattern = "-CYP.tsv")) {
  result <- read_tsv(paste0("./results/", file)) %>% mutate(source = file)
  allResults <- bind_rows(result, allResults)
}
allResults %<>% tidyr::separate(source, sep="[-]", into=c("age", "region"))

allResults %>% group_by(ID) %>% summarize(maxSpec = max(specificity), maxCount = max(betterAUCCount))

allResults %>% group_by(ID, age) %>% summarize(maxSpec = max(specificity), maxCount = max(betterAUCCount))
allResults %>% group_by(ID, region) %>% summarize(maxSpec = max(specificity), maxCount = max(betterAUCCount))
