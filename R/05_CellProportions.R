library(tidyr)
library(readr)
library(dplyr)
library(magrittr)

#needs updating
cellProp <- read_csv("./results/Adult.Regional.cellEstimates.csv")
cellProp %<>% mutate_if(is.numeric, funs(. * -1)) %>% mutate_if(is.numeric, rank)

cellProp %<>% filter(grepl("habenula", structure_name_left_right_stripped))
cellProp <- melt(cellProp)
cellProp <- as_tibble(cellProp %>% arrange(value))
cellProp %<>% spread(structure_name_left_right_stripped, value)
#arrange by the lower of the two regions
cellProp %<>% mutate(min = pmin(`lateral habenular nucleus`, `medial habenular nucleus`)) %>% arrange(min) %>% select(-min)

write_csv(cellProp, "./results/Adult.Habenla.relative.cellEstimates.csv")

cellProp %>% filter(grepl("NeuroExpresso", variable)) %>% arrange(`lateral habenular nucleus`)
cellProp %>% filter(grepl("NeuroExpresso", variable)) %>% arrange(`medial habenular nucleus`)
cellProp %>% filter(grepl("Darm", variable))
