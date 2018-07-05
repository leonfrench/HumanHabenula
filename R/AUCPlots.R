library(plotROC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyr)
source("./R/GeneSetBuilders.R")

matrixAdult <- read_tsv("./data/processed/allen_HBA_brainarea_vs_genes_exp_qcNames.tsv", guess_max=22000)
matrixFetal <- read_tsv("./data/processed/allen_human_fetal_brain_brainarea_vs_genes_exp_qcNames.tsv", guess_max = 10000) 

matrixFetal %<>% select(gene_symbol, contains("habenula")) %>% 
  gather(key="region", value="rank", contains("habenula")) %>% 
  mutate(source = "Fetal")
matrixAdult %<>% select(gene_symbol, contains("habenula")) %>% 
  gather(key="region", value="rank", contains("habenula")) %>% 
  mutate(source = "Adult")
combined <- bind_rows(matrixAdult, matrixFetal)

#redo rank to end at max gene count (no negative values)
combined %<>% group_by(source, region) %>% mutate(rank = rank(rank)) %>% arrange(-rank)

#genesOfInterest <- loadCypGenes(unique(combined$gene_symbol))
#genesOfInterest <- genesOfInterest$MODULES2GENES$`CYP3A genes`

if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { #assume it's already loaded - needs a fix to see if the variable is declared
} else {
  geneSetsGO <- loadGOSets(allGenes)
}
#targetGroupID <- dplyr::filter(tbl_df(geneSetsGO$MODULES), Title == "MHC protein complex")$ID
targetGroupID <- dplyr::filter(tbl_df(geneSetsGO$MODULES), Title == "cilium movement")$ID
targetGroupID <- dplyr::filter(tbl_df(geneSetsGO$MODULES), Title == "cytokine receptor activity")$ID
targetGroupID <- dplyr::filter(tbl_df(geneSetsGO$MODULES), Title == "immunoglobulin binding")$ID


genesOfInterest <- unlist(geneSetsGO$MODULES2GENES[targetGroupID])

#geneSetsToUse <- loadFileSets("Extra")
#genesOfInterest <- unlist(geneSetsToUse$MODULES2GENES["hedonic.down"])



combined %<>% mutate(present = gene_symbol %in% genesOfInterest *1)
combined %<>% mutate(combinedLabel = paste(source, region))

cbbPalette <- c("#000000", "#E69F00", "#000000", "#E69F00")

(AUCPlot <- ggplot(combined, aes(d = present, m = rank, color=combinedLabel)) + ylab("") + 
    #geom_roc(n.cuts=0, linetype = "dashed") + 
    geom_roc(data = combined %>% filter(source == "Fetal"), n.cuts=0, linetype = "dotted") + 
    geom_roc(data = combined %>% filter(source == "Adult"), n.cuts=0, linetype = "solid") + 
    style_roc() + coord_cartesian(expand=F) +
    theme(legend.position = c(1,0), legend.justification = c(1, 0), legend.background= element_rect(fill = "transparent", colour = "transparent"), plot.margin=unit(c(.5,.5,.5,.5),"cm")) + 
    labs(color='Gene Group')  + 
    #facet_grid(dummy ~ ., switch="y")  + ylab("") +
    theme(strip.background = element_blank(), strip.placement = "inside", strip.text = element_blank()) +
    geom_abline(slope = 1, intercept = 0, colour = "grey90", size = 0.2) +
    scale_fill_manual(values=cbbPalette) + scale_colour_manual(values=cbbPalette) +
    guides(color = guide_legend(override.aes = list(linetype = c("solid","solid", "dotted","dotted")))) +
    theme(legend.key.width = unit(2, "line"))
    )

