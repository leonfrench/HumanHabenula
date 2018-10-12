library(cowplot)
library(reshape2)
library(dplyr)
library(magrittr)
library(cowplot)
library(ggplot2)
library(readr)

forPoster<-FALSE
#load fetal and human - files written by limma.allRegions.R
expressionMatrixAdult <- read_tsv("./data/processed/allen_HBA_sampleMatrix_qcNames.tsv")
expressionMatrixFetal <- read_tsv("./data/processed/allen_human_fetal_brain_sampleMatrix_qcNames.tsv")

sampleAnnotAdult <- read_tsv(paste0("./data/processed/allen_HBA_sampleMatrix_colAnnotations_qcNames.tsv"))
sampleAnnotFetal <- read_tsv(paste0("./data/processed/allen_human_fetal_brain_sampleMatrix_colAnnotations_qcNames.tsv"))

#grab it from the limma results file
probesInfo <- read_csv(paste0("./results/limma/lateral habenular nucleus.allen_HBA.csv")) %>% select(-t, -p.value, -is_qc_pass)
genes_of_interest <- c("GPR151")
#genes_of_interest <- c("CYP3A4", "CYP3A5")
#genes_of_interest <- c("UGT2B15")
#genes_of_interest <- c("TRPM8")


#genes_of_interest <- c("CYP3A7")

#genes_of_interest <- c("CYP3A4", "CYP3A43", "CYP3A5","CYP3A7")
#genes_of_interest <- c("CYP3A4", "CYP3A5","CYP3A7")
#genes_of_interest <- c("GPR151", "CYP3A4", "CYP3A5","CHRNB3","TRPM8")





#genes_of_interest <- c("UGT2B15")

#poster todo - plot only significant probes - show p-values? set order

#filter for probes for gene of interest - where is the probe info? 
probesInfo %<>% filter(gene_symbol %in% genes_of_interest)

fetalFocus <- expressionMatrixFetal %>% filter(probe_name %in% probesInfo$probe_name)
adultFocus <- expressionMatrixAdult %>% filter(probe_name %in% probesInfo$probe_name)

fetalFocus <- as_tibble(melt(fetalFocus, value.name = "Expression", variable.name = "uniqueID"))
adultFocus <- as_tibble(melt(adultFocus, value.name = "Expression", variable.name = "uniqueID"))

if (forPoster) {
  adultFocus <- inner_join(adultFocus, probesInfo) 
  fetalFocus <- inner_join(fetalFocus, probesInfo) 
  fetalFocus %>% select(probe_name, gene_symbol) %>% distinct()
  probesToUse <- c("A_23_P8801", "A_24_P389251", "CUST_189_PI416573500", "CUST_194_PI416379584", "CUST_13913_PI416261804", "CUST_13911_PI416261804")
  fetalFocus %<>% filter(probe_name %in% probesToUse)
  adultFocus %<>% filter(probe_name %in% probesToUse)
  nameTable <- fetalFocus %>% select(probe_name, gene_symbol) %>% distinct() %>% rename(name = gene_symbol)
  nameTable %<>% mutate(name = if_else(probe_name == "CUST_189_PI416573500", "GPR151 (probe 1)", name))
  nameTable %<>% mutate(name = if_else(probe_name == "CUST_194_PI416379584", "GPR151 (probe 2)", name))
  adultFocus <- inner_join(adultFocus , nameTable) %>% mutate(probe_name = name)
  fetalFocus <- inner_join(fetalFocus , nameTable) %>% mutate(probe_name = name)
  adultFocus$probe_name <- factor(adultFocus$probe_name, levels= c("GPR151 (probe 1)", "GPR151 (probe 2)","CHRNB3", "CYP3A4","CYP3A5","TRPM8"))
  fetalFocus$probe_name <- factor(fetalFocus$probe_name, levels= c("GPR151 (probe 1)", "GPR151 (probe 2)","CHRNB3", "CYP3A4","CYP3A5","TRPM8"))
} else if (length(genes_of_interest > 1)) {
  probesInfo %<>% mutate(joinedName = paste0(probe_name, " (", gene_symbol, ")"))
  adultFocus <- inner_join(adultFocus, probesInfo) %>% mutate(probe_name = joinedName)
  fetalFocus <- inner_join(fetalFocus, probesInfo) %>% mutate(probe_name = joinedName)
  fetalFocus$probe_name <- factor(fetalFocus$probe_name, levels = unique((fetalFocus %>% arrange(gene_symbol))$probe_name))
  adultFocus$probe_name <- factor(adultFocus$probe_name, levels = unique((adultFocus %>% arrange(gene_symbol))$probe_name))
}



fetalFocus %<>% inner_join(sampleAnnotFetal %>% select(uniqueID, Region = structure_name_left_right_stripped, Donor=donorID))
adultFocus %<>% inner_join(sampleAnnotAdult %>% select(uniqueID, Region = structure_name_left_right_stripped, Donor=donorID))

adultFocus %<>% mutate(Donor = gsub("normalized_microarray_donor", "", Donor))
fetalFocus %<>% mutate(Donor = gsub("lmd_matrix_", "", Donor))

adultFocus %<>% mutate(Region = if_else(grepl("habenula", Region), Region, "remaining structures"))
fetalFocus %<>% mutate(Region = if_else(grepl("habenula", Region), Region, "remaining structures"))

#are they seperable?
adultFocus %>% group_by(probe_name, Region) %>% summarize(maxE = max(Expression))

adultFocus %<>% mutate(Donor = paste0("ID", Donor))
fetalFocus %<>% mutate(Donor = paste0("ID", Donor))

#plot - medial and lateral habenula, then all other Regions, different symbols for probes? 
ggplot(adultFocus, aes(y=Expression, x=probe_name, color=Region)) + geom_point() + facet_wrap( ~Donor)
ggplot(fetalFocus, aes(y=Expression, x=probe_name, color=Region)) + geom_point() + facet_wrap( ~Donor)

#color bind freindly palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

(adultPlot <- ggplot(adultFocus, aes(y=Expression, x=Donor, fill = Region, color=Region)) + 
  geom_violin(data = filter(adultFocus, Region == "remaining structures")) + theme_bw() +
    geom_point(data = filter(adultFocus, Region != "remaining structures")) + 
    facet_wrap( ~probe_name, nrow=1) + #, scales="free_y", nrow=1) +
  scale_fill_manual(values=cbbPalette) + scale_colour_manual(values=cbbPalette) + xlab("") + guides(fill=FALSE, color=FALSE) 
)

(fetalPlot <- ggplot(fetalFocus, aes(y=Expression, x=Donor, fill = Region, color=Region)) + 
  geom_violin(data = filter(fetalFocus, Region == "remaining structures")) + theme_bw() +
  geom_point(data = filter(fetalFocus, Region != "remaining structures")) + 
  facet_wrap( ~probe_name, nrow=1) + #, scales="free_y", nrow=1) +
  scale_fill_manual(values=cbbPalette) + scale_colour_manual(values=cbbPalette) + theme(legend.position="bottom"))

fetalPlot <- fetalPlot + theme(plot.margin = margin(.85, 0.1, 0, 0.1, "cm"))
adultPlot <- adultPlot + theme(plot.margin = margin(.85, 0.1, 0, 0.1, "cm"))

plot_grid(adultPlot, fetalPlot, rel_heights = c(8.5,10), labels = c("A) Adult brain", "B) Fetal brain"), nrow=2, hjust=0)
#save as for GPR151 7x6.6 PDF
#for CYP 11x6.6

if(forPoster) {
  plot_grid(adultPlot + theme(strip.text.x = element_text(size = 14)), 
            fetalPlot+theme(strip.text.x = element_text(size = 14)), #+ guides(fill=FALSE, color=FALSE) , 
            rel_heights = c(8.5,10), labels = c("Adult brain", "Fetal brain"), nrow=2, label_x=0.5, hjust=0.5)
  
  #Tyler Rinker and Luciano Selzer via stackoverflow
  g_legend <- function(a.gplot){ 
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    return(legend)} 
  legend <- g_legend(fetalPlot+theme(strip.text.x = element_text(size = 14)) + theme(legend.position="right")) 
  #plot(legend)
}