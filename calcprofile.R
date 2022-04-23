# Script to analyse phylogenetic profiles of orthogroups
# related to nifH
# By student: B195515
# Course: Comparative and Evolutionary Genomics 2022

setwd("~/CEG_HOME/ICA2/")

iv1.3 <- read.table("phyProfiles/phyprofiles_iv1.3.txt", row.names=1, header=TRUE)
nifH <- "OG0002449"

#----------Obtain Hamming distance

# Modified from https://www.geeksforgeeks.org/how-to-calculate-hamming-distance-in-r/
computeHammingDistance <- function(v1, v2){
  ans <- 0 #initialize value
  for (index in 1:length(v1)){ #iterate over vector length
    if (v1[index] != v2[index]){
      ans <- ans+1 # when not equal, increase distance count
    }
  }
  return (ans)
}

iv1.3out <- iv1.3 # make a copy to not interrupt the calculation
iv1.3out[,'Hamming_distance'] <- NA # create empty column for the score
iv1.3out[,'Hd_grep'] <- NA # create empty column for easy grep search

count <- 0
# for each row ie. orthogroup
for (row in 1:nrow(iv1.3)){
  v1 <- iv1.3[row,] # vector to compare
  v2 <- iv1.3[nifH,] # reference vector
  hd <- computeHammingDistance(v1, v2) # get the HD score
  iv1.3out$Hamming_distance[row] <- hd # fill in the columns
  iv1.3out$Hd_grep[row] <- paste0("hd0",hd) 
  if (hd == 0){
    print(paste(rownames(v1),nifH)) # show output to screen
    count <- count+1
  }
}
print(paste(count, "profiles matched with the NifH profile."))

write.table(iv1.3out, "phyProfiles/iv1.3.out.txt", sep="\t")

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

col_drop <- c('Hamming_distance', 'Hd_grep') # drop these cols
col_orig <- c("AM","LM","CW","RO","GV","L","MA","NS","TA","S","Annotation","Status") # abbrev species name
col_order <- c("CW","TA","L","RO","NS","AM","GV","LM","MA","S","Annotation","Status") # rearrange by grouping
brewer.pal <- (3, "PuOr") # get color palette

plotme <- function(x){
ggplot(melt(x),
aes(variable, Annotation, fill=Status, alpha=value)) + 
geom_tile(colour="gray50") + 
scale_alpha_identity(guide="none") +
scale_fill_manual(values=c("Known"="#dcdee3","Novel"="#F1A340","Hypothetical"="#998EC3")) +
coord_equal(expand=0) + 
theme_bw() + 
theme(panel.grid.major=element_blank(),
axis.text.x=element_text(angle=45, hjust=1),
text=element_text(size = 11)) + 
labs(x="Species", y="Annotation / Group size")
}

saveme <- function(x, plotname){
ggsave(filename=plotname,
  plot= plotme(x),
  width = 20, units = "cm",
  device = "eps")
  }

hd0Annot <- read.table("hd0annot.tsv", sep="\t", header=TRUE)
hd0 <- dplyr::filter(iv1.3out, grepl('hd00$', Hd_grep)) %>% select(-one_of(col_drop))
hd0 <- cbind(hd0, hd0Annot$Annotation, hd0Annot$Status)
colnames(hd0) <- col_orig
hd0 <- hd0[, col_order]
saveme(hd0, "hd0_plot.eps")

hd1Annot <- read.table("hd1annot.tsv", sep="\t", header=TRUE)
hd1 <- dplyr::filter(iv1.3out, grepl('hd01$', Hd_grep)) %>% select(-one_of(col_drop))
hd1 <- cbind(hd1, hd1Annot$Annotation, hd1Annot$Status)
colnames(hd1) <- col_orig
hd1 <- hd1[, col_order]
#hd1$Annotation[1] <- "PEP-CTERM sorting domain-containing protein (13)"
#hd1$Annotation[2] <- "tetratricopeptide repeat protein, CHAT domain*, HP** (11)"
#hd1$Annotation[6] <- "relaxase/mobilization nuclease domain-containing protein, HP** (7)"
saveme(hd1, "hd1_plot.eps")

save.image("workspace.Rdata")
