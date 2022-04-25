# Script to analyse phylogenetic profiles of orthogroups related to nifH
# By student: B195515
# Course: Comparative and Evolutionary Genomics 2022

setwd("~/CEG_HOME/ICA2/")

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

#----------Obtain Hamming distance
# Modified from https://www.geeksforgeeks.org/how-to-calculate-hamming-distance-in-r/
computeHammingDistance <- function(v1, v2){
  ans <- 0 #initialize value
  for (index in 1:length(v1)){ #iterate over vector length
    if (v1[index] != v2[index]){
      ans <- ans+1} # when not equal, increase distance count
  }
  return (ans)
}

# Process HD for each orthogroup compared to the NifH orthogroup
getOrthogroup <- function(original, modified){
  count <- 0
  for (row in 1:nrow(original)){
    v1 <- original[row,] # vector to compare
    v2 <- original[nifH,] # reference vector
    hd <- computeHammingDistance(v1, v2) # get the HD score
    modified$Hamming_distance[row] <- hd # fill in the columns
    modified$Hd_grep[row] <- paste0("hd",hd,"n") 
    if (hd == 0){
      print(paste(rownames(v1),nifH)) # show output to screen
      count <- count+1
    }
  }
  print(paste(count, "profiles matched with the NifH profile."))
  flename <- paste0("phyProfiles/",modified,".txt")
  write.table(modified, flename, sep="\t")
}

# plot of binary values, colored by orthogroup status
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

# save plot in eps format
saveme <- function(x, plotname){
  ggsave(filename=plotname,
    plot = plotme(x),
    width = 20, units = "cm",
    device = "eps")
}

# output from makephyprofiles.sh
iv1.3 <- read.table("phyProfiles/phyprofiles_iv1.3.txt", row.names=1, header=TRUE)
nifH <- "OG0002449"
iv1.3out <- iv1.3 # make a copy to not interrupt the calculation
iv1.3out[,'Hamming_distance'] <- NA # create empty column for the score
iv1.3out[,'Hd_grep'] <- NA # create empty column for easy grep search
iv1.3out <- getOrthogroup(iv1.3, iv1.3out) # compute HD

col_drop <- c('Hamming_distance', 'Hd_grep') # drop these cols
col_orig <- c("AM","LM","CW","RO","GV","L","MA","NS","TA","S","Annotation","Status") # abbrev species name
col_order <- c("CW","TA","L","RO","NS","AM","GV","LM","MA","S","Annotation","Status") # rearrange by grouping
brewer.pal <- (3, "PuOr") # get color palette

# Annotations have been manually examined/modified in the annot.tsv file before loading in R
# Based on sequence searches in UniProt/PFAM/InterPro/BLAST

hd0Annot <- read.table("hd0annot.tsv", sep="\t", header=TRUE)
hd0 <- dplyr::filter(iv1.3out, grepl('hd0n$', Hd_grep)) %>% select(-one_of(col_drop))
hd0 <- cbind(hd0, hd0Annot$Annotation, hd0Annot$Status)
colnames(hd0) <- col_orig
hd0 <- hd0[, col_order]
saveme(hd0, "hd0_plot.eps")

hd1Annot <- read.table("hd1annot.tsv", sep="\t", header=TRUE)
hd1 <- dplyr::filter(iv1.3out, grepl('hd1n$', Hd_grep)) %>% select(-one_of(col_drop))
hd1 <- cbind(hd1, hd1Annot$Annotation, hd1Annot$Status)
colnames(hd1) <- col_orig
hd1 <- hd1[, col_order]
#hd1$Annotation[1] <- "PEP-CTERM sorting domain-containing protein (13)"
#hd1$Annotation[2] <- "tetratricopeptide repeat protein, CHAT domain*, HP** (11)"
#hd1$Annotation[6] <- "relaxase/mobilization nuclease domain-containing protein, HP** (7)"
saveme(hd1, "hd1_plot.eps")

#save.image("workspace.Rdata")

#sessionInfo()
#R version 4.1.0 (2021-05-18)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Scientific Linux 7.9 (Nitrogen)
#
#Matrix products: default
#BLAS/LAPACK: /usr/lib64/libopenblasp-r0.3.3.so
#
#locale:
# [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C
# [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8
# [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8
# [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C
#[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C
#
#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base
#
#other attached packages:
#[1] RColorBrewer_1.1-2 reshape2_1.4.4     ggplot2_3.3.5      data.table_1.14.2
#[5] tidyr_1.2.0        dplyr_1.0.8
#
#loaded via a namespace (and not attached):
# [1] Rcpp_1.0.8       magrittr_2.0.2   tidyselect_1.1.2 munsell_0.5.0
# [5] colorspace_2.0-3 R6_2.5.1         rlang_1.0.1      fansi_1.0.2
# [9] stringr_1.4.0    plyr_1.8.6       tools_4.1.0      grid_4.1.0
#[13] gtable_0.3.0     utf8_1.2.2       cli_3.2.0        DBI_1.1.2
#[17] withr_2.4.3      ellipsis_0.3.2   assertthat_0.2.1 tibble_3.1.6
#[21] lifecycle_1.0.1  crayon_1.5.0     purrr_0.3.4      vctrs_0.3.8
#[25] glue_1.6.2       stringi_1.7.6    compiler_4.1.0   pillar_1.7.0
#[29] generics_0.1.2   scales_1.1.1     pkgconfig_2.0.3
