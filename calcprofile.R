# Script to analyse phylogenetic profiles of orthogroups
# related to nifH
# By student: B195515
# Course: Comparative and Evolutionary Genomics 2022

setwd("~/CEG_HOME/ICA2/phyProfiles")

iv1.3 <- read.table("phyprofiles_iv1.3.txt", row.names=1, header=TRUE)
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

iv1.3out <- iv1.3
iv1.3out[,'Hamming_distance'] <- NA # create empty column for the score
iv1.3out[,'Hd_grep'] <- NA

count <- 0
for (row in 1:nrow(iv1.3)){
  v1 <- iv1.3[row,]
  v2 <- iv1.3[nifH,]
  hd <- computeHammingDistance(v1, v2)
  iv1.3out$Hamming_distance[row] <- hd
  iv1.3out$Hd_grep[row] <- paste0("hd",hd) 
  if (hd == 0){
    print(paste(rownames(v1),nifH))
    count <- count+1
  }
}
print(paste(count, "profiles matched with the nifH profile."))

# save the output
write.table(iv1.3out, "iv1.3.out.txt", sep="\t")
