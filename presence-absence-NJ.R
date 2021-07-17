#this script was used to produce figure S7


library (ape)
library(phangorn)

#set your working directory on the roary output
setwd("ICEscore_40")

#upload your data
mat <- read.csv("gene_presence_absence.Rtab", sep = "\t", header = TRUE)

#remove the first column, you do not need it
mat <- mat[, -c(1)]

#transpose the matrix
tmat <- t(mat)

#calculate distance matrix and build a Neighbor-Joining tree 
tmat.nj = nj(dist.gene(tmat, method = "pairwise"))

#optional, you can root the tree a midpoint
rttree <-midpoint(tmat.nj)

#check your tree
plot(rttree, "phylo",  x.lim = 1100, y.lim = 55)

#all good? save it
write.tree(rttree, 'ROARY_presence-absence_NJ.newick')



# This was not used in the paper, however I wanted to check I had similar outputs so I did this, and here it is

#instead of using the output of roary, you can construct the same tree using as input the "proteinortho.tsv" file produced by proteinortho
# FYI, the ICEs were clustered similarly (the tree separated ICESyms from other ICEs and plasmids, and also separated ICESyms based on their host-strain legume association).

#set your working directory on the proteinortho output 
setwd("proteinortho")

mat <- read.csv("UniqueICEs_prortho.proteinortho.tsv", sep = "\t", header = TRUE)

#remove the columns you do not need
mat <- mat[, -c(1:3)]

#we are now creating a binary matrix, where gene absence is 0 and presence is 1
mat[ mat == "*" ]<- 0
mat[mat != 0] <- 1

#transpose matrix, calculate distance and build tree
tmat <- t(mat)
tmat.nj = nj(dist.gene(tmat, method = "pairwise"))

rttree <-midpoint(tmat.nj)


plot(rttree, "phylo",  x.lim = 1100, y.lim = 55)

write.tree(rttree, 'PRORTHO_presence-absence_NJ.newick')



