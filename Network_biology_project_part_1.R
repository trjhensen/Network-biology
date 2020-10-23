#Set working directory
getwd()
setwd("C:/Users/Tim/Documents")

#load packages
install.packages("readr")
install.packages("tidyverse", INSTALL_opts="--no-multiarch")
library(tidyverse)
library(SNFtool)
library(readr)
#Clear global workspace
rm(list = ls())

#Import metatranscriptomics dataset
METATRANS <- read_csv('pathabundances_metatrans.csv') 
#Set samples as rows and bacterial clades as columns
METATRANS <- t(METATRANS)
colnames(METATRANS) <- as.matrix(METATRANS[1, ])
METATRANS <- METATRANS[-1,];
METATRANS <- as.data.frame(METATRANS)

#Import metagenomics dataset
METAGEN <- read_tsv('metagen_path.tsv') 
#Set samples as rows and bacterial clades as columns
METAGEN <- t(METAGEN)
colnames(METAGEN) <- as.matrix(METAGEN[1, ])
METAGEN <- METAGEN[-1,];
METAGEN <- as.data.frame(METAGEN)

#Make sure the same samples are present in both datasets by taking the intersection
library(dplyr)

#Determine the overlapping sample names
METAGEN1 <- rownames(METAGEN)
METATRANS1 <- rownames(METATRANS)
overlap <- intersect(METAGEN1, METATRANS1)
length(overlap)

#Selecting the overlapping sample names in both datasets
METAGENT <- t(METAGEN)
METATRANST <- t(METATRANS)
METAGENov <- t(METAGENT[ , overlap])
METATRANSov <- t(METAGENT[ , overlap])


#From here on, the code is adapted from code provided by 
#Martina Kutmon, Maastricht University.

## First, set all the parameters:
K = 20;		# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
T = 50; 	# Number of Iterations, usually (10~20)

#Calculate the distances using the euclidian distance metric
DGEN <- dist(METAGENov, method = "euclidean", diag = TRUE, upper = TRUE);
DGEN <- as.matrix(DGEN)
DTRANS <- dist(METATRANSov, method = "euclidean", diag = TRUE, upper = TRUE);
DTRANS <- as.matrix(DTRANS)

## next, construct similarity graphs
W1 = affinityMatrix(DGEN, K, alpha)
W2 = affinityMatrix(DTRANS, K, alpha)

truelabel = c(matrix(1,41,1),matrix(2,41,1)); ##the ground truth of the simulated data;
truelabel1 <- as.matrix(truelabel)

## These similarity graphs have complementary information about clusters.
displayClusters(W1,truelabel);
displayClusters(W2,truelabel);

## next, we fuse all the graphs
W = SNF(list(W1,W2), K, T)

## With this unified graph W of size n x n, you can do either spectral clustering or Kernel NMF. 
## for example, spectral clustering

#See which amount of clusters results in the highest NMI value
NMI <- vector()
for (i in 2:82){
  C = i 
  group = spectralClustering(W, C); 	# the final subtypes information
  group1 <- as.data.frame(group)
  SNF_NMI = calNMI(group, truelabel)
  NMI[i] <- SNF_NMI
}
# Simple  Plot
NMI1 <- as.matrix(NMI)
NMI1 <- NMI1[-1,];
plot(NMI1, main="NMI score per number of clusters",
     xlab = "number of clusters",
     ylab = 'NMI score') 

#Lets start with 20 clusters
group2 = spectralClustering(W, 11); # the final subtypes information
group3 <- as.data.frame(group2)

Wcol <- colnames(W)
group3$name <- Wcol
group3 <- group3 %>% arrange(group2)
loadTableData(group3, data.key.column = "name")
write.csv(group3,"11clusters.csv", row.names = TRUE)

## visualization in Cytoscape (abritory cut-off required to choose which edges should be included)
W_cyto <- group2
diag(W_cyto) <- 0

quantile(W_cyto)
W_cyto[W_cyto < 0.01] <- 0
library(igraph)
graph <- graph_from_adjacency_matrix(W_cyto, weighted=TRUE, mode = "undirected")
createNetworkFromIgraph(graph)







