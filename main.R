source("cell_cluster_annotation.R")
library(data.table)


# Load pre-established marker lists --------------------------------------------------------------
markerList <- fread("testData/cell_type_marker.csv",header = T)
markerList <- lapply(markerList, function(z){ z[!is.na(z) & z != ""]})


# Read LDA result ----------------------------------------------------------------
cellTopic <- read_cellTopic("testData/PF25");
topicGene <- read_topicGene("testData/PF25");


# Cell clustering based on LDA ----------------------------------------------------------------
clusterResult <- LDAFindClusters(cellTopic, clusterNumber = 6)


# LDA automatic cell type prediction using pre-established marker lists ----------------------------------------------------------------
LDA_HGT <- clusterHGT(cellTopic, topicGene, clusterResult, pathways = markerList, n.features = 200, minSize = 0, log.trans = TRUE, p.adjust = TRUE)
LDA_prediction <- rownames(LDA_HGT)[apply(LDA_HGT, 2, which.max)]
LDA_prediction_signif <- ifelse(2 < apply(LDA_HGT, 2, max), LDA_prediction, "unassigned")


# Identify representative PFs(putative functions) for cell clusters ----------------------------------------------------------------
represent_topic <- representTopicCluster(cellTopic, clusterResult);


# Identify representative genes for cell clusters ----------------------------------------------------------------
represent_gene <- representGeneCluster(cellTopic, topicGene, clusterResult, n.features = 200);

