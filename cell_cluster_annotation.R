if(!require("Seurat")) install.packages("Seurat")
if(!require("cluster")) install.packages("cluster")
if(!require("topicmodels")) install.packages("topicmodels")
if(!require("Matrix")) install.packages("Matrix")
if(!require("data.table")) install.packages("data.table")

library(cluster);
library(topicmodels);
library(Matrix)
library(Seurat)
library(data.table);


#' Identify representative PFs (putative functions) for cell clusters
#'
#' @param cellTopic Cell-PF(Document-Topic) matrix.
#' @param clusterResult a vector of clustering result
#' @return A cluster named list of representative PFs
#' @export
#' @examples
#' cellTopic <- read_cellTopic("result/PF52");
#' clusterResult <- LDAFindClusters(cellTopic,clusterNumber = 9)
#' represent_topic <- representTopicCluster(cellTopic,clusterResult);
representTopicCluster <- function (cellTopic, clusterResult) 
{
  clusterResult <- as.character(clusterResult);
  clusterName <- unique(clusterResult);
  #every topic value of every cluster
  cluster_topic <- matrix(0,nrow=dim(cellTopic)[1],ncol=length(clusterName),
                                dimnames = list(1:dim(cellTopic)[1],clusterName));
  for(m in clusterName) {
    cluster_pos <- which(clusterResult==m);
    for(n in 1:dim(cellTopic)[1]){
      temp <- cellTopic[n,cluster_pos];
      cluster_topic[n,m] <- sum(as.numeric(temp))/length(cluster_pos);
    }
  }
  
  distinctiveness <- matrix(0,nrow=dim(cellTopic)[1],ncol=length(clusterName),
                            dimnames = list(1:dim(cellTopic)[1],clusterName));
  represent_topic <-list();
  for(m in clusterName) {
    distinctiveness[,m] <- cluster_topic[,m]/apply(cluster_topic,1,sum);
    orderDistinct <- sort(distinctiveness[,m],decreasing = TRUE)
    difference = c();
    for (i in 1:(length(orderDistinct)-1)) {
      difference[names(orderDistinct[i])] = orderDistinct[i]-orderDistinct[i+1];
    }
    sig_topic <- names(difference[1:which(difference==max(difference))])
    if (length(sig_topic)>1){
      x <- cluster_topic[sig_topic,m]
      if (all(x <= 0.1)) {
        represent_topic[[m]] <- sig_topic
      } else {
        represent_topic[[m]] <- names(x[x > 0.1])
      }
    } else {
      represent_topic[[m]] <- sig_topic
    }
  }
  return(represent_topic)
}



#' Identify representative genes for cell clusters
#' 
#' @description representGeneCluster identifies representative genes for each cell clusters by Identifying their representative PFs (putative functions).
#'
#' @param cellTopic Cell-PF(Document-Topic) matrix.
#' @param topicGene PF-Gene(Topic-Term) matrix.
#' @param clusterResult a vector of clustering result
#' @param n.features A single integer specifying how many top representative genes will be extracted.
#' @return A cluster named list of representative genes
#' @export
#' @examples
#' cellTopic <- read_cellTopic("result/PF52");
#' topicGene <- read_topicGene("result/PF52");
#' clusterResult <- LDAFindClusters(cellTopic,clusterNumber = 9);
#' represent_gene <- representGeneCluster(cellTopic,topicGene,clusterResult,n.features = 200);
representGeneCluster <- function (cellTopic, topicGene, clusterResult, n.features) 
{
  represent_topic <- representTopicCluster(cellTopic,clusterResult);
  clusterResult <- as.character(clusterResult);
  clusterName <- unique(clusterResult);
  #every gene value of 18 cluster
  clusterAllgene <- matrix(0,nrow=dim(topicGene)[1],ncol=length(clusterName),
                                dimnames = list(rownames(topicGene),clusterName));
  for(m in clusterName) {
    if (length(represent_topic[[m]]) > 1) {
      clusterAllgene[,m]<- apply(topicGene[,represent_topic[[m]]],1,mean);
    } else {
      clusterAllgene[,m] <- topicGene[,represent_topic[[m]]];
    }
  }
  clusterAllgene <- clusterAllgene[apply(clusterAllgene,1,sum)!=0,];
  
  #distinctiveness
  distinctiveness <- matrix(0,nrow=dim(clusterAllgene)[1],ncol=length(clusterName),
                     dimnames = list(rownames(clusterAllgene),clusterName));
  for(m in clusterName) {
    distinctiveness[,m] <- clusterAllgene[,m]/apply(clusterAllgene,1,sum);
  }
  

  represent_gene <- list();
  for (m in clusterName){
    accum_value <- cumsum(sort(clusterAllgene[,m],decreasing = T))
    informa_gene <- names(which(accum_value <= 0.85))
    if (length(informa_gene) >= n.features){
    } else {
      informa_gene <- names(accum_value[1:n.features])
    }
    orderAllgene <- matrix(0,nrow=length(informa_gene),ncol=3,
                          dimnames = list(informa_gene,c("value","Dgk","sum")));
    orderAllgene[,1] <- as.numeric(factor(clusterAllgene[informa_gene,m]));
    orderAllgene[,2] <- as.numeric(factor(distinctiveness[informa_gene,m]));
    orderAllgene[,3] <- apply(orderAllgene[,1:2],1,sum);
    A <- sort(orderAllgene[,3],decreasing = T)[1:n.features];
    represent_gene[[m]] <- names(A);
  }
  return(represent_gene)
}


#' LDA automatic cell type prediction
#' 
#' @description clusterHGT identifies representative genes for each cell clusters and performs hypergeometric 
#' test against pre-established cell type marker lists.
#' It returns a score of enrichment in the form of -log10 pvalue. 
#' It can notably be used with cell type signatures to predict cell types or with functionnal pathways 
#'
#' @param cellTopic Cell-PF(Document-Topic) matrix.
#' @param topicGene PF-Gene(Topic-Term) matrix.
#' @param clusterResult a vector of clustering result
#' @param pathways pre-established cell type marker lists
#' @param n.features integer of top n representative genes to consider for hypergeometric test
#' @param minSize minimum number of overlapping genes in pathways
#' @param log.trans if TRUE tranform the pvalue matrix with -log10 and convert it to sparse matrix
#' @param p.adjust if TRUE apply Benjamini Hochberg correction to p-value
#' @importFrom stats phyper
#'
#' @return a matrix of benjamini hochberg adjusted pvalue pvalue or a sparse matrix of (-log10) benjamini hochberg adjusted pvalue
#' @export
#'
#' @examples
#' cellTopic <- read_cellTopic("result/PF52");
#' topicGene <- read_topicGene("result/PF52");
#' clusterResult <- LDAFindClusters(cellTopic,clusterNumber = 9)
clusterHGT <- function (cellTopic, topicGene, clusterResult, pathways, n.features = n.features, minSize = 0, log.trans = TRUE, p.adjust = TRUE) 
{
  represent_gene <- representGeneCluster(cellTopic,topicGene,clusterResult,n.features)
  clusterResult <- as.character(clusterResult);
  message("ranking genes")
  features <- rownames(topicGene)
  cells <- colnames(cellTopic)
  genePos <- pbapply::pbsapply(represent_gene, function(x) which(features %in% x))
  
  i <- as.vector(pbapply::pbsapply(clusterResult,function(x) genePos[,x]))
  j <- rep(seq(length(cells)), each = n.features)
  TargetMatrix <- sparseMatrix(i, j, x = 1, dims = c(length(features), length(cells)), dimnames = list(features, cells))
  
  pathways <- lapply(pathways, function(x) x[x %in% features])
  pathways <- pathways[sapply(pathways, function(x) length(x) >= minSize)]
  message("calculating number of success\n")
  PathwayMat <- pbapply::pbsapply(pathways, function(x) which(features %in% x), simplify = F)
  PathwayLen <- unlist(lapply(PathwayMat, length))
  j <- rep(seq(length(PathwayMat)), times = PathwayLen)
  PathwayMatrix <- sparseMatrix(unlist(PathwayMat), j, x = 1, dims = c(length(features), length(PathwayMat)), dimnames = list(features, names(PathwayMat)))
  
  # Hypergeo ----------------------------------------------------------------
  q <- as.data.frame(as.matrix((t(TargetMatrix) %*% PathwayMatrix) - 1))
  m <- sapply(pathways, function(x) sum(x %in% features))
  n <- sapply(m, function(x) length(features) - x)
  k <- n.features
  message("performing hypergeometric test\n")
  A <- pbapply::pbmapply(FUN = function(q, m, n, k) {
    listhyper <- phyper(seq(-1, max(q)), m, n, k, lower.tail = F)[q + 2]
    return(listhyper)
  }, q = q, m = m, n = n, k = k)
  rownames(A) <- rownames(q)
  A <- t(A)
  if (p.adjust) {
    A <- apply(A, 2, function(x) p.adjust(x, "BH"))
  }
  if (log.trans) {
    A <- as.sparse(-log10(A))
  }
  return(A)
}



#' Identify clusters of cells by LDA and Hellinger distance based K-medoids algorithm
#'
#' @param cellTopic Cell-PF(Document-Topic) matrix.
#' @param clusterNumber positive integer specifying the number of clusters.
#' @return a vector of clustering result
#' @export
#' @examples
#' cellTopic <- read_cellTopic("result/PF52");
#' clusterResult <- LDAFindClusters(cellTopic,clusterNumber = 9)
LDAFindClusters <- function (cellTopic, clusterNumber)
{
  cellTopic <- t(cellTopic)
  cellTopic <- as.matrix(cellTopic)
  delta <- distHellinger(cellTopic);
  delta <- as.dist(delta);
  kc <- pam(delta,clusterNumber);
  clusterResult <- kc$cluster;
  names(clusterResult) <- rownames(cellTopic);
  return(clusterResult)
}


#' Read Cell-PF(Document-Topic) matrix.
#'
#' @param outPrefix Path and Prefix of LDA results
#' @return PF-Gene(Topic-Term) matrix(putative functions as rows, cells as columns)
#' @export
#' @examples
#' cellTopic <- read_cellTopic("result/PF52");
read_cellTopic <- function (outPrefix)
{
  cellTopic <- fread(file=paste(outPrefix,".c2k.txt",sep = ""),header= F,sep="\t",data.table = F);
  cellTopic <- t(cellTopic)
  colnames(cellTopic) <- unlist(c(fread(file=paste(outPrefix,".cell.txt",sep = ""),header= F,sep="\n",data.table = F)));
  rownames(cellTopic) <- as.character(1:dim(cellTopic)[1])
  return(cellTopic)
}


#' Read PF-Gene(Topic-Term) matrix.
#'
#' @param outPrefix Path and Prefix of LDA results
#' @return PF-Gene(Topic-Term) matrix(genes as rows, putative functions as columns)
#' @export
#' @examples
#' topicGene <- read_topicGene("result/PF52");
read_topicGene <- function (outPrefix)
{
  topicGene <- fread(file=paste(outPrefix,".k2g.txt",sep = ""),header= F,sep="\t",data.table = F);
  topicGene <- t(topicGene)
  colnames(topicGene) <- as.character(1:dim(topicGene)[2])
  rownames(topicGene) <- unlist(c(fread(file=paste(outPrefix,".gene.txt",sep = ""),header= F,sep="\n",data.table = F)));
  return(topicGene)
}
