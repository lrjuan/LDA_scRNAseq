library(cluster);
library(topicmodels);
library(Matrix)
library(Seurat)
library(SeuratObject)


# Get a representative topic for each cluster  ----------------------------------------------------------------
representTopicCluster <- function (cellTopic,clusterResult) 
{
  clusterResult <- as.character(clusterResult);
  clusterName <- unique(clusterResult);
  #####every topic value of every cluster(mean/median)
  cluster_topic <- matrix(0,nrow=dim(cellTopic)[1],ncol=length(clusterName),
                                dimnames = list(1:dim(cellTopic)[1],clusterName));
  for(m in clusterName) {
    cluster_pos <- which(clusterResult==m);
    for(n in 1:dim(cellTopic)[1]){
      temp <- cellTopic[n,cluster_pos];
      
      #cluster_topic[n,m] <- median(as.numeric(temp), na.rm = FALSE);
      cluster_topic[n,m] <- sum(as.numeric(temp))/length(cluster_pos);
    }
  }
  
  ##########Get a representative topic for each cluster
  distinctiveness <- matrix(0,nrow=dim(cellTopic)[1],ncol=length(clusterName),
                            dimnames = list(1:dim(cellTopic)[1],clusterName));
  represent_topic <-list();
  for(m in clusterName) {
    distinctiveness[,m] <- cluster_topic[,m]/apply(cluster_topic,1,sum);
    ###difference
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


# Get 20 representative gene for each cluster  ----------------------------------------------------------------
representGeneCluster <- function (topicGene,represent_topic,clusterResult,n.features) 
{
  clusterResult <- as.character(clusterResult);
  clusterName <- unique(clusterResult);
  #every gene value of 18 cluster(mean)
  clusterAllgene <- matrix(0,nrow=dim(topicGene)[1],ncol=length(clusterName),
                                dimnames = list(rownames(topicGene),clusterName));
  for(m in clusterName) {
    if (length(represent_topic[[m]]) > 1) {
      clusterAllgene[,m]<- apply(topicGene[,represent_topic[[m]]],1,mean);
    } else {
      clusterAllgene[,m] <- topicGene[,represent_topic[[m]]];###
    }
  }
  clusterAllgene <- clusterAllgene[apply(clusterAllgene,1,sum)!=0,];
  
  #distinctiveness
  distinctiveness <- matrix(0,nrow=dim(clusterAllgene)[1],ncol=length(clusterName),
                     dimnames = list(rownames(clusterAllgene),clusterName));
  for(m in clusterName) {
    distinctiveness[,m] <- clusterAllgene[,m]/apply(clusterAllgene,1,sum);
  }
  
  #Get 20 representative gene for each cluster
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



# hypergeometric test  ----------------------------------------------------------------
clusterHGT <- function (cellTopic, topicGene, clusterResult, pathways, n.features = n.features, features = NULL, minSize = 5, log.trans = TRUE, p.adjust = TRUE) 
{
  clusterResult <- as.character(clusterResult);
  represent_topic <- representTopicCluster(cellTopic,clusterResult)
  represent_gene <- representGeneCluster(topicGene,represent_topic,clusterResult,n.features) 
  message("ranking genes")
  features <- rownames(topicGene)
  cells <- colnames(cellTopic)
  genePos <- pbapply::pbsapply(represent_gene, function(x) which(features %in% x))
  i <- as.vector(pbapply::pbsapply(clusterResult,function(x) genePos[,x]))#行基因
  j <- rep(seq(length(cells)), each = n.features)#列 细胞
  #TargetMatrix返回每个细胞排名前200个的基因分别是哪些，标注为1
  TargetMatrix <- sparseMatrix(i, j, x = 1, dims = c(length(features), length(cells)), dimnames = list(features, cells))
  
  pathways <- lapply(pathways, function(x) x[x %in% features])
  pathways <- pathways[sapply(pathways, function(x) length(x) >= minSize)]
  message("calculating number of success\n")
  PathwayMat <- pbapply::pbsapply(pathways, function(x) which(features %in% x), simplify = F)
  PathwayLen <- unlist(lapply(PathwayMat, length))
  j <- rep(seq(length(PathwayMat)), times = PathwayLen)
  #PathwayMatrix返回是每个参考细胞类型的marker基因在1万个基因的那里
  PathwayMatrix <- sparseMatrix(unlist(PathwayMat), j, x = 1, dims = c(length(features), length(PathwayMat)), dimnames = list(features, names(PathwayMat)))
  #intersect(names(which(TargetMatrix[,3]>0)),names(which(PathwayMatrix[,4]>0)))
  
  q <- as.data.frame(as.matrix((t(TargetMatrix) %*% PathwayMatrix) - 1))#某个细胞前200个与某个参考细胞类型marker交集有几个基因-1
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

