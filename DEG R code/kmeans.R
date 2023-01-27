###K-means clusterting
#filtering with mean and variance
krld <- as.matrix(assay(rld))
krld_var <- apply(krld, 1, var)
krld_mean <- apply(krld, 1, mean)
plot(log2(krld_mean), log2(krld_var), pch='.')
abline(h=log2(0.08), col='red')
abline(v=log2(8.5), col='red')
krld <- krld[which(krld_var > 0.16 & krld_mean > 8), 1:8]
#scale data
scaledata.k <- t(scale(t(krld)))
##choose K
#SSE elbow
wss <- (nrow(scaledata.k)-1)*sum(apply(scaledata.k,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(scaledata.k,
                                     centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#average silhouette width
library(cluster)
sil <- rep(0, 20)
#repeat k-means for 1:20 and extract silhouette:
for(i in 2:20){
  k1to20 <- kmeans(scaledata.k, centers = i, nstart = 25, iter.max = 20)
  ss <- silhouette(k1to20$cluster, dist(scaledata.k))
  sil[i] <- mean(ss[, 3])
}
# Plot the  average silhouette width
plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
abline(v = which.max(sil), lty = 2)
cat("Average silhouette width optimal number of clusters:", which.max(sil), "\n")

#Calinsky criterion
library(vegan)
fit <- cascadeKM(scaledata.k, 1, 20, iter = 100)
plot(fit, sortg = TRUE, grpmts.plot = TRUE)
calinski.best <- as.numeric(which.max(fit$results[2,]))
cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")

##clustering with k=4
set.seed(20)
kClust <- kmeans(scaledata.k, centers=4, nstart = 1000, iter.max = 20)
kClusters <- kClust$cluster
# function to find centroid in cluster i
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}
kClustcentroids <- sapply(levels(factor(kClusters)), clust.centroid, scaledata.k, kClusters)
library(ggplot2)
library(reshape)
#get in long form for plotting
Kmolten <- melt(kClustcentroids)
colnames(Kmolten) <- c('sample','cluster','value')
#plot
p1 <- ggplot(Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_point() + 
  geom_line() +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster Expression by Time",color = "Cluster")
p1
#a posteriori
cor(kClustcentroids)

##clustering with k=3
set.seed(20)
kClust <- kmeans(scaledata.k, centers=3, nstart = 1000, iter.max = 20)
kClusters <- kClust$cluster
# function to find centroid in cluster i
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}
kClustcentroids <- sapply(levels(factor(kClusters)), clust.centroid, scaledata.k, kClusters)
library(ggplot2)
library(reshape)
#get in long form for plotting
Kmolten <- melt(kClustcentroids)
colnames(Kmolten) <- c('sample','cluster','value')
#plot
p1 <- ggplot(Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_point() + 
  geom_line() +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster Expression by Time",color = "Cluster")
p1
#a posteriori
cor(kClustcentroids)

##clustering with k=2
set.seed(20)
kClust <- kmeans(scaledata.k, centers=2, nstart = 1000, iter.max = 20)
kClusters <- kClust$cluster
# function to find centroid in cluster i
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}
kClustcentroids <- sapply(levels(factor(kClusters)), clust.centroid, scaledata.k, kClusters)
library(ggplot2)
library(reshape)
#get in long form for plotting
Kmolten <- melt(kClustcentroids)
colnames(Kmolten) <- c('sample','cluster','value')
#plot
p1 <- ggplot(Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_point() + 
  geom_line() +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster Expression by Time",color = "Cluster")
p1
#a posteriori
cor(kClustcentroids)

##finding core gene
#Subset the cores molten dataframe so we can plot the core
core2 <- Kmolten[Kmolten$cluster=="1",]
#get cluster 2
K2 <- (scaledata.k[kClusters==1,])
#calculate the correlation with the core
corscore <- function(x){cor(x,core2$value)}
score <- apply(K2, 1, corscore)
#get the data frame into long format for plotting
K2molten <- melt(K2)
colnames(K2molten) <- c('gene','sample','value')
#add the score
K2molten <- merge(K2molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K2molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K2molten$order_factor <- 1:length(K2molten$gene)
#order the dataframe by score
K2molten <- K2molten[order(K2molten$score),]
#set the order by setting the factors
K2molten$order_factor <- factor(K2molten$order_factor , levels = K2molten$order_factor)
# Everything on the same plot
p2 <- ggplot(K2molten, aes(x=sample,y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  #this adds the core 
  geom_line(data=core2, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster 2 Expression by Time",color = "Score")
p2

##Comparing cluster methode
#we made the hr and TreeR objects above.
hr <- hclust(as.dist(1-cor(t(scaledata.k), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hclustk4 = cutree(hr, k=4) #cut tree to find 4 clusters
library(dendextend)
TreeR = as.dendrogram(hr, method="complete")
plot(TreeR,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")
#this plots the bar below:
colored_bars(hclustk4, TreeR, sort_by_labels_order = T, y_shift=-0.1, rowLabels = c("k=4"),cex.rowLabels=0.7)
#line up with K-means and hierarchical
plot(TreeR,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")
the_bars <- cbind(hclustk4, kClusters)
colored_bars(the_bars, TreeR, sort_by_labels_order = T, y_shift=-0.1, rowLabels = c("Treecut",'K-means'),cex.rowLabels=0.7)

##comparing k-means and hierarchical
#these functions from the WCGNA package are great for this:
source('https://raw.githubusercontent.com/cran/WGCNA/master/R/matchLabels.R')
source('https://raw.githubusercontent.com/cran/WGCNA/master/R/accuracyMeasures.R')
hclustk4 <- paste0('H-',hclustk4)
kClusters <- paste0('K-',kClusters)
OT<- overlapTable(hclustk4, kClusters)
#get rid of 0 values...
OT$pTable[OT$pTable == 0] <- 2e-300
textMatrix= paste(signif(OT$countTable, 2), "\n(",
                  signif(OT$pTable, 1), ")", sep= "")
dim(textMatrix)= dim(OT$countTable)
par(mar=c(10,10,10,10))
library(gplots)
heatmap.2(x= -log(OT$pTable),
          dendrogram = "none",
          Colv =F,
          Rowv = F,
          scale = c("none"),
          col="heat.colors",
          na.rm=TRUE,
          cellnote = textMatrix,
          notecol="grey30",
          notecex=0.6,
          trace=c("none"),
          cexRow = 0.8,
          cexCol = 0.8,
          main = "Cluster-Cluster Overlap",
          xlab = "K-means (k=4)",
          ylab = "TreeCut (k=4)"
)

##K-means clustering heatmap
set.seed(5)
km <- kmeans(Kmatrix, 5)
m.kmeans <- cbind(scaledata.k, kClust$cluster)
dim(m.kmeans)
o<-order(m.kmeans[,9])
m.kmeans<-m.kmeans[o,]
write.csv(m.kmeans, file="C:/Users/BM/OneDrive/BM/Project/NGS/data/pkh201911/HTseq_count/rename/pkh201911_kmeans.csv")
heatmap.2(m.kmeans[,1:8], cexCol = 0.8, scale="row",trace="none",dendrogram="column",col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
pheatmap(m.kmeans[,1:8], cluster_rows=F,cluster_cols=F, col=brewer.pal(4,"YlOrRd"),border_color=NA, labels_row = NA)
