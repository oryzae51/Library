##Pull down data, data preprocessing
#import data as datafram(df)
setwd("/Users/hmkim/data/quant_data/mm9/mm9_072120")
MLL1_rn=read.table("SRR1030491.txt", sep="\t", header=TRUE) #for extranct geneid
MLL1_1df=read.table("SRR1030491.txt", sep="\t", header=TRUE, row.names="Geneid")
MLL1_2df=read.table("SRR1030500.txt", sep="\t", header=TRUE, row.names="Geneid")
MLL1_3df=read.table("SRR1030501.txt", sep="\t", header=TRUE, row.names="Geneid")
MLL1_4df=read.table("SRR1030502.txt", sep="\t", header=TRUE, row.names="Geneid")

#extract read count column and integrate all samples counts in cnts
nc1 <- c(MLL1_1df$X.media.bm.ETL4TiB.KHM.align_data.SAMfile.MLL1test_M10.SRR1030491)
mm1 <- c(MLL1_2df$X.media.bm.ETL4TiB.KHM.align_data.SAMfile.MLL1test_M10.SRR1030500)
nc2 <- c(MLL1_3df$X.media.bm.ETL4TiB.KHM.align_data.SAMfile.MLL1test.temp.SRR1030501)
mm2 <- c(MLL1_4df$X.media.bm.ETL4TiB.KHM.align_data.SAMfile.MLL1test.temp.SRR1030502)
Geneid <- c(MLL1_rn$Geneid)
cnts<-data.frame(nc1, mm1, nc2, mm2, Geneid, row.names = "Geneid")
cnts <- cnts[rowSums(cnts > 1) >=4,]#drop genes with low counts, this is necessary for rld quality
#make condition table for downstream processing
condition <- c("nc", "mm", "nc", "mm")
sampleName <- c("nc1", "mm1", "nc2", "mm2")
sampleTable <- data.frame(condition, sampleName, row.names="sampleName")
sampleTable$condition <- as.factor(sampleTable$condition)
all(rownames(sampleTable) == colnames(cnts)) #check condition table is correct

##Conduct DESeq2
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cnts,
                              colData = sampleTable,
                              design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "nc")#Make WT sample level as reference
dds <- DESeq(dds)
dds$condition #check dds

##Draw MA plots
#nc vs mm
res1 <- results(dds, contrast = c("condition", "nc", "mm"))
resultsNames(dds)
resLFC1 <- lfcShrink(dds, coef="condition_mm_vs_nc", type="apeglm")
plotMA(resLFC1, ylim=c(-5, 5))
resOrdered <- res1[order(res1$pvalue),]
summary(res1)
#MLL1 vs WT
res1 <- results(dds, contrast = c("condition", "MLL1", "WT"))
resultsNames(dds)
resLFC1 <- lfcShrink(dds, coef="condition_MLL1_vs_WT", type="apeglm")
plotMA(resLFC1, ylim=c(-10, 10))
resOrdered <- res1[order(res1$pvalue),]
summary(res1)
#MLL2 vs WT
res2 <- results(dds, contrast = c("condition", "MLL2", "WT"))
resultsNames(dds)
resLFC2 <- lfcShrink(dds, coef="condition_MLL2_vs_WT", type="apeglm")
plotMA(resLFC2, ylim=c(-10, 10))
resOrdered <- res2[order(res2$pvalue),]
summary(res2)
#MLL4 vs WT
res4 <- results(dds, contrast = c("condition", "MLL4", "WT"))
resultsNames(dds)
resLFC4 <- lfcShrink(dds, coef="condition_MLL4_vs_WT", type="apeglm")
plotMA(resLFC4, ylim=c(-10, 10))
resOrdered <- res4[order(res4$pvalue),]
summary(res4)

##Count data transformation/normalization for downstream analyze
#Normalizae with rlog in DESeq2
library("vsn")
rld <- rlog(dds, blind=FALSE)
meanSdPlot(assay(rld))
#Compare with normTrnaform. You don't need to do when you are making data. 
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))

##Data quality assessment by sample clustering and visualization
#It is kind of QC
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])
pheatmap(assay(ntd)[select,])
pheatmap(assay(rld)[select,])

##PCA
sampleDists1 <- dist(t(assay(rld)))
sampleDists <- dist(scale(t(assay(rld))))#recommand this for consistancy
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, sep="-")
#colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA(rld, intgroup=c("condition"))


##hirachy heatmap
library("gplots")
library("RColorBrewer")
library("genefilter")
#get deviation(rowVars()) from assay(rld), and get high deviation genes
# topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), 50)
# topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), 500)
topVarGenes <- head(assay(rld), decreasing=TRUE, 1000)
#scaling dataset
scaledata <- t(scale(t(topVarGenes))) # Centers and scales data.
scaledata <- scaledata[complete.cases(scaledata),]
#clustering row(gene) and column(sample) with correlation analysis
hr <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete") # Clusters columns by Spearman correlation.
#pdf("~/data/asdfasdfasdftest.pdf")#For export. Not necessary if you are using Rstudio
heatmap.2(topVarGenes,
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hc),
          scale="row",
          trace="none",
          dendrogram="both",
          cexRow = 0.4,
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255)
)
#dev.off()#For export. Not necessary if you are using Rstudio
write.csv(assay(rld)[topVarGenes,],file="~/data/testasdf.csv")
write.csv(topVarGenes,file="~/data/test.csv")
write.csv(assay(rld),file="~/data/test.csv")
geneex <- assay(rld)[topVarGenes, ]
#Cluster samples to check for outlier
TreeC = as.dendrogram(hc, methods="average")
plot(TreeC, 
     main = "Sample Clustering", 
     ylab = "Height")
TreeR = as.dendrogram(hr, method="average")
plot(TreeR,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")
#extract cluster and show as a dendragram
hclusth1.5 = cutree(hr, h=1.5) #cut tree at height of 1.5
hclusth1.0 = cutree(hr, h=1.0) #cut tree at height of 1.0
hclusth0.5 = cutree(hr, h=0.5) #cut tree at height of 0.5
library(dendextend)
#plot the tree
plot(TreeR,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")
#add the three cluster vectors
the_bars <- cbind(hclusth0.5, hclusth1.0, hclusth1.5)
#this makes the bar
colored_bars(the_bars, TreeR, sort_by_labels_order = T, y_shift=-0.1, rowLabels = c("h=0.5","h=1.0","h=1.5"),cex.rowLabels=0.7)
#this will add lines showing the cut heights
abline(h=1.5, lty = 2, col="grey")
abline(h=1.0, lty = 2, col="grey")
abline(h=0.5, lty = 2, col="grey")
#Designate the numbers of cluster
hclustk4 = cutree(hr, k=4)
plot(TreeR,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")
colored_bars(hclustk4, TreeR, sort_by_labels_order = T, y_shift=-0.1, rowLabels = c("k=4"),cex.rowLabels=0.7)
#Dynamic cut tree method
library(dynamicTreeCut)
clusDyn <- cutreeDynamic(hr, distM = as.matrix(as.dist(1-cor(t(topVarGenes)))))
plot(TreeR,
    leaflab = "none",
    main = "Gene Clustering",
    ylab = "Height")
colored_bars(clusDyn, TreeR, sort_by_labels_order = T, y_shift=-0.1, rowLabels = c("Dynamic"),cex.rowLabels=0.7)

##k-means clustsering
# topVarGenes<-order(-rowVars(assay(rld)))[0:1000]
# Kmatrix <- assay(rld)[topVarGenes, ]
# Kmatrix <- Kmatrix - rowMeans(Kmatrix)
# set.seed(5)
# km <- kmeans(Kmatrix, 5)
# m.kmeans <- cbind(Kmatrix, km$cluster)
# dim(m.kmeans)
# o<-order(m.kmeans[,9])
# m.kmeans<-m.kmeans[o,]
# write.csv(m.kmeans, file="~/data/kmean-test.csv")
# pdf("C:/Users/BM/OneDrive/BM/Project/NGS/data/pkh201911/HTseq_count/rename/pkh201911_top1000_Kmeans.pdf")
# heatmap.2(m.kmeans[,1:8], cexCol = 0.8, scale="row",trace="none",dendrogram="column",col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
# dev.off()