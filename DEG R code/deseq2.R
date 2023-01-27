#import data as datafram(df)
##Pull down data, data preprocessing
setwd("/Users/hmkim/data/quant_data/MLL3trim_072720")
setwd("/Users/hmkim/data/quant_data/MLL3_s1_072820")
setwd("/Users/hmkim/data/quant_data/MLL-KO_073120")
setwd("/Users/hmkim/data/quant_data/LJK")
setwd("/Users/hmkim/data/quant_data/MLL3_star_confirm_082320")

##Automate count data.frame 
#Load data directory
src_dir <- c()
src_dir <- c("/Users/hmkim/data/quant_data/MLL-KO_073120")
src_files <- list.files(src_dir)
src_files <- src_files[!src_files %in% "summary"]
src_files
dirVec <- c()
for (i in src_files){
    direc <- paste(getwd(), "/", i, sep = "")
    dirVec[length(dirVec) + 1] <- direc
}
tableList <- list()
for (i in dirVec){
    tableList[[length(tableList) + 1]] <- read.table(i, sep = "\t", header = TRUE)
}
#count dataframe
cnts <- data.frame()
for (i in 1:length(tableList)){
    if (i == 1){
        cnts <- data.frame(c(tableList[[i]][, 7]))
    }else{
        cnts <- cbind(cnts, c(tableList[[i]][, 7]))
    }
}
rownames(cnts) <- c(tableList[[1]]$Geneid)
cnts <- cnts[rowSums(cnts > 1) >=length(cnts),]#drop genes with low counts, this is necessary for rld quality
names(cnts) <- c("MLL1_1", "MLL1_2", "MLL2_1", "MLL2_2", "MLL4_1", "MLL4_2", "WT1", "WT2") #change col name
names(cnts) <- c("MLL1_1", "MLL1_2", "MLL2_1", "MLL2_2", "MLL3_1", "MLL3_2", "WT1", "WT2") #change col name

# MLL1_1rn=read.table("MLL1-KO-RNA1Aligned.out.sam.txt", sep="\t", header=TRUE)#for extract geneid
# 
# MLL1_1df=read.table("MLL1-KO-RNA1Aligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
# MLL1_2df=read.table("MLL1-KO-RNA2Aligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
# MLL2_1df=read.table("MLL2-KO-RNA2Aligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
# MLL2_2df=read.table("MLL2-KO-RNA3Aligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
# MLL4_1df=read.table("MLL4-KO-RNA1Aligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
# MLL4_2df=read.table("MLL4-KO-RNA2Aligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
# WT1df=read.table("WT-HCT116-R1Aligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
# WT2df=read.table("WT-HCT116-R2Aligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
# 
# #extract read count column and integrate all samples counts in cnts
# MLL1_1 <- c(MLL1_1df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL_KO_star_hg38.MLL1.KO.RNA1Aligned.out.sam)
# MLL1_2 <- c(MLL1_2df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL_KO_star_hg38.MLL1.KO.RNA2Aligned.out.sam)
# MLL2_1 <- c(MLL2_1df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL_KO_star_hg38.MLL2.KO.RNA2Aligned.out.sam)
# MLL2_2 <- c(MLL2_2df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL_KO_star_hg38.MLL2.KO.RNA3Aligned.out.sam)
# MLL4_1 <- c(MLL4_1df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL_KO_star_hg38.MLL4.KO.RNA1Aligned.out.sam)
# MLL4_2 <- c(MLL4_2df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL_KO_star_hg38.MLL4.KO.RNA2Aligned.out.sam)
# WT1 <- c(WT1df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL_KO_star_hg38.WT.HCT116.R1Aligned.out.sam)
# WT2 <- c(WT2df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL_KO_star_hg38.WT.HCT116.R2Aligned.out.sam)
# Geneid <- c(MLL1_1rn$Geneid)
# cnts<-data.frame(MLL1_1, MLL1_2, MLL2_1, MLL2_2, MLL4_1, MLL4_2, WT1, WT2, Geneid, row.names = "Geneid")

#make condition table for downstream processing
condition <- c("MLL1", "MLL1", "MLL2", "MLL2", "MLL4", "MLL4", "WT", "WT")
sampleName <- c("MLL1_1", "MLL1_2", "MLL2_1", "MLL2_2", "MLL4_1", "MLL4_2", "WT1", "WT2")
condition <- c("MLL1", "MLL1", "MLL2", "MLL2", "MLL3", "MLL3", "WT", "WT")
sampleName <- c("MLL1_1", "MLL1_2", "MLL2_1", "MLL2_2", "MLL3_1", "MLL3_2", "WT1", "WT2")
sampleTable <- data.frame(condition, sampleName, row.names="sampleName")
sampleTable$condition <- as.factor(sampleTable$condition)
all(rownames(sampleTable) == colnames(cnts)) #check condition table is correct

##Conduct DESeq2
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cnts,
                              colData = sampleTable,
                              design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "WT")#Make WT sample level as reference
dds <- DESeq(dds)
dds$condition #check dds



##Draw MA plots
#MLL1-KO
res_MLL1 <- results(dds, contrast = c("condition", "MLL1", "WT"), alpha = 0.01)
resultsNames(dds)
resLFC_MLL1 <- lfcShrink(dds, coef="condition_MLL1_vs_WT", type="apeglm", res=res_MLL1)
plotMA(resLFC_MLL1, ylim=c(-5, 5))
#resOrdered_lfc_MLL1 <- resLFC_MLL1[order(resLFC_MLL1$padj),]
#resOrdered_MLL1_pval <- res_MLL1[order(res_MLL1$pvalue),]
#resOrdered_MLL1_lfc <- res1[order(res1$log2FoldChange),]
#summary(res_MLL1)
summary(resLFC_MLL1)
resLFC_MLL1 <- subset(resLFC_MLL1, padj<0.01)
summary(resLFC_MLL1)
resLFC_MLL1_up <- subset(resLFC_MLL1, log2FoldChange > 0)
summary(resLFC_MLL1_up)
resLFC_MLL1_down <- subset(resLFC_MLL1, log2FoldChange < 0)
summary(resLFC_MLL1_down)
#write.csv(as.data.frame(resOrdered1p), 
#          file="~/data/sigtest1.csv")
write.csv(as.data.frame(resLFC_MLL1_up), 
          file = "~/data/deg_data/public_data_confirm/resLFC_MLL1_up.csv")
write.csv(as.data.frame(resLFC_MLL1_down), 
          file = "~/data/deg_data/public_data_confirm/resLFC_MLL1_down.csv")
write.csv(as.data.frame(resLFC_MLL1_lfc_down2), 
          file = "~/data/deg_data/resLFC_MLL1_padj01_down2.csv")
MLL1_up_rowname <- c(rownames(resLFC_MLL1_up))
MLL1_down_rowname <- c(rownames(resLFC_MLL1_down))
write.csv(as.data.frame(MLL1_up_rowname), 
          file = "~/data/deg_data/public_data_confirm/MLL1_up_rowname.csv")
write.csv(as.data.frame(MLL1_down_rowname), 
          file = "~/data/deg_data/public_data_confirm/MLL1_down_rowname.csv")

#resSig1 <- subset(res1, padj < 0.01)
#resSig1 <- subset(resSig1, log2FoldChange>1 | log2FoldChange< -1)
#resSig1 <- subset(resSig1, log2FoldChange>1)


#MLL2-KO
res_MLL2 <- results(dds, contrast = c("condition", "MLL2", "WT"), alpha = 0.01)
resultsNames(dds)
resLFC_MLL2 <- lfcShrink(dds, coef="condition_MLL2_vs_WT", type="apeglm", res=res_MLL2)
plotMA(resLFC_MLL2, ylim=c(-5, 5))
#resOrdered_lfc_MLL2 <- resLFC_MLL2[order(resLFC_MLL2$padj),]
#resOrdered_MLL1_pval <- res_MLL1[order(res_MLL1$pvalue),]
#resOrdered_MLL1_lfc <- res1[order(res1$log2FoldChange),]
#summary(res_MLL1)
summary(resLFC_MLL2)
resLFC_MLL2 <- subset(resLFC_MLL2, padj<0.01)
summary(resLFC_MLL2)
resLFC_MLL2_up <- subset(resLFC_MLL2, log2FoldChange > 0)
summary(resLFC_MLL2_up)
resLFC_MLL2_down <- subset(resLFC_MLL2, log2FoldChange < 0)
summary(resLFC_MLL2_down)
#write.csv(as.data.frame(resOrdered1p), 
#          file="~/data/sigtest1.csv")
write.csv(as.data.frame(resLFC_MLL2_up), 
          file = "~/data/deg_data/public_data_confirm/resLFC_MLL2_up.csv")
write.csv(as.data.frame(resLFC_MLL2_down), 
          file = "~/data/deg_data/public_data_confirm/resLFC_MLL2_down.csv")
write.csv(as.data.frame(resLFC_MLL1_lfc_down2), 
          file = "~/data/deg_data/resLFC_MLL2_padj01_down2.csv")
MLL2_up_rowname <- c(rownames(resLFC_MLL2_up))
MLL2_down_rowname <- c(rownames(resLFC_MLL2_down))
write.csv(as.data.frame(MLL2_up_rowname), 
          file = "~/data/deg_data/public_data_confirm/MLL2_up_rowname.csv")
write.csv(as.data.frame(MLL2_down_rowname), 
          file = "~/data/deg_data/public_data_confirm/MLL2_down_rowname.csv")
#resSig1 <- subset(res1, padj < 0.01)
#resSig1 <- subset(resSig1, log2FoldChange>1 | log2FoldChange< -1)
#resSig1 <- subset(resSig1, log2FoldChange>1)


#MLL4-KO
res_MLL3 <- results(dds, contrast = c("condition", "MLL3", "WT"), alpha = 0.01)
resultsNames(dds)
resLFC_MLL3 <- lfcShrink(dds, coef="condition_MLL3_vs_WT", type="apeglm", res=res_MLL4)
plotMA(resLFC_MLL3, ylim=c(-5, 5))
#resOrdered_lfc_MLL4 <- resLFC_MLL4[order(resLFC_MLL4$padj),]
#resOrdered_MLL1_pval <- res_MLL1[order(res_MLL1$pvalue),]
#resOrdered_MLL1_lfc <- res1[order(res1$log2FoldChange),]
#summary(res_MLL1)
summary(resLFC_MLL3)
resLFC_MLL3 <- subset(resLFC_MLL3, padj<0.01)
summary(resLFC_MLL3)
resLFC_MLL3_up <- subset(resLFC_MLL3, log2FoldChange > 0)
summary(resLFC_MLL3_up)
resLFC_MLL3_down <- subset(resLFC_MLL3, log2FoldChange < 0)
summary(resLFC_MLL3_down)
#write.csv(as.data.frame(resOrdered1p), 
#          file="~/data/sigtest1.csv")
write.csv(as.data.frame(resLFC_MLL3_up), 
          file = "~/data/deg_data/public_data_confirm/resLFC_MLL3_up.csv")
write.csv(as.data.frame(resLFC_MLL3_down), 
          file = "~/data/deg_data/public_data_confirm/resLFC_MLL3_down.csv")
write.csv(as.data.frame(resLFC_MLL4_lfc_down2), 
          file = "~/data/deg_data/resLFC_MLL4_padj01_down2.csv")
MLL3_up_rowname <- c(rownames(resLFC_MLL3_up))
MLL3_down_rowname <- c(rownames(resLFC_MLL3_down))
write.csv(as.data.frame(MLL3_up_rowname), 
          file = "~/data/deg_data/public_data_confirm/MLL3_up_rowname.csv")
write.csv(as.data.frame(MLL3_down_rowname), 
          file = "~/data/deg_data/public_data_confirm/MLL3_down_rowname.csv")
#resSig1 <- subset(res1, padj < 0.01)
#resSig1 <- subset(resSig1, log2FoldChange>1 | log2FoldChange< -1)
#resSig1 <- subset(resSig1, log2FoldChange>1)

library(gplots)
data <- list("MLL1_up" = MLL1_up2, "MLL2_up" = MLL2_up2, "MLL4_up" = MLL4_up2)
venn(data)
data <- list("MLL1_down" = MLL1_down2, "MLL2_down" = MLL2_down2, "MLL4_down" = MLL4_down2)
venn(data)

MLL1_KO_diff_up2 <- setdiff(MLL1_up2, MLL2_up2)
MLL1_KO_diff_up2 <- setdiff(MLL1_KO_diff_up2, MLL4_up2)
MLL1_KO_diff_down2 <- setdiff(MLL1_down2, MLL2_down2)
MLL1_KO_diff_down2 <- setdiff(MLL1_KO_diff_down2, MLL4_down2)
MLL2_KO_diff_up2 <- setdiff(MLL2_up2, MLL1_up2)
MLL2_KO_diff_up2 <- setdiff(MLL2_KO_diff_up2, MLL4_up2)
MLL2_KO_diff_down2 <- setdiff(MLL2_down2, MLL1_down2)
MLL2_KO_diff_down2 <- setdiff(MLL2_KO_diff_down2, MLL4_down2)
MLL4_KO_diff_up2 <- setdiff(MLL4_up2, MLL1_up2)
MLL4_KO_diff_up2 <- setdiff(MLL4_KO_diff_up2, MLL2_up2)
MLL4_KO_diff_down2 <- setdiff(MLL4_down2, MLL1_down2)
MLL4_KO_diff_down2 <- setdiff(MLL4_KO_diff_down2, MLL2_down2)

MLL_KO_intersect_down <- intersect(MLL1_down_rowname, MLL2_down_rowname)
MLL_KO_intersect_down <- intersect(MLL_KO_intersect_down, MLL3_down_rowname)
MLL_KO_intersect_up <- intersect(MLL1_up_rowname, MLL2_up_rowname)
MLL_KO_intersect_up <- intersect(MLL_KO_intersect_up, MLL3_up_rowname)


write.csv(as.data.frame(MLL1_KO_diff_up2), 
          file = "~/data/deg_data/MLL1_KO_diff_up2.csv")
write.csv(as.data.frame(MLL1_KO_diff_down2), 
          file = "~/data/deg_data/MLL1_KO_diff_down2.csv")
write.csv(as.data.frame(MLL2_KO_diff_up2), 
          file = "~/data/deg_data/MLL2_KO_diff_up2.csv")
write.csv(as.data.frame(MLL2_KO_diff_down2), 
          file = "~/data/deg_data/MLL2_KO_diff_down2.csv")
write.csv(as.data.frame(MLL4_KO_diff_up2), 
          file = "~/data/deg_data/MLL4_KO_diff_up2.csv")
write.csv(as.data.frame(MLL4_KO_diff_down2), 
          file = "~/data/deg_data/MLL4_KO_diff_down2.csv")
write.csv(as.data.frame(MLL_KO_intersect_down2), 
          file = "~/data/deg_data/MLL_KO_intersect_down2.csv")

resLFC_MLL1df <- read.csv("~/data/deg_data/resLFC_MLL1_padj01.csv", header = TRUE, sep = ",")
names(resLFC_MLL1df)[1] = c("Geneid")
resLFC_MLL2df <- read.csv("~/data/deg_data/resLFC_MLL2_padj01.csv", header = TRUE, sep = ",")
names(resLFC_MLL2df)[1] = c("Geneid")
resLFC_MLL4df <- read.csv("~/data/deg_data/resLFC_MLL4_padj01.csv", header = TRUE, sep = ",")
names(resLFC_MLL4df)[1] = c("Geneid")

################################################
##venn diagram
up_public <- read.csv("~/data/deg_data/public_data_confirm/up_public.csv", header = F, sep = ",")
up_public <- c(up_public$V1)
down_public <- read.csv("~/data/deg_data/public_data_confirm/down_public.csv", header = F, sep = ",")
down_public <- c(down_public$V1)
up_inter_sym <- read.csv("~/data/deg_data/public_data_confirm/MLL_KO_intersect_up_genesym.csv", header = F, sep = ",")
up_inter_sym <- c(up_inter_sym$V1)
down_inter_sym <- read.csv("~/data/deg_data/public_data_confirm/MLL_KO_intersect_down_genesym.csv", header = F, sep = ",")
down_inter_sym <- c(down_inter_sym$V1)
library(gplots)
data <- list("Public_up" = up_public, "Pipeline_up" = up_inter_sym)
venn(data)
data <- list("Public_down" = down_public, "Pipeline_down" = down_inter_sym)
venn(data)
data <- list("Public_up" = up_public, "Pipeline_down" = down_inter_sym)
venn(data)
data <- list("Public_down" = down_public, "Pipeline_up" = up_inter_sym)
venn(data)
##venn diagram end
##################################################


MLLintersect <- intersect(c(resLFC_MLL1df$Geneid), c(resLFC_MLL2df$Geneid))
MLLintersect <- intersect(MLLintersect, c(resLFC_MLL4df$Geneid))
resLFC_MLLKO_padj01_intersect <- subset(resLFC_MLL1df, Geneid %in% MLLintersect)
write.csv(as.data.frame(MLL_KO_intersect_down), 
          file = "~/data/deg_data/public_data_confirm/MLL_KO_intersect_down.csv")
write.csv(as.data.frame(MLL_KO_intersect_up), 
          file = "~/data/deg_data/public_data_confirm/MLL_KO_intersect_up.csv")

rld <- rlog(resLFC_MLLKO_padj01_intersect, blind = FALSE)

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
scaledata <- t(scale(t(rldTable)))
scaledata <- scaledata[complete.cases(scaledata),]
#clustering row(gene) and column(sample) with correlation analysis
hr <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete") # Clusters columns by Spearman correlation.
#pdf("~/data/asdfasdfasdftest.pdf")#For export. Not necessary if you are using Rstudio
heatmap.2(rldTable,
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hc),
          scale="row",
          trace="none",
          dendrogram="col",
          cexRow = 0.4,
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255), 
          labRow = NA
)
#dev.off()#For export. Not necessary if you are using Rstudio
write.csv(assay(rld)[topVarGenes,],file="~/data/testMLL3trim.csv")
write.csv(topVarGenes,file="~/data/testMLL3trim.csv")
write.csv(assay(rld),file="~/data/testMLL3trim.csv")
geneex <- assay(rld)[topVarGenes, ]
library(venn)
library(gplots)

prot1 <- read.csv("~/data/prot1.csv", header = FALSE)
prot2 <- read.csv("~/data/prot2.csv", header = FALSE)
data <- list(M = prot1, P = prot2)

venn(data)



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
