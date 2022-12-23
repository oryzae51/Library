#import data as datafram(df)
##Pull down data, data preprocessing
cnt_dir <- "/Users/kimhyoungmin/data/060922_PKH_ZWC/quant/cnt"
results_dir <- "/Users/kimhyoungmin/Library/Mobile\ Documents/com~apple~CloudDocs/ETL/Experiment data/Cowork/060922_PKH_ZWC/"
setwd("/Users/kimhyoungmin/data/060922_PKH_ZWC/quant/cnt")

##Automate count data.frame 
#Load data directory
src_dir <- c()
src_dir <- c("/Users/kimhyoungmin/data/060922_PKH_ZWC/quant/cnt")
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
names(cnts) <- c("F_ZC3H4_0D_1", "F_ZC3H4_0D_2", "F_ZC3H4_2D_1", "F_ZC3H4_2D_2", "F_ZC3H4_3D_1", "F_ZC3H4_3D_2", "Vector_0D_1", "Vector_0D_2", "Vector_2D_1", "Vector_2D_2", "Vector_3D_1", "Vector_3D_2") #change col name

#make condition table for downstream processing
#condition <- c("F_ZC3H4_0D", "F_ZC3H4_0D", "F_ZC3H4_2D", "F_ZC3H4_2D", "F_ZC3H4_3D", "F_ZC3H4_3D", "Vector_0D", "Vector_0D", "Vector_2D", "Vector_2D", "Vector_3D", "Vector_3D")
treatment <- c("F_ZC3H4", "F_ZC3H4", "F_ZC3H4", "F_ZC3H4", "F_ZC3H4", "F_ZC3H4", "vector", "vector", "vector", "vector", "vector", "vector")
days <- c("0D", "0D", "2D", "2D", "3D", "3D", "0D", "0D", "2D", "2D", "3D", "3D")
sampleName <- c("F_ZC3H4_0D_1", "F_ZC3H4_0D_2", "F_ZC3H4_2D_1", "F_ZC3H4_2D_2", "F_ZC3H4_3D_1", "F_ZC3H4_3D_2", "Vector_0D_1", "Vector_0D_2", "Vector_2D_1", "Vector_2D_2", "Vector_3D_1", "Vector_3D_2")
sampleTable <- data.frame(treatment, days, sampleName, row.names="sampleName")
sampleTable$treatment <- as.factor(sampleTable$treatment)
sampleTable$days <- as.factor(sampleTable$days)
all(rownames(sampleTable) == colnames(cnts)) #check condition table is correct

##Conduct DESeq2
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cnts,
                              colData = sampleTable,
                              design = ~ treatment + days + treatment:days)
#dds$condition <- relevel(dds$condition, ref = "Vector_0D")#Make WT sample level as reference
#dds$condition <- relevel(dds$condition)#test
dds$treatment <- relevel(dds$treatment, ref = "vector")
dds$days <- relevel(dds$days, ref = "0D")
ddsTP <- DESeq(dds, test="LRT", reduced = ~ treatment + days)
resTP <- results(ddsTP)
head(resTP[order(resTP$log2FoldChange), ], 4)
#dds$condition #check dds

##########################
#Normalized counts for a gene with condition-specific changes over times
data <- plotCounts(ddsTP, which.min(resTP$padj), 
                   intgroup = c("treatment", "days"), returnData = TRUE)
data <- plotCounts(ddsTP, "ENSG00000126016.15",
                   intgroup = c("treatment", "days"), returnData = TRUE)
ggplot(data, aes(x=days, y=count, color=treatment, group=treatment)) + 
  geom_point() + stat_smooth(se=FALSE, method="loess") + scale_y_log10()
#############################

###########################
#Wald tests for the log2 fold changes at individual time points can be investigated using the test argument to results
resultsNames(ddsTP)
res2D <- results(ddsTP, name="treatmentF_ZC3H4.days2D", alpha = 0.05)
summary(res2D)
res2D <- lfcShrink(ddsTP, coef="treatmentF_ZC3H4.days2D", type="apeglm", res=res2D)
plotMA(res2D, ylim=c(-5, 5))
#GIve symbol
res2D_padj <- subset(res2D, padj<0.05)
res2D_f2_padj05 <- subset(res2D_padj, abs(log2FoldChange) > 1)

resultsNames(ddsTP)
res3D <- results(ddsTP, name="treatmentF_ZC3H4.days3D", alpha = 0.05)
summary(res3D)
res3D <- lfcShrink(ddsTP, coef="treatmentF_ZC3H4.days3D", type="apeglm", res=res3D)
plotMA(res3D, ylim=c(-5, 5))
res3D_padj <- subset(res3D, padj<0.05)
res3D_f2_padj05 <- subset(res3D_padj, abs(log2FoldChange) > 1)

#Raw deg
write.csv(as.data.frame(res2D), 
          file = paste(results_dir, "res2D_treatmentF_ZC3H4_days2D_2.csv", sep = ""))

res3D <- results(ddsTP, name="treatmentF_ZC3H4.days3D", alpha = 0.05)
summary(res3D)
res3D <- lfcShrink(ddsTP, coef="treatmentF_ZC3H4.days3D", type="apeglm", res=res3D)
plotMA(res3D, ylim=c(-5, 5))
write.csv(as.data.frame(res3D), 
          file = paste(results_dir, "res2D_treatmentF_ZC3H4_days3D.csv", sep = ""))

###########################
#We can further more cluster significant genes by their profiles, 
# #Load files
# intersect_up <- read.csv(paste(results_dir, "2D_vs_3D_intersect_up.csv", sep=""), header = TRUE)
# intersect_down <- read.csv(paste(results_dir, "2D_vs_3D_intersect_down.csv", sep=""), header = TRUE)
# resTP3 <- resTP
# resTP_up <- resTP3[(rownames(resTP3) %in% intersect_up$x), ]
# resTP_down <- resTP3[(rownames(resTP3) %in% intersect_down$x), ]

#We can further more cluster significant genes by their profiles, 
betas <- coef(ddsTP)
colnames(betas)
library("pheatmap")
mat_up <- betas[(rownames(betas) %in% rownames(resTP_up)), -c(1, 2, 3, 4)]
mat_down <- betas[(rownames(betas) %in% rownames(resTP_down)), -c(1, 2, 3, 4)]
mat <- rbind(mat_up, mat_down)
pheatmap(mat, cluster_cols = FALSE, show_rownames = FALSE, breaks=seq(from=-6, to=6, length=101))
mat_2D <- rbind(mat_up_2D, mat_down_2D)
pheatmap(mat_2D, cluster_cols = FALSE, show_rownames = FALSE, breaks=seq(from=-4.1, to=7, length=101))
mat_3D <- rbind(mat_up_3D, mat_down_3D)
pheatmap(mat_3D, cluster_cols = FALSE, show_rownames = FALSE, breaks=seq(from=-4.1, to=7, length=101))


####################################################################################
##Draw MA plots
resultsNames(dds)
res <- results(dds, contrast = c("condition", "Vector_3D", "Vector_0D"), alpha = 0.05)
resLFC <- lfcShrink(dds, coef="condition_Vector_3D_vs_Vector_0D", type="apeglm", res=res)
plotMA(resLFC, ylim=c(-5, 5))
png(paste(results_dir, "MAplot_Vector_3D_vs_Vector_0D.png", sep = ""), width = 130, height = 130, units='mm', res = 600)
plotMA(resLFC, ylim=c(-5, 5))
dev.off()

summary(resLFC)
#########################################################################################

##assign gene symbol to result####################################################
library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl_id <- c(row.names(resTP))
ensembl_id <- substr(ensembl_id, 1, 15)
res_genesym <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id", 
                     values = ensembl_id, 
                     mart = ensembl)
#remove duplicate -> 여기 if 문 넣어서 0일때는 넘어가게 해야함
if (length(which(duplicated(res_genesym$ensembl_gene_id))) == 0){
  res_genesym <- res_genesym
}else{
  res_genesym <- res_genesym[-which(duplicated(res_genesym$ensembl_gene_id)),]
}

#add empty values to not assigned id in res_genesym
notassigned <- setdiff(ensembl_id, res_genesym$ensembl_gene_id)
empty_col <- rep(" ", time=length(notassigned))
notassigned_df <- data.frame(notassigned, empty_col)
names(notassigned_df) <- c("ensembl_gene_id", "hgnc_symbol")
res_genesym <- rbind(res_genesym, notassigned_df)
#add genesym column to res
resTP$ensembl_gene_id <- c(ensembl_id)
resTP <- resTP[c(order(resTP$ensembl_gene_id)), ]
res_genesym <- res_genesym[c(order(res_genesym$ensembl_gene_id)), ]
resTP$gene_sym <- res_genesym$hgnc_symbol
resTP <- resTP[!(resTP$ensembl_gene_id %in% notassigned), ]

length(resTP$gene_sym)

#res Cutoff
summary(resLFC)
res_padj <- subset(resLFC, padj<0.05)
summary(res_padj)
res_f2_padj05 <- subset(res_padj, abs(log2FoldChange) > 1)
summary(res_f2_padj05)
# res_f2_1.5_padj05 <- subset(res_padj, abs(log2FoldChange) > log2(1.5))
# summary(res_f2_1.5_padj05)

##Save deg_data
#Raw deg
write.csv(as.data.frame(resLFC), 
          file = paste(results_dir, "resLFC_condition_Vector_3D_vs_Vector_0D.csv", sep = ""))
#padj<0.05
write.csv(as.data.frame(res_padj), 
          file = paste(results_dir, "resLFC_padj05_condition_Vector_3D_vs_Vector_0D.csv", sep = ""))
#padj<0.05, foldchange>2
write.csv(as.data.frame(res_f2_padj05), 
          file = paste(results_dir, "resLFC_f2_padj05.csv", sep = ""))
#padj<0.05, foldchange>log2(1.5)
# write.csv(as.data.frame(res_f2_1.5_padj05), 
#           file = "/Users/hmkim/data/buffer/EP400NL/restet_f2_1.5_padj05.csv")

##assign gene symbol to MLL2 result#####################################################

###############################################################################
##Volcano plot MLL2
library(EnhancedVolcano)
res_df <- data.frame(resLFC)
GLUT <- c("SLC2A1", "SLC2A2", "SLC2A3", "SLC2A4", "SLC2A10")
png("/Users/hmkim/data/buffer/PolII/Volcanoplot_PolII.png", 
    width = 260, 
    height = 260, 
    units='mm', 
    res = 300)
EnhancedVolcano(res_df,
                lab = res_df$gene_sym,
                selectLab = "",
                x = 'log2FoldChange',
                y = 'padj', 
                pCutoff = 0.05, 
                FCcutoff = log2(1.5), 
                #title = "",
                #drawConnectors = TRUE, 
                labSize = 6,
                pointSize = 1,
                boxedLabels = FALSE, 
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = FALSE,
                gridlines.major = FALSE, 
                gridlines.minor = FALSE
)

EnhancedVolcano(res_df_test,
                lab = res_df_test$gene_sym,
                selectLab = "",
                x = 'log2FoldChange',
                y = 'padj', 
                pCutoff = 0.05, 
                FCcutoff = log2(1.5), 
                #title = "",
                #drawConnectors = TRUE, 
                labSize = 6,
                pointSize = 1,
                boxedLabels = FALSE, 
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = FALSE,
                gridlines.major = FALSE, 
                gridlines.minor = FALSE
)
dev.off()
###############################################################################

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
#pheatmap(assay(ntd)[select,])
pheatmap(assay(rld)[select,])

##PCA
#sampleDists1 <- dist(t(assay(rld)))
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


################################################################################
##hirachy heatmap
#heatmap with top genes
library("gplots")
library("RColorBrewer")
library("genefilter")
#get deviation(rowVars()) from assay(rld), and get high deviation genes
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 1000)
#cluster gene and sample with pearson correlation
library("pheatmap")
# Pairwise correlation between samples (columns)
cols.cor <- cor(assay(rld)[ topVarGenes, ], use = "pairwise.complete.obs", method = "pearson")
# Pairwise correlation between rows (genes)
rows.cor <- cor(t(assay(rld)[ topVarGenes, ]), use = "pairwise.complete.obs", method = "pearson")
#Draw heatmap
png("/Users/hmkim/data/Cowork/HKB/033121_ClusteringHeatmap_top1000.png", 
    width = 130, 
    height = 130, 
    units='mm', 
    res = 300)
hm <- pheatmap(assay(rld)[ topVarGenes, ], 
               scale = "row", 
               clustering_distance_cols = as.dist(1 - cols.cor), 
               clustering_distance_rows = as.dist(1 - rows.cor), 
               col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), 
               show_rownames = FALSE
)
dev.off()

#heatmap with padj<0.05, abs(foldchange)>2
rld_f2_padj05 <- assay(rld)
rld_f2_padj05 <- rld_f2_padj05[rownames(rld_f2_padj05) %in% rownames(restet_f2_1.5_padj05), ]
#cluster gene and sample with pearson correlation
library("pheatmap")
# Pairwise correlation between samples (columns)
cols.cor <- cor(rld_f2_padj05, use = "pairwise.complete.obs", method = "pearson")
# Pairwise correlation between rows (genes)
rows.cor <- cor(t(rld_f2_padj05), use = "pairwise.complete.obs", method = "pearson")
#Draw heatmap
png("/Users/hmkim/data/buffer/EP400NL/EP400NL_ClusteringHeatmap_f15_padj05.png", 
    width = 200, 
    height = 200, 
    units='mm',
    res = 150)
pheatmap(rld_f2_padj05, 
         scale = "row", 
         clustering_distance_cols = as.dist(1 - cols.cor), 
         clustering_distance_rows = as.dist(1 - rows.cor), 
         col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), 
         show_rownames = FALSE, 
         treeheight_row = 0
)
dev.off()
###############################################################################
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
