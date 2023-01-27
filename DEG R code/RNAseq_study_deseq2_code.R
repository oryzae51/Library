##Pull down data, data preprocessing
#import data as datafram(df)
setwd("/Users/hmkim/data/quant_data/MLL3_star_confirm_082320")
MLL3_1_1df=read.table("SRR5664471.1.fastqAligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
MLL3_1_2df=read.table("SRR5664472.1.fastqAligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
MLL3_2_1df=read.table("SRR5664473.1.fastqAligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
MLL3_2_2df=read.table("SRR5664474.1.fastqAligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
MLL3_3_1df=read.table("SRR5664475.1.fastqAligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
MLL3_3_2df=read.table("SRR5664476.1.fastqAligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
WT_1df=read.table("SRR5664477.1.fastqAligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
WT_2df=read.table("SRR5664478.1.fastqAligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
MLL3_1_1rn=read.table("SRR5664471.1.fastqAligned.out.sam.txt", sep="\t", header=TRUE)#for extranct geneid

#extract read count column and integrate all samples counts in cnts
MLL3_1_1 <- c(MLL3_1_1df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL3.SAli_star_hg38.SRR5664471.1.fastqAligned.out.sam)
MLL3_1_2 <- c(MLL3_1_2df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL3.SAli_star_hg38.SRR5664472.1.fastqAligned.out.sam)
MLL3_2_1 <- c(MLL3_2_1df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL3.SAli_star_hg38.SRR5664473.1.fastqAligned.out.sam)
MLL3_2_2 <- c(MLL3_2_2df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL3.SAli_star_hg38.SRR5664474.1.fastqAligned.out.sam)
MLL3_3_1 <- c(MLL3_3_1df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL3.SAli_star_hg38.SRR5664475.1.fastqAligned.out.sam)
MLL3_3_2 <- c(MLL3_3_2df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL3.SAli_star_hg38.SRR5664476.1.fastqAligned.out.sam)
WT_1 <- c(WT_1df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL3.SAli_star_hg38.SRR5664477.1.fastqAligned.out.sam)
WT_2 <- c(WT_2df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL3.SAli_star_hg38.SRR5664478.1.fastqAligned.out.sam)
Geneid <- c(MLL3_1_1rn$Geneid)#save Geneid
#make count matrix cnts
cnts<-data.frame(MLL3_1_1, MLL3_1_2, MLL3_2_1, MLL3_2_2, MLL3_3_1, MLL3_3_2, WT_1, WT_2, Geneid, row.names = "Geneid")
cnts <- cnts[rowSums(cnts > 1) >=8,]#drop genes with low counts, this is necessary for rld quality

#make condition table for downstream processing
condition <- c("MLL3_1", "MLL3_1", "MLL3_2", "MLL3_2", "MLL3_3", "MLL3_3", "WT", "WT")
sampleName <- c("MLL3_1_1", "MLL3_1_2", "MLL3_2_1", "MLL3_2_2", "MLL3_3_1", "MLL3_3_2", "WT_1", "WT_2")
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
res_MLL3_1 <- results(dds, contrast = c("condition", "MLL3_1", "WT"), alpha = 0.01)
resultsNames(dds)
plotMA(res_MLL3_1, ylim = c(-5, 5))
resLFC_MLL3_1 <- lfcShrink(dds, coef="condition_MLL3_1_vs_WT", type="apeglm", res=res_MLL3_1)
plotMA(resLFC_MLL3_1, ylim = c(-5, 5))

res_MLL3_2 <- results(dds, contrast = c("condition", "MLL3_2", "WT"), alpha = 0.01)
resultsNames(dds)
resLFC_MLL3_2 <- lfcShrink(dds, coef="condition_MLL3_2_vs_WT", type="apeglm", res=res_MLL3_2)
plotMA(resLFC_MLL3_2, ylim = c(-5, 5))

res_MLL3_3 <- results(dds, contrast = c("condition", "MLL3_3", "WT"), alpha = 0.01)
resultsNames(dds)
resLFC_MLL3_3 <- lfcShrink(dds, coef="condition_MLL3_3_vs_WT", type="apeglm", res=res_MLL3_3)
plotMA(resLFC_MLL3_3, ylim = c(-5, 5))


##Save the DESeq2 result
write.csv(as.data.frame(resLFC_MLL3_1), file = "/Users/hmkim/rseqstudy/resLFC_MLL3_1.csv")
write.csv(as.data.frame(resLFC_MLL3_2), file = "/Users/hmkim/rseqstudy/resLFC_MLL3_2.csv")
write.csv(as.data.frame(resLFC_MLL3_3), file = "/Users/hmkim/rseqstudy/resLFC_MLL3_3.csv")


##Extract DESeq2 result with condition
summary(resLFC_MLL3_1)
summary(resLFC_MLL3_2)
#subset sample adjusted p value
resLFC_MLL3_1_padj <- subset(resLFC_MLL3_1, padj<0.01)
summary(resLFC_MLL3_1_padj)
resLFC_MLL3_1_padj_twofold <- subset(resLFC_MLL3_1_padj, log2FoldChange < -1 | log2FoldChange > 1)
summary(resLFC_MLL3_1_padj_twofold)

resLFC_MLL3_2_padj_twofold <- subset(resLFC_MLL3_2, (log2FoldChange < -1 | log2FoldChange > 1) & padj < 0.01 )
summary(resLFC_MLL3_2_padj_twofold)

##Count data transformation/normalization for downstream analyze
#Normalizae with rlog in DESeq2
library("vsn")
rld <- rlog(dds, blind=FALSE)
meanSdPlot(assay(rld))


##PCA(Principle component analysis)
sampleDists <- dist(t(assay(rld)))
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


##hierachy heatmap
#get deviation(rowVars()) from assay(rld), and get high deviation genes
library("gplots")
library("RColorBrewer")
library("genefilter")
topVarGenes <- head(assay(rld), decreasing=TRUE, 500)
#scaling dataset
scaledata <- t(scale(t(topVarGenes))) # Centers and scales data.
scaledata <- scaledata[complete.cases(scaledata),]
#clustering row(gene) and column(sample) with correlation analysis
hr <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete") # Clusters columns by Spearman correlation.
heatmap.2(topVarGenes,
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hc),
          scale="row",
          trace="none",
          dendrogram="both",
          cexRow = 0.4,
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255)
)

library(EnhancedVolcano)
MLL3_1_vol=read.table("/Users/hmkim/rseqstudy/resLFC_MLL3_1.csv", sep=",", header=TRUE)
EnhancedVolcano(MLL3_1_vol,
                lab = rownames(MLL3_1_vol),
                x = 'log2FoldChange',
                y = 'padj', 
                pCutoff = 0.01, 
                FCcutoff = 1, 
                title = "MLL3-KO",
                labSize = 0.0,
                boxedLabels = FALSE
)
