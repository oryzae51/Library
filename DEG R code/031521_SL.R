#import data as datafram(df)
##Pull down data, data preprocessing
setwd("/Users/hmkim/data/quant_data/031521_SL/quant_data")

##Automate count data.frame 
#Load data directory
src_dir <- c()
src_dir <- c("/Users/hmkim/data/quant_data/031521_SL/quant_data")
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
names(cnts) <- c("shNT_1", "Klf2_KD_1", "Snd1_KD_1", "Klf2_Snd1_dKD_1",
                 "shNT_2", "Klf2_KD_2", "Snd1_KD_2", "Klf2_Snd1_dKD_2") #change col name



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
condition <- c("shNT", "Klf2_KD", "Snd1_KD", "Klf2_Snd1_dKD",
               "shNT", "Klf2_KD", "Snd1_KD", "Klf2_Snd1_dKD")
sampleName <- c("shNT_1", "Klf2_KD_1", "Snd1_KD_1", "Klf2_Snd1_dKD_1",
                "shNT_2", "Klf2_KD_2", "Snd1_KD_2", "Klf2_Snd1_dKD_2")
sampleTable <- data.frame(condition, sampleName, row.names="sampleName")
sampleTable$condition <- as.factor(sampleTable$condition)
all(rownames(sampleTable) == colnames(cnts)) #check condition table is correct

##Conduct DESeq2
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cnts,
                              colData = sampleTable,
                              design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "shNT")#Make WT sample level as reference
dds <- DESeq(dds)
dds$condition #check dds



##Draw MA plots
#dKO
res_dKD <- results(dds, 
                   contrast = c("condition", "Klf2_Snd1_dKD", "shNT"), 
                   alpha = 0.05)
resultsNames(dds)
resLFC_dKD <- lfcShrink(dds, 
                        coef="condition_Klf2_Snd1_dKD_vs_shNT", 
                        type="apeglm", 
                        res=res_dKD)
plotMA(resLFC_dKD, ylim=c(-8, 8))

#resOrdered_lfc_MLL1 <- resLFC_MLL1[order(resLFC_MLL1$padj),]
#resOrdered_MLL1_pval <- res_MLL1[order(res_MLL1$pvalue),]
#resOrdered_MLL1_lfc <- res1[order(res1$log2FoldChange),]
#summary(res_MLL1)
summary(resLFC_dKD)
resLFC_dKD_padj <- subset(resLFC_dKD, padj<0.05)
summary(resLFC_dKD_padj)
resLFC_dKD_up <- subset(resLFC_dKD, log2FoldChange > 0)
summary(resLFC_dKD_up)
resLFC_dKD_down <- subset(resLFC_dKD, log2FoldChange < 0)
summary(resLFC_dKD_down)
#write.csv(as.data.frame(resOrdered1p), 
#          file="~/data/sigtest1.csv")
write.csv(as.data.frame(resLFC_dKD), 
          file = "~/data/resLFC_dKD.csv")
write.csv(as.data.frame(resLFC_dKD), 
          file = "~/data/resLFC_dKD_old.csv")
write.csv(as.data.frame(resLFC_dKD_down), 
          file = "~/data/resLFC_dKD_old_down.csv")
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


#Klf2-KD
res_Klf2 <- results(dds, 
                    contrast = c("condition", "Klf2_KD", "shNT"), 
                    alpha = 0.05)
resultsNames(dds)
resLFC_Klf2 <- lfcShrink(dds, 
                         coef="condition_Klf2_KD_vs_shNT", 
                         type="apeglm", 
                         res=res_Klf2)
plotMA(resLFC_Klf2, ylim=c(-5, 5))
#resOrdered_lfc_MLL2 <- resLFC_MLL2[order(resLFC_MLL2$padj),]
#resOrdered_MLL1_pval <- res_MLL1[order(res_MLL1$pvalue),]
#resOrdered_MLL1_lfc <- res1[order(res1$log2FoldChange),]
#summary(res_MLL1)
summary(resLFC_Klf2)
resLFC_Klf2_padj <- subset(resLFC_Klf2, padj<0.05)
summary(resLFC_Klf2_padj)
resLFC_Klf2_fc_up <- subset(resLFC_Klf2_padj, log2FoldChange > 0)
resLFC_Klf2_fc_down <- subset(resLFC_Klf2_padj, log2FoldChange < 0)
resLFC_Klf2_fc125up <- subset(resLFC_Klf2_padj, log2FoldChange > 0.32)
summary(resLFC_Klf2_fc125up)
resLFC_Klf2_fc125down <- subset(resLFC_Klf2_padj, log2FoldChange < -0.32)
summary(resLFC_Klf2_fc125down)
#write.csv(as.data.frame(resOrdered1p), 
#          file="~/data/sigtest1.csv")
write.csv(as.data.frame(resLFC_Klf2_fc_up), 
          file = "/Users/hmkim/data/Cowork/SL/DEG_list/fc125/resLFC_Klf2_fc_up.csv")
write.csv(as.data.frame(resLFC_Klf2_fc_down), 
          file = "/Users/hmkim/data/Cowork/SL/DEG_list/fc125/resLFC_Klf2_fc_down.csv")

write.csv(as.data.frame(resLFC_Klf2_fc125up), 
          file = "/Users/hmkim/data/Cowork/SL/DEG_list/fc125/resLFC_Klf2_fc125up.csv")
write.csv(as.data.frame(resLFC_Klf2_fc125down), 
          file = "/Users/hmkim/data/Cowork/SL/DEG_list/fc125/resLFC_Klf2_fc125down.csv")


#resSig1 <- subset(res1, padj < 0.01)
#resSig1 <- subset(resSig1, log2FoldChange>1 | log2FoldChange< -1)
#resSig1 <- subset(resSig1, log2FoldChange>1)


#MLL4-KO
res_Snd1 <- results(dds, 
                    contrast = c("condition", "Snd1_KD", "shNT"), 
                    alpha = 0.05)
resultsNames(dds)
resLFC_Snd1 <- lfcShrink(dds, 
                         coef="condition_Snd1_KD_vs_shNT", 
                         type="apeglm", 
                         res=res_Snd1)
plotMA(resLFC_Snd1, ylim=c(-5, 5))
#resOrdered_lfc_MLL4 <- resLFC_MLL4[order(resLFC_MLL4$padj),]
#resOrdered_MLL1_pval <- res_MLL1[order(res_MLL1$pvalue),]
#resOrdered_MLL1_lfc <- res1[order(res1$log2FoldChange),]
#summary(res_MLL1)
summary(resLFC_Snd1)
resLFC_Snd1_padj <- subset(resLFC_Snd1, padj<0.05)
summary(resLFC_Snd1_padj)
resLFC_Snd1_fc_up <- subset(resLFC_Snd1_padj, log2FoldChange > 0)
resLFC_Snd1_fc_down <- subset(resLFC_Snd1_padj, log2FoldChange < 0)
resLFC_Snd1_fc125up <- subset(resLFC_Snd1_padj, log2FoldChange > 0.32)
summary(resLFC_Snd1_fc125up)
resLFC_Snd1_fc125down <- subset(resLFC_Snd1_padj, log2FoldChange < -0.32)
summary(resLFC_Snd1_fc125down)
#write.csv(as.data.frame(resOrdered1p), 
#          file="~/data/sigtest1.csv")
write.csv(as.data.frame(resLFC_Snd1_fc_up), 
          file = "/Users/hmkim/data/Cowork/SL/DEG_list/fc125/resLFC_Snd1_fc_up.csv")
write.csv(as.data.frame(resLFC_Snd1_fc_down), 
          file = "/Users/hmkim/data/Cowork/SL/DEG_list/fc125/resLFC_Snd1_fc_down.csv")

write.csv(as.data.frame(resLFC_Snd1_fc125up), 
          file = "/Users/hmkim/data/Cowork/SL/DEG_list/fc125/resLFC_Snd1_fc125up.csv")
write.csv(as.data.frame(resLFC_Snd1_fc125down), 
          file = "/Users/hmkim/data/Cowork/SL/DEG_list/fc125/resLFC_Snd1_fc125down.csv")

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


##################################################3
##venn diagram
resLFC_Klf2_up_df <- as.data.frame(resLFC_Klf2_up)
Klf2_down_geneid <- rownames(resLFC_Klf2_down)
Klf2_up_geneid <- rownames(resLFC_Klf2_up)
Snd1_down_geneid <- rownames(resLFC_Snd1_down)
Snd1_up_geneid <- rownames(resLFC_Snd1_up)

library(gplots)
data_up <- list("shKlf2 up-regulated genes (n=4819)" = Klf2_up_geneid, "shSnd1 up-rgulated genes (n=4018)" = Snd1_up_geneid)
png("~/klf2_Snd1_up_venn.png", width = 130, height = 130, units='mm', res = 1200)
venn(data_up)
dev.off()

data_down <- list("shKlf2 down-regulated genes (n=4949)" = Klf2_down_geneid, "shSnd1 down-regulated genes (n=4166)" = Snd1_down_geneid)
png("~/klf2_Snd1_down_venn.png", width = 130, height = 130, units='mm', res = 1200)
venn(data_down)
dev.off()

# MLLintersect <- intersect(c(resLFC_MLL1df$Geneid), c(resLFC_MLL2df$Geneid))
# MLLintersect <- intersect(MLLintersect, c(resLFC_MLL4df$Geneid))
# resLFC_MLLKO_padj01_intersect <- subset(resLFC_MLL1df, Geneid %in% MLLintersect)
# write.csv(as.data.frame(MLL_KO_intersect_down), 
#           file = "~/data/deg_data/public_data_confirm/MLL_KO_intersect_down.csv")
# write.csv(as.data.frame(MLL_KO_intersect_up), 
#           file = "~/data/deg_data/public_data_confirm/MLL_KO_intersect_up.csv")

##venn diagram end
################################################


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
sampleDists <- dist(t(assay(rld)))
sampleDists <- dist(scale(t(assay(rld))))#recommand this for consistancy
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, sep="-")
#colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png("/Users/hmkim/data/Cowork/SL/distance_heatmap.png", width = 198, height = 170, units='mm', res = 1200)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

sampleDists <- dist(scale(t(assay(rld))))#recommand this for consistancy
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, sep="-")
hrsDM <- hclust(as.dist(1-cor(t(sampleDistMatrix), method="spearman")), 
                method="complete") # Cluster rows by Pearson correlation.
hcsDM <- hclust(as.dist(1-cor(t(sampleDistMatrix), method="spearman")), 
                method="complete") # Clusters columns by Spearman correlation.
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

png("/Users/hmkim/data/Cowork/SL/032621_Distance_heatmap_spearman.png", width = 198, height = 170, units='mm', res = 1200)
heatmap.2(sampleDistMatrix,
          Rowv=as.dendrogram(hrsDM), 
          Colv=as.dendrogram(hcsDM),
          scale = "none",
          trace = "none", 
          dendrogram = "both", 
          col = colors, 
          cexRow = 1,
          cexCol = 0.8, margins=c(7, 5),
          density.info="none", 
          lhei=c(1.5, 7)
)
dev.off()

png("/Users/hmkim/data/Cowork/SL/PCA_SL.png", width = 198, height = 170, units='mm', res = 1200)
plotPCA(rld, intgroup=c("condition"))
dev.off()


##hirachy heatmap
library("gplots")
library("RColorBrewer")
library("genefilter")
#get deviation(rowVars()) from assay(rld), and get high deviation genes
# topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), 50)
# topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), 500)
topVarGenes <- head(assay(rld), decreasing=TRUE, 1000)
topVarGenes <- head(order(assay(rld), decreasing=TRUE), 1000)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 50)
write.csv(as.data.frame(topVarGenes), 
          file = "/Users/hmkim/data/Cowork/SL/032921 topVarGenes_1k.csv")
#scaling dataset
scaledata <- t(scale(t(topVarGenes))) # Centers and scales data.
scaledata <- scaledata[complete.cases(scaledata),]

rldTable <- assay(rld)
scaledata <- t(scale(t(rldTable)))
scaledata <- scaledata[complete.cases(scaledata),]
#clustering row(gene) and column(sample) with correlation analysis
hr <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete") # Clusters columns by Spearman correlation.
#pdf("~/data/asdfasdfasdftest.pdf")#For export. Not necessary if you are using Rstudio
rowlab <- read.csv("~/data/rowlab.csv", header = FALSE)
rowlab <- rowlab$V1
png("~/testtest.png", width = 198, height = 170, units='mm', res = 1200)
heatmap.2(rldTable,
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hc),
          scale="row",
          trace="none",
          dendrogram="col",
          cexRow = 1,
          cexCol = 0.8,
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255), 
          labRow = rowlab,
          margins=c(7, 5), 
          density.info="none", 
          lhei=c(1.5, 7)
)

heatmap.2(assay(rld)[ topVarGenes, ],
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hc),
          scale="row",
          trace="none",
          dendrogram="col",
          cexRow = 1,
          cexCol = 0.8,
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255), 
          density.info="none"
)
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
           )
hm <- heatmap.2( assay(rld)[ topVarGenes, ], scale="row", 
                 trace="none", dendrogram="both", 
                 col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
)
hcluster <- as.hclust( hm$rowDendrogram )
cutr_hcluster<-cutree( hcluster, h=20 )
write.csv(as.data.frame(cutr_hcluster), 
          file = "/Users/hmkim/data/Cowork/SL/cutr_hcluster.csv")
dev.off()
asdf <- assay(rld)[ topVarGenes, ]
write.csv(as.data.frame(asdf), 
          file = "/Users/hmkim/data/Cowork/SL/032921 topVarGenes_1k_2.csv")

test <- row.names(rldTable)

##########################################
####Heatmap-pvalue
rldTable_df <- as.data.frame(rldTable)
head(rldTable_df)
#Keep only the significantly differentiated genes
resLFC_dKD_df <- as.data.frame(resLFC_dKD)
sigGenes <- rownames(resLFC_dKD_df[resLFC_dKD_df$padj <= 0.05 & abs(resLFC_dKD_df$log2FoldChange>1),])
rldTable_df <- rldTable_df[rldTable_df$]


##########################################
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
