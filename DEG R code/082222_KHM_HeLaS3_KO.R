#import data as datafram(df)
##Pull down data, data preprocessing
cnt_dir <- "/Users/kimhyoungmin/data/082222_KHM_HeLa/quant/cnt"
results_dir <- "/Users/kimhyoungmin/Library/Mobile Documents/com~apple~CloudDocs/ETL/Experiment\ data/2022\ data/082222_KHM_HeLaS3_KO_results/"
setwd("/Users/kimhyoungmin/data/082222_KHM_HeLa/quant/cnt")

##Automate count data.frame 
#Load data directory
src_dir <- c()
src_dir <- c("/Users/kimhyoungmin/data/082222_KHM_HeLa/quant/cnt")
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
names(cnts) <- c("FK1_1", "FK1_2", "FK2_1", "FK2_2", "FKD_1", "FKD_2", "MLL2_1", "MLL2_2", "WT_HL_1", "WT_HL_2") #change col name

#make condition table for downstream processing
condition <- c("FK1", "FK1", "FK2", "FK2", "FKD", "FKD", "MLL2", "MLL2", "WT_HL", "WT_HL")
sampleName <- c("FK1_1", "FK1_2", "FK2_1", "FK2_2", "FKD_1", "FKD_2", "MLL2_1", "MLL2_2", "WT_HL_1", "WT_HL_2")
sampleTable <- data.frame(condition, sampleName, row.names="sampleName")
sampleTable$condition <- as.factor(sampleTable$condition)
all(rownames(sampleTable) == colnames(cnts)) #check condition table is correct

##Conduct DESeq2
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cnts,
                              colData = sampleTable,
                              design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "WT_HL")#Make WT sample level as reference
dds <- DESeq(dds)
dds$condition #check dds

##PCA plot
library("vsn")
rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup=c("condition"))

##Draw MA plots - FK1
resultsNames(dds)
res_FK1 <- results(dds, contrast = c("condition", "FK1", "WT_HL"), alpha = 0.05)
resLFC_FK1 <- lfcShrink(dds, coef="condition_FK1_vs_WT_HL", type="apeglm", res=res_FK1)
plotMA(resLFC_FK1, ylim=c(-5, 5))
# png(paste(results_dir, "MAplot_Vector_3D_vs_Vector_0D.png", sep = ""), width = 130, height = 130, units='mm', res = 600)
# plotMA(resLFC, ylim=c(-5, 5))
# dev.off()
summary(resLFC_FK1)

##Draw MA plots - FK2
resultsNames(dds)
res_FK2 <- results(dds, contrast = c("condition", "FK2", "WT_HL"), alpha = 0.05)
resLFC_FK2 <- lfcShrink(dds, coef="condition_FK2_vs_WT_HL", type="apeglm", res=res_FK2)
plotMA(resLFC_FK2, ylim=c(-5, 5))
# png(paste(results_dir, "MAplot_Vector_3D_vs_Vector_0D.png", sep = ""), width = 130, height = 130, units='mm', res = 600)
# plotMA(resLFC, ylim=c(-5, 5))
# dev.off()
summary(resLFC_FK2)

##Draw MA plots - FKD
resultsNames(dds)
res_FKD <- results(dds, contrast = c("condition", "FKD", "WT_HL"), alpha = 0.05)
resLFC_FKD <- lfcShrink(dds, coef="condition_FKD_vs_WT_HL", type="apeglm", res=res_FKD)
plotMA(resLFC_FKD, ylim=c(-5, 5))
# png(paste(results_dir, "MAplot_Vector_3D_vs_Vector_0D.png", sep = ""), width = 130, height = 130, units='mm', res = 600)
# plotMA(resLFC, ylim=c(-5, 5))
# dev.off()
summary(resLFC_FKD)

##Draw MA plots - MLL2
resultsNames(dds)
res_MLL2 <- results(dds, contrast = c("condition", "MLL2", "WT_HL"), alpha = 0.05)
resLFC_MLL2 <- lfcShrink(dds, coef="condition_MLL2_vs_WT_HL", type="apeglm", res=res_MLL2)
plotMA(resLFC_MLL2, ylim=c(-5, 5))
# png(paste(results_dir, "MAplot_Vector_3D_vs_Vector_0D.png", sep = ""), width = 130, height = 130, units='mm', res = 600)
# plotMA(resLFC, ylim=c(-5, 5))
# dev.off()
summary(resLFC_MLL2)

###########################################################
#res Cutoff
res_padj_FK1 <- subset(resLFC_FK1, padj<0.05)
summary(res_padj_FK1)
res_f2_padj05_FK1 <- subset(res_padj_FK1, abs(log2FoldChange) > 1)
summary(res_f2_padj05_FK1)

res_padj_FK2 <- subset(resLFC_FK2, padj<0.05)
summary(res_padj_FK2)
res_f2_padj05_FK2 <- subset(res_padj_FK2, abs(log2FoldChange) > 1)
summary(res_f2_padj05_FK2)

res_padj_FKD <- subset(resLFC_FKD, padj<0.05)
summary(res_padj_FKD)
res_f2_padj05_FKD <- subset(res_padj_FKD, abs(log2FoldChange) > 1)
summary(res_f2_padj05_FKD)

res_padj_MLL2 <- subset(resLFC_MLL2, padj<0.05)
summary(res_padj_MLL2)
res_f2_padj05_MLL2 <- subset(res_padj_MLL2, abs(log2FoldChange) > 1)
summary(res_f2_padj05_MLL2)

res_f15_padj05_FK1 <- subset(res_padj_FK1, abs(log2FoldChange) > log2(1.5))
res_f15_padj05_FK2 <- subset(res_padj_FK2, abs(log2FoldChange) > log2(1.5))
res_f15_padj05_FKD <- subset(res_padj_FKD, abs(log2FoldChange) > log2(1.5))
res_f15_padj05_MLL2 <- subset(res_padj_MLL2, abs(log2FoldChange) > log2(1.5))

########Save res
#Raw deg
write.csv(as.data.frame(resLFC_FK1), 
          file = paste(results_dir, "resLFC_FK1.csv", sep = ""))
#padj<0.05, foldchange>2
write.csv(as.data.frame(res_f2_padj05_FK1), 
          file = paste(results_dir, "res_f2_padj05_FK1.csv", sep = ""))

#Raw deg
write.csv(as.data.frame(resLFC_FK2), 
          file = paste(results_dir, "resLFC_FK2.csv", sep = ""))
#padj<0.05, foldchange>2
write.csv(as.data.frame(res_f2_padj05_FK2), 
          file = paste(results_dir, "res_f2_padj05_FK2.csv", sep = ""))

#Raw deg
write.csv(as.data.frame(resLFC_FKD), 
          file = paste(results_dir, "resLFC_FKD.csv", sep = ""))
#padj<0.05, foldchange>2
write.csv(as.data.frame(res_f2_padj05_FKD), 
          file = paste(results_dir, "res_f2_padj05_FKD.csv", sep = ""))

#Raw deg
write.csv(as.data.frame(resLFC_MLL2), 
          file = paste(results_dir, "resLFC_MLL2.csv", sep = ""))
#padj<0.05, foldchange>2
write.csv(as.data.frame(res_f2_padj05_MLL2), 
          file = paste(results_dir, "res_f2_padj05_MLL2.csv", sep = ""))

write.csv(as.data.frame(res_f15_padj05_FK1), file = paste(results_dir, "res_f15_padj05_FK1.csv", sep = ""))
write.csv(as.data.frame(res_f15_padj05_FK2), file = paste(results_dir, "res_f15_padj05_FK2.csv", sep = ""))
write.csv(as.data.frame(res_f15_padj05_FKD), file = paste(results_dir, "res_f15_padj05_FKD.csv", sep = ""))
write.csv(as.data.frame(res_f15_padj05_MLL2), file = paste(results_dir, "res_f15_padj05_MLL2.csv", sep = ""))

write.csv(as.data.frame(res_padj_FK1), file = paste(results_dir, "res_padj05_FK1.csv", sep = ""))
write.csv(as.data.frame(res_padj_FK2), file = paste(results_dir, "res_padj05_FK2.csv", sep = ""))
write.csv(as.data.frame(res_padj_FKD), file = paste(results_dir, "res_padj05_FKD.csv", sep = ""))
write.csv(as.data.frame(res_padj_MLL2), file = paste(results_dir, "res_padj05_MLL2.csv", sep = ""))
