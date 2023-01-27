#import data as datafram(df)
##Pull down data, data preprocessing
cnt_dir <- "/Users/kimhyoungmin/data/082222_ADJ_PRMT/quant/cnt"
results_dir <- "/Users/kimhyoungmin/Library/Mobile Documents/com~apple~CloudDocs/ETL/Experiment\ data/Cowork/082222_ADJ_KO/low_test/"
setwd("/Users/kimhyoungmin/data/082222_ADJ_PRMT/quant/cnt")

##Automate count data.frame 
#Load data directory
src_dir <- c()
src_dir <- c(cnt_dir)
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
#cnts <- cnts[rowSums(cnts > 1) >=length(cnts),]#drop genes with low counts, this is necessary for rld quality
names(cnts) <- c("NT_M2_1", "NT_M2_2", "PRMT_M2_1", "PRMT_M2_2", "PRMT_NT_1", "PRMT_NT_2", "WT_293_1", "WT_293_2") #change col name

#make condition table for downstream processing
condition <- c("NT_M2", "NT_M2", "PRMT_M2", "PRMT_M2", "PRMT_NT", "PRMT_NT", "WT_293", "WT_293")
sampleName <- c("NT_M2_1", "NT_M2_2", "PRMT_M2_1", "PRMT_M2_2", "PRMT_NT_1", "PRMT_NT_2", "WT_293_1", "WT_293_2")
sampleTable <- data.frame(condition, sampleName, row.names="sampleName")
sampleTable$condition <- as.factor(sampleTable$condition)
all(rownames(sampleTable) == colnames(cnts)) #check condition table is correct

##Conduct DESeq2
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cnts,
                              colData = sampleTable,
                              design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "WT_293")#Make WT sample level as reference
####pre-filtering the dataset
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)
####
dds <- DESeq(dds)
dds$condition #check dds

##PCA plot
library("vsn")
rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup=c("condition"))

##Draw MA plots - NT_M2
resultsNames(dds)
res_NT_M2 <- results(dds, contrast = c("condition", "NT_M2", "WT_293"), alpha = 0.05)
resLFC_NT_M2 <- lfcShrink(dds, coef="condition_NT_M2_vs_WT_293", type="apeglm", res=res_NT_M2)
plotMA(resLFC_NT_M2, ylim=c(-5, 5))
# png(paste(results_dir, "MAplot_Vector_3D_vs_Vector_0D.png", sep = ""), width = 130, height = 130, units='mm', res = 600)
# plotMA(resLFC, ylim=c(-5, 5))
# dev.off()
summary(resLFC_NT_M2)

##Draw MA plots -  PRMT_NT
resultsNames(dds)
res_PRMT_NT <- results(dds, contrast = c("condition", "PRMT_NT", "WT_293"), alpha = 0.05)
resLFC_PRMT_NT <- lfcShrink(dds, coef="condition_PRMT_NT_vs_WT_293", type="apeglm", res=res_PRMT_NT)
plotMA(resLFC_PRMT_NT, ylim=c(-5, 5))
# png(paste(results_dir, "MAplot_Vector_3D_vs_Vector_0D.png", sep = ""), width = 130, height = 130, units='mm', res = 600)
# plotMA(resLFC, ylim=c(-5, 5))
# dev.off()
summary(resLFC_PRMT_NT)

##Draw MA plots -  PRMT_M2
resultsNames(dds)
res_PRMT_M2 <- results(dds, contrast = c("condition", "PRMT_M2", "WT_293"), alpha = 0.05)
resLFC_PRMT_M2_asdf <- lfcShrink(dds, coef="condition_PRMT_M2_vs_WT_293", type="apeglm", res=res_PRMT_M2)
plotMA(resLFC_PRMT_M2, ylim=c(-5, 5))
# png(paste(results_dir, "MAplot_Vector_3D_vs_Vector_0D.png", sep = ""), width = 130, height = 130, units='mm', res = 600)
# plotMA(resLFC, ylim=c(-5, 5))
# dev.off()
summary(resLFC_PRMT_M2_asdf)
###########################################################
#res Cutoff
res_padj_NT_M2 <- subset(resLFC_NT_M2, padj<0.05)
res_f2_padj05_NT_M2 <- subset(res_padj_NT_M2, abs(log2FoldChange) > 1)

res_padj_PRMT_NT <- subset(resLFC_PRMT_NT, padj<0.05)
res_f2_padj05_PRMT_NT <- subset(res_padj_PRMT_NT, abs(log2FoldChange) > 1)

res_padj_PRMT_M2 <- subset(resLFC_PRMT_M2, padj<0.05)
res_f2_padj05_PRMT_M2 <- subset(res_padj_PRMT_M2, abs(log2FoldChange) > 1)

res_f15_padj05_NT_M2 <- subset(res_padj_NT_M2, abs(log2FoldChange) > log2(1.5))
#summary(res_f15_padj05_FK1)
res_f15_padj05_PRMT_NT <- subset(res_padj_PRMT_NT, abs(log2FoldChange) > log2(1.5))
#summary(res_f15_padj05_FK2)
res_f15_padj05_PRMT_M2 <- subset(res_padj_PRMT_M2, abs(log2FoldChange) > log2(1.5))
#summary(res_f15_padj05_FKD)


########Save res
#Raw deg
write.csv(as.data.frame(resLFC_NT_M2), file = paste(results_dir, "resLFC_NT_M2.csv", sep = ""))
#padj<0.05, foldchange>2
write.csv(as.data.frame(res_f2_padj05_NT_M2), file = paste(results_dir, "res_f2_padj05_NT_M2.csv", sep = ""))

#Raw deg
write.csv(as.data.frame(resLFC_PRMT_NT), file = paste(results_dir, "resLFC_PRMT_NT.csv", sep = ""))
#padj<0.05, foldchange>2
write.csv(as.data.frame(res_f2_padj05_PRMT_NT), file = paste(results_dir, "res_f2_padj05_PRMT_NT.csv", sep = ""))

#Raw deg
write.csv(as.data.frame(resLFC_PRMT_M2), file = paste(results_dir, "resLFC_PRMT_M2.csv", sep = ""))
#padj<0.05, foldchange>2
write.csv(as.data.frame(res_f2_padj05_PRMT_M2), file = paste(results_dir, "res_f2_padj05_PRMT_M2.csv", sep = ""))

write.csv(as.data.frame(res_f15_padj05_NT_M2), file = paste(results_dir, "res_f15_padj05_NT_M2.csv", sep = ""))
write.csv(as.data.frame(res_f15_padj05_PRMT_NT), file = paste(results_dir, "res_f15_padj05_PRMT_NT.csv", sep = ""))
write.csv(as.data.frame(res_f15_padj05_PRMT_M2), file = paste(results_dir, "res_f15_padj05_PRMT_M2.csv", sep = ""))

write.csv(as.data.frame(res_padj_NT_M2), file = paste(results_dir, "res_padj05_NT_M2.csv", sep = ""))
write.csv(as.data.frame(res_padj_PRMT_NT), file = paste(results_dir, "res_padj05_PRMT_NT.csv", sep = ""))
write.csv(as.data.frame(res_padj_PRMT_M2), file = paste(results_dir, "res_padj05_PRMT_M2.csv", sep = ""))
