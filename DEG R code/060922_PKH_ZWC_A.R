#make conditionA table for downstream processing
conditionA <- c("F_ZC3H4_0D", "F_ZC3H4_0D", "Vector_0D", "Vector_0D")
sampleNameA <- c("F_ZC3H4_0D_1", "F_ZC3H4_0D_2", "Vector_0D_1", "Vector_0D_2")
sampleTableA <- data.frame(conditionA, sampleNameA, row.names="sampleNameA")
sampleTableA$conditionA <- as.factor(sampleTableA$conditionA)
all(rownames(sampleTableA) == colnames(cntsA)) #check conditionA table is correct

##Conduct DESeq2
library("DESeq2")
ddsA <- DESeqDataSetFromMatrix(countData = cntsA,
                                colData = sampleTableA,
                                design = ~ conditionA)
ddsA$conditionA <- relevel(ddsA$conditionA, ref = "Vector_0D")#Make WT sample level as reference
ddsA <- DESeq(ddsA)
ddsA$conditionA #check dds


####################################################################################
##Draw MA plots
resA <- results(ddsA, contrast = c("conditionA", "F_ZC3H4_0D", "Vector_0D"), alpha = 0.05)
resultsNames(ddsA)
resLFCA <- lfcShrink(ddsA, coef="conditionA_F_ZC3H4_0D_vs_Vector_0D", type="apeglm", res=resA)
plotMA(resLFCA, ylim=c(-5, 5))
png(paste(results_dir, "MAplot_vsRing1bKD.png", sep = ""), width = 130, height = 130, units='mm', res = 600)
plotMA(resLFCA, ylim=c(-5, 5))
dev.off()

summary(resLFCA)

#########################################################################
##assign gene symbol to result####################################################
library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl_id <- c(row.names(resLFCA))
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
resLFCA$ensembl_gene_id <- c(ensembl_id)
resLFCA <- resLFCA[c(order(resLFCA$ensembl_gene_id)), ]
res_genesym <- res_genesym[c(order(res_genesym$ensembl_gene_id)), ]
resLFCA$gene_sym <- res_genesym$hgnc_symbol
resLFCA <- resLFCA[!(resLFCA$ensembl_gene_id %in% notassigned), ]

#res Cutoff
summary(resLFCA)
res_padjA <- subset(resLFCA, padj<0.05)
summary(res_padjA)
res_f2_padj05A <- subset(res_padjA, abs(log2FoldChange) > 1)
summary(res_f2_padj05A)

#Raw deg
write.csv(as.data.frame(resLFCA), 
          file = paste(results_dir, "resLFCA.csv", sep = ""))
#padj<0.05, foldchange>2
write.csv(as.data.frame(res_f2_padj05A), 
          file = paste(results_dir, "resLFC_f2_padj05A.csv", sep = ""))
#foldchange>2
write.csv(as.data.frame(res_padjA), 
          file = paste(results_dir, "res_padjA.csv", sep = ""))

###############################################################################

##Count data transformation/normalization for downstream analyze
#Normalizae with rlog in DESeq2
library("vsn")
rldA <- rlog(ddsA, blind=FALSE)
meanSdPlot(assay(rldA))
#Compare with normTrnaform. You don't need to do when you are making data. 
ntdA <- normTransform(ddsA)
meanSdPlot(assay(ntdA))


##PCA
#sampleDistsA1 <- dist(t(assay(rldA)))
sampleDistsA <- dist(scale(t(assay(rldA))))#recommand this for consistancy
library("RColorBrewer")
sampleDistMatrixA <- as.matrix(sampleDistsA)
rownames(sampleDistMatrixA) <- paste(rldA$conditionA, sep="-")
#colnames(sampleDistMatrixA) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrixA,
         clustering_distance_rows=sampleDistsA,
         clustering_distance_cols=sampleDistsA,
         col=colors)
plotPCA(rldA, intgroup=c("conditionA"))
