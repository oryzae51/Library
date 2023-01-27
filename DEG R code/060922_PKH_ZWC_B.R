#make conditionB table for downstream processing
conditionB <- c("F_ZC3H4_2D", "F_ZC3H4_2D", "Vector_2D", "Vector_2D")
sampleNameB <- c("F_ZC3H4_2D_1", "F_ZC3H4_2D_2", "Vector_2D_1", "Vector_2D_2")
sampleTableB <- data.frame(conditionB, sampleNameB, row.names="sampleNameB")
sampleTableB$conditionB <- as.factor(sampleTableB$conditionB)
all(rownames(sampleTableB) == colnames(cntsB)) #check conditionB table is correct

##Conduct DESeq2
library("DESeq2")
ddsB <- DESeqDataSetFromMatrix(countData = cntsB,
                               colData = sampleTableB,
                               design = ~ conditionB)
ddsB$conditionB <- relevel(ddsB$conditionB, ref = "Vector_2D")#Make WT sample level as reference
ddsB <- DESeq(ddsB)
ddsB$conditionB #check dds


####################################################################################
##Draw MA plots
resB <- results(ddsB, contrast = c("conditionB", "F_ZC3H4_2D", "Vector_2D"), alpha = 0.05)
resultsNames(ddsB)
resLFCB <- lfcShrink(ddsB, coef="conditionB_F_ZC3H4_2D_vs_Vector_2D", type="apeglm", res=resB)
plotMA(resLFCB, ylim=c(-5, 5))
png(paste(results_dir, "MAplot_vsRing1bKD.png", sep = ""), width = 130, height = 130, units='mm', res = 600)
plotMA(resLFCB, ylim=c(-5, 5))
dev.off()

summary(resLFCB)

#########################################################################
##assign gene symbol to result####################################################
library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl_id <- c(row.names(resLFCB))
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
resLFCB$ensembl_gene_id <- c(ensembl_id)
resLFCB <- resLFCB[c(order(resLFCB$ensembl_gene_id)), ]
res_genesym <- res_genesym[c(order(res_genesym$ensembl_gene_id)), ]
resLFCB$gene_sym <- res_genesym$hgnc_symbol
resLFCB <- resLFCB[!(resLFCB$ensembl_gene_id %in% notassigned), ]

length(resLFCB$gene_sym)

#res Cutoff
summary(resLFCB)
res_padjB <- subset(resLFCB, padj<0.05)
summary(res_padjB)
res_f2_padj05B <- subset(res_padjB, abs(log2FoldChange) > 1)
summary(res_f2_padj05B)

#Raw deg
write.csv(as.data.frame(resLFCB), 
          file = paste(results_dir, "resLFCB.csv", sep = ""))
#padj<0.05, foldchange>2
write.csv(as.data.frame(res_f2_padj05B), 
          file = paste(results_dir, "resLFC_f2_padj05B.csv", sep = ""))
#foldchange>2
write.csv(as.data.frame(res_padjB), 
          file = paste(results_dir, "res_padjB.csv", sep = ""))

###############################################################################

##Count data transformation/normalization for downstream analyze
#Normalizae with rlog in DESeq2
library("vsn")
rldB <- rlog(ddsB, blind=FALSE)
meanSdPlot(assay(rldB))
#Compare with normTrnaform. You don't need to do when you are making data. 
ntdB <- normTransform(ddsB)
meanSdPlot(assay(ntdB))


##PCB
#sampleDistsB1 <- dist(t(assay(rldB)))
sampleDistsB <- dist(scale(t(assay(rldB))))#recommand this for consistancy
library("RColorBrewer")
sampleDistMatrixB <- as.matrix(sampleDistsB)
rownames(sampleDistMatrixB) <- paste(rldB$conditionB, sep="-")
#colnames(sampleDistMatrixB) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrixB,
         clustering_distance_rows=sampleDistsB,
         clustering_distance_cols=sampleDistsB,
         col=colors)
plotPCA(rldB, intgroup=c("conditionB"))