#make conditionC table for downstream processing
conditionC <- c("F_ZC3H4_3D", "F_ZC3H4_3D", "Vector_3D", "Vector_3D")
sampleNameC <- c("F_ZC3H4_3D_1", "F_ZC3H4_3D_2", "Vector_3D_1", "Vector_3D_2")
sampleTableC <- data.frame(conditionC, sampleNameC, row.names="sampleNameC")
sampleTableC$conditionC <- as.factor(sampleTableC$conditionC)
all(rownames(sampleTableC) == colnames(cntsC)) #check conditionC table is correct

##Conduct DESeq2
library("DESeq2")
ddsC <- DESeqDataSetFromMatrix(countData = cntsC,
                               colData = sampleTableC,
                               design = ~ conditionC)
ddsC$conditionC <- relevel(ddsC$conditionC, ref = "Vector_3D")#Make WT sample level as reference
ddsC <- DESeq(ddsC)
ddsC$conditionC #check dds


####################################################################################
##Draw MA plots
resC <- results(ddsC, contrast = c("conditionC", "F_ZC3H4_3D", "Vector_3D"), alpha = 0.05)
resultsNames(ddsC)
resLFCC <- lfcShrink(ddsC, coef="conditionC_F_ZC3H4_3D_vs_Vector_3D", type="apeglm", res=resC)
plotMA(resLFCC, ylim=c(-5, 5))
png(paste(results_dir, "MAplot_vsRing1bKD.png", sep = ""), width = 130, height = 130, units='mm', res = 600)
plotMA(resLFC, ylim=c(-5, 5))
dev.off()

summary(resLFCC)

#########################################################################
##assign gene symbol to result####################################################
library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl_id <- c(row.names(resLFCC))
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
resLFCC$ensembl_gene_id <- c(ensembl_id)
resLFCC <- resLFCC[c(order(resLFCC$ensembl_gene_id)), ]
res_genesym <- res_genesym[c(order(res_genesym$ensembl_gene_id)), ]
resLFCC$gene_sym <- res_genesym$hgnc_symbol
resLFCC <- resLFCC[!(resLFCC$ensembl_gene_id %in% notassigned), ]

length(resLFCC$gene_sym)

#res Cutoff
summary(resLFCC)
res_padjC <- subset(resLFCC, padj<0.05)
summary(res_padjC)
res_f2_padj05C <- subset(res_padjC, abs(log2FoldChange) > 1)
summary(res_f2_padj05C)

#Raw deg
write.csv(as.data.frame(resLFCC), 
          file = paste(results_dir, "resLFCC.csv", sep = ""))
#padj<0.05, foldchange>2
write.csv(as.data.frame(res_f2_padj05C), 
          file = paste(results_dir, "resLFC_f2_padj05C.csv", sep = ""))
#foldchange>2
write.csv(as.data.frame(res_padjC), 
          file = paste(results_dir, "res_padjC.csv", sep = ""))