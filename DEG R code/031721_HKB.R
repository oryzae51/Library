#import data as datafram(df)
##Pull down data, data preprocessing
setwd("/Users/hmkim/data/quant_data/HKB/quant")

##Automated count data.frame 
#Load data directory
src_dir <- c()
src_dir <- c("/Users/hmkim/data/quant_data/HKB/quant")
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
names(cnts) <- c("DPF2_1", "DPF2_2", "WT_1", "WT_2") #change col name

#make condition table for downstream processing
condition <- c("DPF2", "DPF2", "WT", "WT")
sampleName <- c("DPF2_1", "DPF2_2", "WT_1", "WT_2")
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

##lfc
#DPF2-KO
res <- results(dds, contrast = c("condition", "DPF2", "WT"), alpha = 0.05)
resultsNames(dds)
res <- lfcShrink(dds, coef="condition_DPF2_vs_WT", type="apeglm", res=res)

##assign gene symbol to result
library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl_id <- c(row.names(res))
ensembl_id <- substr(ensembl_id, 1, 15)
res_genesym <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id", 
                     values = ensembl_id, 
                     mart = ensembl)
#remove duplicate
res_genesym <- res_genesym[-which(duplicated(res_genesym$ensembl_gene_id)),]
#add empty values to not assigned id in res_genesym
notassigned <- setdiff(ensembl_id, res_genesym$ensembl_gene_id)
empty_col <- rep(" ", time=length(notassigned))
notassigned_df <- data.frame(notassigned, empty_col)
names(notassigned_df) <- c("ensembl_gene_id", "hgnc_symbol")
res_genesym <- rbind(res_genesym, notassigned_df)
#add genesym column to res
res$ensembl_gene_id <- c(ensembl_id)
res <- res[c(order(res$ensembl_gene_id)), ]
res_genesym <- res_genesym[c(order(res_genesym$ensembl_gene_id)), ]
res$gene_sym <- res_genesym$hgnc_symbol
res <- res[!(res$ensembl_gene_id %in% notassigned), ]

#res Cutoff
summary(res)
res_padj <- subset(res, padj<0.05)
summary(res_padj)
res_f2_padj05 <- subset(res_padj, abs(log2FoldChange) > 1)
summary(res_f2_padj05)

##Save deg_data
#Raw deg
write.csv(as.data.frame(res), 
          file = "/Users/hmkim/data/Cowork/HKB/deg_data/res.csv")
#padj<0.05, foldchange>2
write.csv(as.data.frame(res_f2_padj05), 
          file = "/Users/hmkim/data/Cowork/HKB/deg_data/res_f2_padj05.csv")
###############################################################################

###############################################################################
#Draw MA plot
png("/Users/hmkim/data/Cowork/HKB/033121_MAplot_padj05.png", 
    width = 130, 
    height = 130, 
    units='mm', 
    res = 300)
plotMA(res, ylim=c(-6, 6))
dev.off()
###############################################################################

###############################################################################
##Volcano plot
library(EnhancedVolcano)
res_df <- data.frame(res)
#tp53 target gene labeling
tp53targetgenes <- read.csv("~/data/tp53TargetGenes.csv", header = FALSE, sep = ",")
tp53targetgenes <- tp53targetgenes$V1
png("/Users/hmkim/data/Cowork/HKB/040021_Volcanoplot.png", 
    width = 260, 
    height = 260, 
    units='mm', 
    res = 300)
EnhancedVolcano(res_df,
                lab = res_df$gene_sym,
                x = 'log2FoldChange',
                y = 'padj', 
                pCutoff = 0.05, 
                FCcutoff = 1, 
                title = "DPF2 KO",
                #drawConnectors = TRUE, 
                labSize = 2,
                pointSize = 1,
                boxedLabels = FALSE, 
                legendLabSize = 10,
                legendIconSize = 3.0
)
dev.off()
###############################################################################
library(gplots)
data <- list("Public_up" = up_public, "Pipeline_up" = up_inter_sym)
venn(data)
data <- list("Public_down" = down_public, "Pipeline_down" = down_inter_sym)
venn(data)
data <- list("Public_up" = up_public, "Pipeline_down" = down_inter_sym)
venn(data)
data <- list("Public_down" = down_public, "Pipeline_up" = up_inter_sym)
venn(data)

##Count data transformation/normalization for downstream analyze
#Normalizae with rlog in DESeq2
library("vsn")
rld <- rlog(dds, blind=FALSE)
meanSdPlot(assay(rld))

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
png("/Users/hmkim/data/Cowork/HKB/033121_PCA.png", width = 130, height = 130, units='mm', res = 300)
plotPCA(rld, intgroup=c("condition"))
dev.off()

################################################################################
##hirachy heatmap
#heatmap with top 1000
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
rld_f2_padj05 <- rld_f2_padj05[rownames(rld_f2_padj05) %in% rownames(res_f2_padj05), ]
#cluster gene and sample with pearson correlation
library("pheatmap")
# Pairwise correlation between samples (columns)
cols.cor <- cor(rld_f2_padj05, use = "pairwise.complete.obs", method = "pearson")
# Pairwise correlation between rows (genes)
rows.cor <- cor(t(rld_f2_padj05), use = "pairwise.complete.obs", method = "pearson")
#Draw heatmap
png("/Users/hmkim/data/Cowork/HKB/033121_ClusteringHeatmap_f2_padj05.png", 
    width = 130, 
    height = 130, 
    units='mm', 
    res = 300)
hm <- pheatmap(rld_f2_padj05, 
               scale = "row", 
               clustering_distance_cols = as.dist(1 - cols.cor), 
               clustering_distance_rows = as.dist(1 - rows.cor), 
               col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), 
               show_rownames = FALSE
)
dev.off()
###############################################################################

#tree cutting and cluster assigning
hcluster <- as.hclust( hm$rowDendrogram )
cut_hcluster<-cutree( hcluster, h=10 )
factor(cut_hcluster)
cut_hcluster<-as.data.frame(cut_hcluster)
cluster_geneid <- row.names(cut_hcluster)
cluster_geneid <- substr(cluster_geneid, 1, 15)
cut_hcluster <- data.frame(cut_hcluster, cluster_geneid, hm$rowInd)
colnames(cut_hcluster) <- c("cluster", "gene_id", "rowInd")
write.csv(as.data.frame(cut_hcluster), 
          file = "/Users/hmkim/data/Cowork/HKB/deg_data/cut_hcluster.csv")


write.csv(assay(rld)[topVarGenes,],file="~/data/testMLL3trim.csv")
write.csv(topVarGenes,file="~/data/testMLL3trim.csv")
write.csv(assay(rld),file="~/data/testMLL3trim.csv")
geneex <- assay(rld)[topVarGenes, ]

macro_up <- read.csv("~/data/macro_up.csv", header = FALSE, sep = ",")
macro_up <- macro_up$V1
macro_down <- read.csv("~/data/macro_down.csv", header = FALSE, sep = ",")
macro_down <- macro_down$V1

me_up <- read.csv("/Users/hmkim/data/Cowork/HKB/deg_data/DPF2_up_f2_padj05.csv", header = FALSE, sep = ",")
me_up <- me_up$V1
me_down <- read.csv("/Users/hmkim/data/Cowork/HKB/deg_data/DPF2_down_f2_padj05.csv", header = FALSE, sep = ",")
me_down <- me_down$V1

library(venn)
library(gplots)
data <- list(Macrogen_up = macro_up, HKB_up = me_up)
venn(data)
data <- list(Macrogen_down = macro_down, HKB_down = me_down)
venn(data)