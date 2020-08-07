setwd("/Users/hmkim/data/quant_data/MLL3trim_072720")
setwd("/Users/hmkim/data/quant_data/MLL3_s1_072820")
setwd("/Users/hmkim/data/quant_data/MLL-KO_073120")
MLL3_11rn=read.table("MLL1-KO-RNA1Aligned.out.sam.txt", sep="\t", header=TRUE)#for extract geneid

MLL3_11df=read.table("MLL1-KO-RNA1Aligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
MLL3_12df=read.table("MLL1-KO-RNA2Aligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
MLL3_21df=read.table("MLL2-KO-RNA2Aligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
MLL3_22df=read.table("MLL2-KO-RNA3Aligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
MLL3_31df=read.table("MLL4-KO-RNA1Aligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
MLL3_32df=read.table("MLL4-KO-RNA2Aligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
WT1df=read.table("WT-HCT116-R1Aligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")
WT2df=read.table("WT-HCT116-R2Aligned.out.sam.txt", sep="\t", header=TRUE, row.names="Geneid")


cnts<-data.frame(MLL3_11df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL_KO_star_hg38.MLL1.KO.RNA1Aligned.out.sam, 
                 MLL3_12df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL_KO_star_hg38.MLL1.KO.RNA2Aligned.out.sam, 
                 MLL3_21df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL_KO_star_hg38.MLL2.KO.RNA2Aligned.out.sam,
                 MLL3_22df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL_KO_star_hg38.MLL2.KO.RNA3Aligned.out.sam, 
                 MLL3_31df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL_KO_star_hg38.MLL4.KO.RNA1Aligned.out.sam, 
                 MLL3_32df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL_KO_star_hg38.MLL4.KO.RNA2Aligned.out.sam, 
                 WT1df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL_KO_star_hg38.WT.HCT116.R1Aligned.out.sam, 
                 WT2df$X.media.bm.790240e4.2887.451f.ad02.1b19c4b4e120.KHM.align_data.SAMfile.MLL_KO_star_hg38.WT.HCT116.R2Aligned.out.sam)
Geneid <- c(MLL3_11rn$Geneid)
rownames(cnts) <- c(MLL3_11rn$Geneid)
colnames(cnts) <- c("MLL1_1", "MLL1_2", "MLL2_1", "MLL2_2", "MLL4_1", "MLL4_2", "WT_1", "WT2")
cnts <- cnts[rowSums(cnts > 1) >=8,]


library(edgeR)
group <- factor(c(2, 2, 3, 3, 4, 4, 1, 1))
y <- DGEList(counts = cnts, group = group)
keep <- filterByExpr(y)
filtered_y <- y[keep, , keep.lib.sizes = FALSE]
norm_y <- calcNormFactors(filtered_y)
design <- model.matrix(~group)
dispers_y <- estimateDisp(norm_y, design)

#Perform quasi-likelihood F-test
fit <- glmQLFit(dispers_y, design)
qlf.2vs1 <- glmQLFTest(fit, coef = 2)
topTags(qlf.2vs1)
FDR <- p.adjust(qlf.2vs1$table$PValue, method = "BH")
sum(FDR<0.01)
qlf.2vs1$table <- cbind(qlf.2vs1$table, FDR)
edge_sig <- subset(qlf.2vs1$table, FDR<0.01)
edge_sig <- subset(qlf.2vs1$table, logFC>1 | logFC< -1)
edge_sig <- subset(edge_sig, logFC< -2)
edge_sig <- subset(edge_sig, logFC>0)
sum(edge_sig$FDR<0.01)
#edge_sig <- subset(qlf.2vs1$table, FDR < 0.01)
write.csv(as.data.frame(edge_sig), 
          file="~/data/deg_data/MLL1vsWT_FDR001.csv")

qlf.3vs1 <- glmQLFTest(fit, coef = 3)
topTags(qlf.3vs1)
FDR <- p.adjust(qlf.3vs1$table$PValue, method = "BH")
sum(FDR<0.01)
qlf.3vs1$table <- cbind(qlf.3vs1$table, FDR)
edge_sig <- subset(qlf.3vs1$table, FDR<0.01)
edge_sig <- subset(qlf.3vs1$table, logFC>1 | logFC< -1)
edge_sig <- subset(edge_sig, logFC< -2)
edge_sig <- subset(edge_sig, logFC>0)
sum(edge_sig$FDR<0.01)
#edge_sig <- subset(qlf.3vs1$table, FDR < 0.01)
write.csv(as.data.frame(edge_sig), 
          file="~/data/deg_data/MLL2vsWT_FDR001.csv")

qlf.4vs1 <- glmQLFTest(fit, coef = 4)
topTags(qlf.4vs1)
FDR <- p.adjust(qlf.4vs1$table$PValue, method = "BH")
sum(FDR<0.01)
qlf.4vs1$table <- cbind(qlf.4vs1$table, FDR)
edge_sig <- subset(qlf.4vs1$table, FDR<0.01)
edge_sig <- subset(qlf.4vs1$table, logFC>1 | logFC< -1)
edge_sig <- subset(edge_sig, logFC< -2)
edge_sig <- subset(edge_sig, logFC>0)
sum(edge_sig$FDR<0.01)
#edge_sig <- subset(qlf.4vs1$table, FDR < 0.01)
write.csv(as.data.frame(edge_sig), 
          file="~/data/deg_data/MLL4vsWT_FDR001.csv")

qlf.anova <- glmQLFTest(fit, coef = 2:4)
topTags(qlf.anova)
FDR <- p.adjust(qlf.anova$table$PValue, method = "BH")
sum(FDR<0.01)
qlf.anova$table <- cbind(qlf.anova$table, FDR)
edge_sig <- subset(qlf.anova$table, FDR<0.01)
edge_sig <- subset(edge_sig, 
                   (logFC.group2>0 & logFC.group3>0) |
                     (logFC.group3>0 & logFC.group4>0) |
                     (logFC.group2>0 & logFC.group4>0))
edge_sig <- subset(edge_sig, logFC.group2>0 & logFC.group3>0 & logFC.group4>0)
edge_sig <- subset(edge_sig, logFC>1)
sum(edge_sig$FDR<0.01)
#edge_sig <- subset(qlf.anova$table, FDR < 0.01)
write.csv(as.data.frame(edge_sig), 
          file="~/data/deg_confirm/edge_mine_up_htseq_sig_anova.csv")
write.csv(as.data.frame(qlf.anova$table), 
          file="~/data/deg_confirm/edge_htseq_sig_anova_fulldata.csv")

MLL1vsWT <- read.csv("~/data/deg_data/MLL1vsWT_FDR001.csv", sep = ",", header = TRUE, row.names = "X")
MLL2vsWT <- read.csv("~/data/deg_data/MLL2vsWT_FDR001.csv", sep = ",", header = TRUE, row.names = "X")
MLL4vsWT <- read.csv("~/data/deg_data/MLL4vsWT_FDR001.csv", sep = ",", header = TRUE, row.names = "X")

MLL1vsWT <- read.csv("~/data/deg_data/MLL1vsWT_FDR001.csv", sep = ",", header = TRUE)
MLL2vsWT <- read.csv("~/data/deg_data/MLL2vsWT_FDR001.csv", sep = ",", header = TRUE)
MLL4vsWT <- read.csv("~/data/deg_data/MLL4vsWT_FDR001.csv", sep = ",", header = TRUE)
MLLintersect <- intersect(c(MLL1vsWT$X), c(MLL2vsWT$X))
MLLintersect <- intersect(MLLintersect, c(MLL4vsWT$X))


library(venn)
library(gplots)
prot1 <- read.csv("~/data/deg_confirm/edge_mine_up.csv", header = FALSE)
prot2 <- read.csv("~/data/deg_confirm/deseq2_mine_up.csv", header = FALSE)
data <- list(A = prot1$V1, B = prot2$V1)
venn(data)

library(venn)
library(gplots)
prot2 <- read.csv("~/data/deg_confirm/edge_mine_up_s2_2_rmdot.csv", header = FALSE)
prot1 <- read.csv("~/data/deg_confirm/edge_mine_up_s2_rmdot.csv", header = FALSE)
prot3 <- read.csv("~/data/deg_confirm/edge_mine_up_s2_3_rmdot.csv", header = FALSE)
prot4 <- read.csv("~/data/deg_confirm/edge_mine_up_s2_rmdot_merge.csv", header = FALSE)

prot_ht_sig1 <- read.csv("~/data/deg_confirm/edge_mine_up_BH/edge_mine_up_htseq_sig_1_rmdot.csv", header = FALSE)
prot_ht_sig2 <- read.csv("~/data/deg_confirm/edge_mine_up_BH/edge_mine_up_htseq_sig_2_rmdot.csv", header = FALSE)
prot_ht_sig3 <- read.csv("~/data/deg_confirm/edge_mine_up_BH/edge_mine_up_htseq_sig_3_rmdot.csv", header = FALSE)
prot_ht_sig4 <- read.csv("~/data/deg_confirm/edge_mine_up_down_BH/edge_mine_up_down_htseq_sig_1_rmdot.csv", header = FALSE)
prot_ht_sig5 <- read.csv("~/data/deg_confirm/edge_mine_up_down_BH/edge_mine_up_down_htseq_sig_2_rmdot.csv", header = FALSE)
prot_ht_sig6 <- read.csv("~/data/deg_confirm/edge_mine_up_down_BH/edge_mine_up_down_htseq_sig_3_rmdot.csv", header = FALSE)
prot_ht_sig7 <- read.csv("~/data/deg_confirm/edge_mine_up_htseq_sig_1_lfd2.csv", header = FALSE)
prot_ht_sig8 <- read.csv("~/data/deg_confirm/edge_mine_up_htseq_sig_2_lfd2.csv", header = FALSE)
prot_ht_sig9 <- read.csv("~/data/deg_confirm/edge_mine_up_htseq_sig_3_lfd2.csv", header = FALSE)
prot_ht_sig10 <- read.csv("~/data/deg_confirm/edge_mine_up_htseq_sig_1_lfd-2.csv", header = FALSE)
prot_ht_sig11 <- read.csv("~/data/deg_confirm/edge_mine_up_htseq_sig_2_lfd-2.csv", header = FALSE)
prot_ht_sig12 <- read.csv("~/data/deg_confirm/edge_mine_up_htseq_sig_3_lfd-2.csv", header = FALSE)
prot_ht_sig13 <- read.csv("~/data/deg_confirm/edge_mine_up_fc_sig_1_rmdot.csv", header = FALSE)
prot_ht_sig14 <- read.csv("~/data/deg_confirm/edge_mine_up_fc_sig_2_rmdot.csv", header = FALSE)
prot_ht_sig15 <- read.csv("~/data/deg_confirm/edge_mine_up_fc_sig_3_rmdot.csv", header = FALSE)

prot_ht_sig_merge <- read.csv("~/data/deg_confirm/edge_mine_up_htseq_sig_merge.csv", header = FALSE)
prot_ht_sig_anova <- read.csv("~/data/deg_confirm/edge_mine_up_htseq_sig_anova_rmdot.csv", header = FALSE)
prot_ht_sig_intersect <- as.matrix(intersect(prot_ht_sig1$V1, prot_ht_sig2$V1))
prot_ht_sig_intersect <- as.matrix(intersect(prot_ht_sig_intersect, prot_ht_sig3$V1))
prot_ht_sig_intersect_u_d <- as.matrix(intersect(prot_ht_sig4$V1, prot_ht_sig5$V1))
prot_ht_sig_intersect_u_d <- as.matrix(intersect(prot_ht_sig_intersect_u_d, prot_ht_sig6$V1))
prot_ht_sig_intersect_d_lfd2 <- as.matrix(intersect(prot_ht_sig10$V1, prot_ht_sig11$V1))
prot_ht_sig_intersect_d_ldf2 <- as.matrix(intersect(prot_ht_sig_intersect_u_d, prot_ht_sig12$V1))
prot_ht_sig_intersect_fc <- as.matrix(intersect(prot_ht_sig13$V1, prot_ht_sig14$V1))
prot_ht_sig_intersect_fc <- as.matrix(intersect(prot_ht_sig_intersect_fc, prot_ht_sig15$V1))


public <- read.csv("~/data/deg_confirm/edge_public_up.csv", header = FALSE)
public_u_d <- read.csv("~/data/deg_confirm/edge_public_up_down.csv", header = FALSE)

#
data <- list(A = prot_ht_sig1$V1, B = prot_ht_sig2$V1)
venn(data)
data <- list(A = prot_ht_sig1$V1, C = prot_ht_sig3$V1)
venn(data)
data <- list(B = prot_ht_sig2$V1, C = prot_ht_sig3$V1)
venn(data)
data <- list(A = prot_ht_sig1$V1, B = prot_ht_sig2$V1, C = prot_ht_sig3$V1)
venn(data)
data <- list(A = prot_ht_sig7$V1, B = prot_ht_sig8$V1, C = prot_ht_sig9$V1)
venn(data)
data <- list(A = prot_ht_sig10$V1, B = prot_ht_sig11$V1, C = prot_ht_sig12$V1)
venn(data)
data <- list(A = prot_ht_sig1$V1, D = prot_ht_sig_anova$V1)
venn(data)
data <- list(A = prot_ht_sig1$V1, B = prot_ht_sig2$V1, C = prot_ht_sig3$V1, D = prot_ht_sig_anova$V1)
venn(data)


data <- list(A = public$V1, B = prot_ht_sig1$V1)
venn(data)
data <- list(A = public$V1, B = prot_ht_sig2$V1)
venn(data)
data <- list(A = public$V1, B = prot_ht_sig3$V1)
venn(data)
data <- list(A = public$V1, B = prot_ht_sig_merge$V1)
venn(data)
data <- list(A = public$V1, B = prot_ht_sig_anova$V1)
venn(data)
data <- list(A = public$V1, B = prot_ht_sig_intersect)
venn(data)
data <- list(A = public$V1, B = prot_ht_sig_intersect_fc)
venn(data)
data <- list(A = public_u_d$V1, B = prot_ht_sig_intersect_u_d)
venn(data)
data <- list(A = prot_ht_sig_intersect, B = prot_ht_sig_intersect_fc)
venn(data)
data <- list(A = public$V1, B = prot_ht_sig_intersect)
venn(data)

data <- list(KHM_pipeline_edge = prot1$V1, Public_pipeline_edge_up = prot3$V1)
venn(data)

data <- list(KHM_pipeline_edge = prot2$V1, Public_pipeline_edge_up = prot3$V1)
venn(data)

data <- list(A = prot2$V1, B = prot3$V1, C = prot1$V1)
venn(data)

data <- list(KHM_pipeline_edge = prot_ht_sig1$V1, Public_pipeline_edge_up = public$V1)
venn(data)

