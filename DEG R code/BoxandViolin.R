###Violin plot, boxplot with TPM
library('ggplot2')
library('ggsignif')
FZC3H4_0D_1 <- na.omit(read.table(paste(results_dir, "TPMcal_result/out/F-ZC3H4-0D-1Aligned.out.bam_genes.out", sep=""), header = TRUE, fill = TRUE, sep="\t"))
FZC3H4_0D_2 <- na.omit(read.table(paste(results_dir, "TPMcal_result/out/F-ZC3H4-0D-2Aligned.out.bam_genes.out", sep=""), header = TRUE, fill = TRUE, sep="\t"))
FZC3H4_2D_1 <- na.omit(read.table(paste(results_dir, "TPMcal_result/out/F-ZC3H4-2D-1Aligned.out.bam_genes.out", sep=""), header = TRUE, fill = TRUE, sep="\t"))
FZC3H4_2D_2 <- na.omit(read.table(paste(results_dir, "TPMcal_result/out/F-ZC3H4-2D-2Aligned.out.bam_genes.out", sep=""), header = TRUE, fill = TRUE, sep="\t"))
FZC3H4_3D_1 <- na.omit(read.table(paste(results_dir, "TPMcal_result/out/F-ZC3H4-3D-1Aligned.out.bam_genes.out", sep=""), header = TRUE, fill = TRUE, sep="\t"))
FZC3H4_3D_2 <- na.omit(read.table(paste(results_dir, "TPMcal_result/out/F-ZC3H4-3D-2Aligned.out.bam_genes.out", sep=""), header = TRUE, fill = TRUE, sep="\t"))
Vector_0D_1 <- na.omit(read.table(paste(results_dir, "TPMcal_result/out/Vector-0D-1Aligned.out.bam_genes.out", sep=""), header = TRUE, fill = TRUE, sep="\t"))
Vector_0D_2 <- na.omit(read.table(paste(results_dir, "TPMcal_result/out/Vector-0D-2Aligned.out.bam_genes.out", sep=""), header = TRUE, fill = TRUE, sep="\t"))
Vector_2D_1 <- na.omit(read.table(paste(results_dir, "TPMcal_result/out/Vector-2D-1Aligned.out.bam_genes.out", sep=""), header = TRUE, fill = TRUE, sep="\t"))
Vector_2D_2 <- na.omit(read.table(paste(results_dir, "TPMcal_result/out/Vector-2D-2Aligned.out.bam_genes.out", sep=""), header = TRUE, fill = TRUE, sep="\t"))
Vector_3D_1 <- na.omit(read.table(paste(results_dir, "TPMcal_result/out/Vector-3D-1Aligned.out.bam_genes.out", sep=""), header = TRUE, fill = TRUE, sep="\t"))
Vector_3D_2 <- na.omit(read.table(paste(results_dir, "TPMcal_result/out/Vector-3D-2Aligned.out.bam_genes.out", sep=""), header = TRUE, fill = TRUE, sep="\t"))

################################################
#padj>05, logfoldchange > abs1 in day 2
resLFC_2D <- na.omit(read.table(paste(results_dir, "resLFC_f2_padj05B.tsv", sep=""), header = TRUE, fill = TRUE, sep="\t"))
res2D_padj05_up <- res2D_padj[(res2D_padj$log2FoldChange >=1), ]
mat_up_2D <- betas[(rownames(betas) %in% rownames(res2D_padj05_up)), -c(1, 2, 3, 4)]
mat_up_2D.m <- reshape2::melt(mat_up_2D, id.vars=NULL)
ggplot(mat_up_2D.m, aes(x=Var2, y=value)) +
  geom_boxplot()+
  geom_signif(comparisons = list(c("treatmentF_ZC3H4.days2D", "treatmentF_ZC3H4.days3D")), 
              map_signif_level = FALSE, step_increase = 0.1, test = "t.test") +
  geom_jitter(shape=16, position=position_jitter(0.1)) + xlab("vs days0D") + ylab("coefficients")

write.csv(mat_up_2D.df, 
          file = paste(results_dir, "TP_f2_padj05_up_coef.csv", sep = ""), quote = FALSE)

ggplot(mat_up_2D.m, aes(x=Var2, y=value, fill=Var2)) +
  geom_violin()+
  geom_signif(comparisons = list(c("treatmentF_ZC3H4.days2D", "treatmentF_ZC3H4.days3D")), 
              map_signif_level = FALSE, step_increase = 0.1, test = "t.test")
  

res2D_padj05_down <- res2D_padj[(res2D_padj$log2FoldChange <=-1), ]
mat_down_2D <- betas[(rownames(betas) %in% rownames(res2D_padj05_down)), -c(1, 2, 3, 4)]
mat_down_2D.m <- reshape2::melt(mat_down_2D, id.vars=NULL)
ggplot(mat_down_2D.m, aes(x=Var2, y=value)) +
  geom_boxplot()+
  geom_signif(comparisons = list(c("treatmentF_ZC3H4.days2D", "treatmentF_ZC3H4.days3D")), 
              map_signif_level = FALSE, step_increase = 0.1, test = "t.test") + 
  geom_jitter(shape=16, position=position_jitter(0.1)) + xlab("vs days0D") + ylab("coefficients")

write.csv(mat_down_2D.df, 
          file = paste(results_dir, "TP_f2_padj05_down_coef.csv", sep = ""), quote = FALSE)
####################################################3
#padj>05, logfoldchange > abs1 in day 3
res3D_padj05_up <- res3D_padj[(res3D_padj$log2FoldChange >=1), ]
mat_up_3D <- betas[(rownames(betas) %in% rownames(res3D_padj05_up)), -c(1, 2, 3, 4)]
mat_up_3D.m <- reshape2::melt(mat_up_3D, id.vars=NULL)
ggplot(mat_up_3D.m, aes(x=Var2, y=value)) +
  geom_boxplot()+
  geom_signif(comparisons = list(c("treatmentF_ZC3H4.days2D", "treatmentF_ZC3H4.days3D")), 
              map_signif_level = FALSE, step_increase = 0.1, test = "t.test") +
  geom_jitter(shape=16, position=position_jitter(0.1)) + xlab("vs days0D") + ylab("coefficients")
write.csv(mat_up_3D, 
          file = paste(results_dir, "TP_3D_f2_padj05_up_coef_n1034.csv", sep = ""), quote = FALSE)

res3D_padj05_down <- res3D_padj[(res3D_padj$log2FoldChange <=-1), ]
mat_down_3D <- betas[(rownames(betas) %in% rownames(res3D_padj05_down)), -c(1, 2, 3, 4)]
mat_down_3D.m <- reshape2::melt(mat_down_3D, id.vars=NULL)
ggplot(mat_down_3D.m, aes(x=Var2, y=value)) +
  geom_boxplot()+
  geom_signif(comparisons = list(c("treatmentF_ZC3H4.days2D", "treatmentF_ZC3H4.days3D")), 
              map_signif_level = FALSE, step_increase = 0.1, test = "t.test") +
  geom_jitter(shape=16, position=position_jitter(0.1)) + xlab("vs days0D") + ylab("coefficients")
write.csv(mat_down_3D, 
          file = paste(results_dir, "TP_3D_f2_padj05_down_coef_n1490.csv", sep = ""), quote = FALSE)

#rowname(FZC3H4_0D_tpm)
#FZC3H4_0D_tpm <- c((FZC3H4_0D_1$TPM+FZC3H4_0D_2$TPM)/2)
#FZC3H4_2D_tpm <- c((FZC3H4_2D_1$TPM+FZC3H4_2D_2$TPM)/2)
#FZC3H4_3D_tpm <- c((FZC3H4_3D_1$TPM+FZC3H4_3D_2$TPM)/2)
#Vector_0D_tpm <- c((Vector_0D_1$TPM+FZC3H4_0D_2$TPM)/2)
#Vector_2D_tpm <- c((Vector_2D_1$TPM+FZC3H4_2D_2$TPM)/2)
#Vector_3D_tpm <- c((Vector_3D_1$TPM+FZC3H4_3D_2$TPM)/2)

#tpm_merge <- data.frame(Vector_0D_tpm, Vector_2D_tpm, Vector_3D_tpm, FZC3H4_0D_tpm, FZC3H4_2D_tpm, FZC3H4_3D_tpm, row.names = (FZC3H4_0D_1$Gene_Id))
#tpm_merge <- tpm_merge[rowSums(tpm_merge > 1) >=length(tpm_merge),]
#tpm_merge <- tpm_merge[rownames(tpm_merge) %in% rownames(mat_down), ]
#tpm_merge <- tpm_merge[!(rownames(tpm_merge)=="ENSG00000278771.1"),]
#tpm_merge.m <- reshape2::melt(tpm_merge, id.vars=NULL)

mat_up.m <- reshape2::melt(mat_up, id.vars=NULL)
mat_down.m <- reshape2::melt(mat_down, id.vars=NULL)

ggplot(tpm_merge.m, aes(x=variable, y=value)) +
  geom_flat_violin()

ggplot(tpm_merge.m, aes(x=variable, y=value)) +
  geom_boxplot()+
  geom_signif()

ggplot(mat_up.m, aes(x=Var2, y=value)) +
  geom_boxplot()+
  geom_signif(comparisons = list(c("treatmentF_ZC3H4.days2D", "treatmentF_ZC3H4.days3D")), 
              map_signif_level = FALSE, step_increase = 0.1, test = "t.test") + 
  geom_jitter(shape=16, position=position_jitter(0.1)) + xlab("vs Vector") + ylab("coefficients")
ggplot(mat_down.m, aes(x=Var2, y=value)) +
  geom_boxplot()+
  geom_signif(comparisons = list(c("treatmentF_ZC3H4.days2D", "treatmentF_ZC3H4.days3D")), 
                             map_signif_level = FALSE, step_increase = 0.1, test = "t.test") + 
  geom_jitter(shape=16, position=position_jitter(0.1)) + xlab("vs Vector") + ylab("coefficients")
  

testdf = data.frame(AAA=rnorm(100,1,1),BBB=rnorm(100,2,1.5),CCC=rnorm(100,1.5,1.2))
df.m <- reshape2::melt(testdf, id.vars = NULL)
ggplot(df.m, aes(x=variable, y=value)) + geom_violin()
