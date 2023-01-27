setwd("/Users/kimhyoungmin/Library/Mobile\ Documents/com~apple~CloudDocs/ETL/Experiment\ data/Cowork/060922_PKH_ZWC/out")
setwd("/Users/kimhyoungmin/data/060922_PKH_ZWC/TPMcal_result/out")
setwd("/Users/kimhyoungmin/Library/Mobile\ Documents/com~apple~CloudDocs/ETL/Experiment data/Cowork/072622_HKB_DPF2mut")
setwd("/Users/kimhyoungmin/data/072622_HKB_DPF2mut/tpm")

a_1 <- read.delim("MuT-1_QCed_bam_genes.out", header = TRUE, sep = "\t")
a_2 <- read.delim("MuT-2_QCed_bam_genes.out", header = TRUE, sep = "\t")
b_1 <- read.delim("MuTDox-1_QCed_bam_genes.out", header = TRUE, sep = "\t")
b_2 <- read.delim("MuTDox-2_QCed_bam_genes.out", header = TRUE, sep = "\t")
c_1 <- read.delim("WT-1_QCed_bam_genes.out", header = TRUE, sep = "\t")
c_2 <- read.delim("WT-2_QCed_bam_genes.out", header = TRUE, sep = "\t")
d_1 <- read.delim("WTDox-1_QCed_bam_genes.out", header = TRUE, sep = "\t")
d_2 <- read.delim("WTDox-2_QCed_bam_genes.out", header = TRUE, sep = "\t")
#e_1 <- read.delim("Vector-2D-1Aligned.out.bam_genes.out", header = TRUE, sep = "\t")
#e_2 <- read.delim("Vector-2D-2Aligned.out.bam_genes.out", header = TRUE, sep = "\t")
#f_1 <- read.delim("Vector-3D-1Aligned.out.bam_genes.out", header = TRUE, sep = "\t")
#f_2 <- read.delim("Vector-3D-2Aligned.out.bam_genes.out", header = TRUE, sep = "\t")

a_1 <- a_1[, c(1, 7)]
a_2 <- a_2[, c(1, 7)]
b_1 <- b_1[, c(1, 7)]
b_2 <- b_2[, c(1, 7)]
c_1 <- c_1[, c(1, 7)]
c_2 <- c_2[, c(1, 7)]
d_1 <- d_1[, c(1, 7)]
d_2 <- d_2[, c(1, 7)]
#e_1 <- e_1[, c(1, 7)]
#e_2 <- e_2[, c(1, 7)]
#f_1 <- f_1[, c(1, 7)]
#f_2 <- f_2[, c(1, 7)]

a <- data.frame(a_1$TPM, a_2$TPM, row.names = a_1$Gene_Id)
b <- data.frame(b_1$TPM, b_2$TPM, row.names = b_1$Gene_Id)
c <- data.frame(c_1$TPM, c_2$TPM, row.names = c_1$Gene_Id)
d <- data.frame(d_1$TPM, d_2$TPM, row.names = d_1$Gene_Id)
#e <- data.frame(e_1$TPM, e_2$TPM, row.names = e_1$Gene_Id)
#f <- data.frame(f_1$TPM, f_2$TPM, row.names = f_1$Gene_Id)

TPM_df_avg <- data.frame(apply(a, 1, mean), apply(b, 1, mean), apply(c, 1, mean), apply(d, 1, mean), apply(e, 1, mean), apply(f, 1, mean), row.names = a_1$Gene_Id)
TPM_df_avg <- data.frame(apply(a, 1, mean), apply(b, 1, mean), apply(c, 1, mean), apply(d, 1, mean))
colnames(TPM_df_avg) <- c("z0_tpm_1", "z2_tpm_1", "z3_tpm_1", "v0_tpm_1", "v2_tpm_1", "v3_tpm_1")
colnames(TPM_df_avg) <- c("Mut", "MutDox", "WT", "WTDox")

TPM_df_avg <- TPM_df_avg[(rownames(TPM_df_avg) %in% resvsDOXWT_intersect$X), ]


library(ggplot2)
library(ggsignif)
TPM_df_log10 <- log10(TPM_df_avg)
TPM_df_log10.m <- reshape2::melt(TPM_df_log10, id.vars=NULL)
ggplot(TPM_df_log10.m, aes(x=variable, y=value)) +
  geom_boxplot()+
  geom_signif(comparisons = list(c("Mut", "MutDox")), 
              map_signif_level = TRUE, step_increase = 0.1, test = "t.test")+
  geom_signif(comparisons = list(c("WT", "WTDox")), 
              map_signif_level = TRUE, step_increase = 0.1, test = "t.test")+
  xlab("") + ylab("log10TPM")

ggplot(TPM_df_log10.m, aes(x=variable, y=value)) +
  geom_boxplot()+
  geom_signif(comparisons = list(c("Mut", "MutDox"), c("Mut", "WT"), c("WT", "WTDox"), c("MutDox", "WTDox")), 
              map_signif_level = TRUE, step_increase = 0.1, test = "t.test")+
  xlab("") + ylab("log10TPM")

write.table(TPM_df_log10, file = "/Users/kimhyoungmin/data/060922_PKH_ZWC/TPMcal_result/TPM_df_log10.tsv", sep="\t", col.names = TRUE, row.names = TRUE)
