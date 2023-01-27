MuTDox_1 <- read.delim("/Users/kimhyoungmin/data/072622_HKB_DPF2mut/tpm/MuTDox-1_QCed_bam_genes.out", header = TRUE, sep = "\t")
WT_1 <- read.delim("/Users/kimhyoungmin/data/072622_HKB_DPF2mut/tpm/WT-1_QCed_bam_genes.out", header = TRUE, sep = "\t")
MuTDox_2 <- read.delim("/Users/kimhyoungmin/data/072622_HKB_DPF2mut/tpm/MuTDox-2_QCed_bam_genes.out", header = TRUE, sep = "\t")
WTDox_2 <- read.delim("/Users/kimhyoungmin/data/072622_HKB_DPF2mut/tpm/WTDox-2_QCed_bam_genes.out", header = TRUE, sep = "\t")
MuT_1 <- read.delim("/Users/kimhyoungmin/data/072622_HKB_DPF2mut/tpm/MuT-1_QCed_bam_genes.out", header = TRUE, sep = "\t")
WTDox_1 <- read.delim("/Users/kimhyoungmin/data/072622_HKB_DPF2mut/tpm/WTDox-1_QCed_bam_genes.out", header = TRUE, sep = "\t")
WT_2 <- read.delim("/Users/kimhyoungmin/data/072622_HKB_DPF2mut/tpm/WT-2_QCed_bam_genes.out", header = TRUE, sep = "\t")
MuT_2 <- read.delim("/Users/kimhyoungmin/data/072622_HKB_DPF2mut/tpm/MuT-2_QCed_bam_genes.out", header = TRUE, sep = "\t")

setwd("/Users/kimhyoungmin/data/060922_PKH_ZWC/TPMcal_result/out")

resvsDox_down_gid_sym <- read.delim("/Users/kimhyoungmin/data/072622_HKB_DPF2mut/resvsDox_down_gid_sym.txt", header = FALSE, sep = "\t")

MuTDox_1 <- MuTDox_1[(MuTDox_1$Gene_Id %in% resvsDox_down_gid_sym$V3), ]
WT_1 <- WT_1[(WT_1$Gene_Id %in% resvsDox_down_gid_sym$V3), ]
MuTDox_2 <- MuTDox_2[(MuTDox_2$Gene_Id %in% resvsDox_down_gid_sym$V3), ]
WTDox_2 <- WTDox_2[(WTDox_2$Gene_Id %in% resvsDox_down_gid_sym$V3), ]
MuT_1 <- MuT_1[(MuT_1$Gene_Id %in% resvsDox_down_gid_sym$V3), ]
WTDox_1 <- WTDox_1[(WTDox_1$Gene_Id %in% resvsDox_down_gid_sym$V3), ]
WT_2 <- WT_2[(WT_2$Gene_Id %in% resvsDox_down_gid_sym$V3), ]
MuT_2 <- MuT_2[(MuT_2$Gene_Id %in% resvsDox_down_gid_sym$V3), ]

MuTDox_1 <- MuTDox_1[,c(1, 7)]
WT_1 <- WT_1[,c(1, 7)]
MuTDox_2 <- MuTDox_2[,c(1, 7)]
WTDox_2 <- WTDox_2[,c(1, 7)]
MuT_1 <- MuT_1[,c(1, 7)]
WTDox_1 <- WTDox_1[,c(1, 7)]
WT_2 <- WT_2[,c(1, 7)]
MuT_2 <- MuT_2[,c(1, 7)]

a <- data.frame(MuT_1$TPM, MuT_2$TPM, row.names = MuT_1$Gene_Id)
b <- data.frame(MuTDox_1$TPM, MuTDox_2$TPM, row.names = MuT_1$Gene_Id)
c <- data.frame(WT_1$TPM, WT_2$TPM, row.names = MuT_1$Gene_Id)
d <- data.frame(WTDox_1$TPM, WTDox_2$TPM, row.names = MuT_1$Gene_Id)

TPM_df <- data.frame(MuT_1$TPM, MuT_2$TPM, MuTDox_1$TPM, MuTDox_2$TPM, WT_1$TPM, WT_2$TPM, WTDox_1$TPM, WTDox_2$TPM)
rownames(TPM_df) <- MuT_1$Gene_Id
write.table(TPM_df_avg, "/Users/kimhyoungmin/data/072622_HKB_DPF2mut/tpm/TPM_df_avg.tsv", sep = "\t", row.names = TRUE)
write.table(TPM_df_log10, "/Users/kimhyoungmin/data/072622_HKB_DPF2mut/tpm/TPM_df_log10.tsv", sep = "\t", row.names = TRUE)

TPM_df_avg <- data.frame(apply(a, 1, mean), apply(b, 1, mean), apply(c, 1, mean), apply(d, 1, mean))
TPM_df_avg_dox_vs_wt <- TPM_df_avg[(rownames(TPM_df_avg) %in% intersect_up), ]
colnames(TPM_df_avg) <- c("MuT", "MuTDox", "WT", "WTDox")

TPM_df_log10 <- log10(TPM_df_avg_dox_vs_wt)
TPM_df_log10.m <- reshape2::melt(TPM_df_log10, id.vars=NULL)
ggplot(TPM_df_log10.m, aes(x=variable, y=value)) +
  geom_boxplot()+
  geom_signif(comparisons = list(c("MuT", "MuTDox")), 
              map_signif_level = TRUE, step_increase = 0.1, test = "t.test")+
  geom_signif(comparisons = list(c("MuTDox", "WTDox")), 
              map_signif_level = TRUE, step_increase = 0.1, test = "t.test")+
  xlab("") + ylab("logTPM")

