library(EnhancedVolcano)
MLL1_vol=read.table("/Users/hmkim/data/deg_data/MLL-KO_volcano/MLL1_volcano.csv", sep=",", header=TRUE)
MLL2_vol=read.table("/Users/hmkim/data/deg_data/MLL-KO_volcano/MLL2_volcano.csv", sep=",", header=TRUE)
MLL4_vol=read.table("/Users/hmkim/data/deg_data/MLL-KO_volcano/MLL4_volcano.csv", sep=",", header=TRUE)
LIN_vol=read.table("/Users/hmkim/data/deg_data/MLL-KO_volcano/190524_FLAG-KLF2_results.csv", sep=",", header=TRUE)
#tp53tg <- read.table("/Users/hmkim/data/deg_data/target_gene_list/tp53TargetGenes.csv", sep = ",", header = FALSE)
e2ftg <- read.table("/Users/hmkim/data/deg_data/target_gene_list/E2F.csv", sep = ",", header = FALSE)
myctg <- read.table("/Users/hmkim/data/deg_data/target_gene_list/MYC.csv", sep = ",", header = FALSE)
foxa1tg <- read.table("/Users/hmkim/data/deg_data/target_gene_list/FOXA1.csv", sep = ",", header = FALSE)
gata1tg <- read.table("/Users/hmkim/data/deg_data/target_gene_list/GATA1.csv", sep = ",", header = FALSE)
nanogtg <- read.table("/Users/hmkim/data/deg_data/target_gene_list/NANOG.csv", sep = ",", header = FALSE)

MLL1_vol_sub <- subset(MLL1_vol, MLL1_vol$log2FoldChange<0)
MLL1_DNArep_sub <- subset(MLL1_vol, MLL1_vol$Genesym %in% DNArep)
MLL1_G2M_sub <- subset(MLL1_vol, MLL1_vol$Genesym %in% G2M)
MLL1_TP53TargetGenes_sub <- subset(MLL1_vol, MLL1_vol$Genesym %in% tp53tg$V1)
MLL2_TP53TargetGenes_sub <- subset(MLL2_vol, MLL2_vol$Genesym %in% tp53tg$V1)
MLL4_TP53TargetGenes_sub <- subset(MLL4_vol, MLL4_vol$Genesym %in% tp53tg$V1)

MLL1_E2FTargetGenes_sub <- subset(MLL1_vol, MLL1_vol$Genesym %in% e2ftg$V1)
MLL2_E2FTargetGenes_sub <- subset(MLL2_vol, MLL2_vol$Genesym %in% e2ftg$V1)
MLL4_E2FTargetGenes_sub <- subset(MLL4_vol, MLL4_vol$Genesym %in% e2ftg$V1)
write.csv(as.data.frame(MLL1_E2FTargetGenes_sub), file = "~/data/deg_data/MLL-KO_volcano/MLL1_E2FTargetGenes_sub.csv")
write.csv(as.data.frame(MLL2_E2FTargetGenes_sub), file = "~/data/deg_data/MLL-KO_volcano/MLL2_E2FTargetGenes_sub.csv")
write.csv(as.data.frame(MLL4_E2FTargetGenes_sub), file = "~/data/deg_data/MLL-KO_volcano/MLL4_E2FTargetGenes_sub.csv")


MLL1_NFKBTargetGenes_sub <- subset(MLL1_vol, MLL1_vol$Genesym %in% nfkbtg$V1)
MLL2_NFKBTargetGenes_sub <- subset(MLL2_vol, MLL2_vol$Genesym %in% nfkbtg$V1)
MLL4_NFKBTargetGenes_sub <- subset(MLL4_vol, MLL4_vol$Genesym %in% nfkbtg$V1)

MLL1_MYCTargetGenes_sub <- subset(MLL1_vol, MLL1_vol$Genesym %in% myctg$V1)
MLL2_MYCTargetGenes_sub <- subset(MLL2_vol, MLL2_vol$Genesym %in% myctg$V1)
MLL4_MYCTargetGenes_sub <- subset(MLL4_vol, MLL4_vol$Genesym %in% myctg$V1)
write.csv(as.data.frame(MLL1_MYCTargetGenes_sub), file = "~/data/deg_data/MLL-KO_volcano/MLL1_MYCTargetGenes_sub.csv")
write.csv(as.data.frame(MLL2_MYCTargetGenes_sub), file = "~/data/deg_data/MLL-KO_volcano/MLL2_MYCTargetGenes_sub.csv")
write.csv(as.data.frame(MLL4_MYCTargetGenes_sub), file = "~/data/deg_data/MLL-KO_volcano/MLL4_MYCTargetGenes_sub.csv")

MLL1_FOXA1TargetGenes_sub <- subset(MLL1_vol, MLL1_vol$Genesym %in% foxa1tg$V1)
MLL2_FOXA1TargetGenes_sub <- subset(MLL2_vol, MLL2_vol$Genesym %in% foxa1tg$V1)
MLL4_FOXA1TargetGenes_sub <- subset(MLL4_vol, MLL4_vol$Genesym %in% foxa1tg$V1)
write.csv(as.data.frame(MLL1_FOXA1TargetGenes_sub), file = "~/data/deg_data/MLL-KO_volcano/MLL1_FOXA1TargetGenes_sub.csv")
write.csv(as.data.frame(MLL2_FOXA1TargetGenes_sub), file = "~/data/deg_data/MLL-KO_volcano/MLL2_FOXA1TargetGenes_sub.csv")
write.csv(as.data.frame(MLL4_FOXA1TargetGenes_sub), file = "~/data/deg_data/MLL-KO_volcano/MLL4_FOXA1TargetGenes_sub.csv")

MLL1_GATA1TargetGenes_sub <- subset(MLL1_vol, MLL1_vol$Genesym %in% gata1tg$V1)
MLL2_GATA1TargetGenes_sub <- subset(MLL2_vol, MLL2_vol$Genesym %in% gata1tg$V1)
MLL4_GATA1TargetGenes_sub <- subset(MLL4_vol, MLL4_vol$Genesym %in% gata1tg$V1)
write.csv(as.data.frame(MLL1_GATA1TargetGenes_sub), file = "~/data/deg_data/MLL-KO_volcano/MLL1_GATA1TargetGenes_sub.csv")
write.csv(as.data.frame(MLL2_GATA1TargetGenes_sub), file = "~/data/deg_data/MLL-KO_volcano/MLL2_GATA1TargetGenes_sub.csv")
write.csv(as.data.frame(MLL4_GATA1TargetGenes_sub), file = "~/data/deg_data/MLL-KO_volcano/MLL4_GATA1TargetGenes_sub.csv")

MLL1_NANOGTargetGenes_sub <- subset(MLL1_vol, MLL1_vol$Genesym %in% nanogtg$V1)
MLL2_NANOGTargetGenes_sub <- subset(MLL2_vol, MLL2_vol$Genesym %in% nanogtg$V1)
MLL4_NANOGTargetGenes_sub <- subset(MLL4_vol, MLL4_vol$Genesym %in% nanogtg$V1)
write.csv(as.data.frame(MLL1_NANOGTargetGenes_sub), file = "~/data/deg_data/MLL-KO_volcano/MLL1_NANOGTargetGenes_sub.csv")
write.csv(as.data.frame(MLL2_NANOGTargetGenes_sub), file = "~/data/deg_data/MLL-KO_volcano/MLL2_NANOGTargetGenes_sub.csv")
write.csv(as.data.frame(MLL4_NANOGTargetGenes_sub), file = "~/data/deg_data/MLL-KO_volcano/MLL4_NANOGTargetGenes_sub.csv")


#FEN1;MCM7;RFC2;POLD2;POLE3;MCM5;MCM6
DNArep <- c("FEN1", "MCM7", "RFC2", "POLD2", "POLE3", "MCM6", "MCM5")
G2M <- c("MYC", "E2F1", "CDC25A")
EnhancedVolcano(MLL4_MYCTargetGenes_sub,
                lab = MLL4_MYCTargetGenes_sub$Genesym,
                x = 'log2FoldChange',
                y = 'padj', 
                pCutoff = 0.01, 
                FCcutoff = 1, 
                title = "MLL4-KO MYC target gene",
                #drawConnectors = TRUE, 
                labSize = 0.0,
                boxedLabels = FALSE
                )

EnhancedVolcano(LIN_vol,
                lab = LIN_vol$Gene.names,
                x = 'Log2.Difference.Klf2_control',
                y = 'X.Log.p.value', 
                pCutoff = 0.01, 
                FCcutoff = 1, 
                title = "KLF2 MASS",
                #drawConnectors = TRUE, 
                labSize = 5.0,
                boxedLabels = FALSE
)
LIN_vol$X.Log.p.value
LIN_vol$Log2.Difference.Klf2_control

summary(resLFC_MLL1)
write.csv(as.data.frame(resLFC_MLL1), file = "~/data/asdfasdf.csv")

library(EnhancedVolcano)
EnhancedVolcano(resLFC_MLL1,
                lab = resLFC_MLL1$Genesym,
                x = 'log2FoldChange',
                y = 'padj', 
                pCutoff = 0.01, 
                FCcutoff = 1, 
                selectLab = c("TP53", "TP73", "NFKB1", "BCL2"), 
                drawConnectors = TRUE)
