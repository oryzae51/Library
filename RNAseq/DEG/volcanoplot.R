library(EnhancedVolcano)
MLL1_vol=read.table("/Users/hmkim/data/deg_data/MLL-KO_volcano/MLL1_volcano.csv", sep=",", header=TRUE)
MLL2_vol=read.table("/Users/hmkim/data/deg_data/MLL-KO_volcano/MLL2_volcano.csv", sep=",", header=TRUE)
MLL4_vol=read.table("/Users/hmkim/data/deg_data/MLL-KO_volcano/MLL4_volcano.csv", sep=",", header=TRUE)
tp53tg <- read.table("/Users/hmkim/data/deg_data/MLL-KO_volcano/tp53TargetGenes.csv", sep = ",", header = FALSE)
e2ftg <- read.table("/Users/hmkim/data/deg_data/MLL-KO_volcano/E2FTargetGenes.csv", sep = ",", header = FALSE)
nfkbtg <- read.table("/Users/hmkim/data/deg_data/MLL-KO_volcano/NFKBTargetGenes.csv", sep = ",", header = FALSE)
MLL1_vol_sub <- subset(MLL1_vol, MLL1_vol$log2FoldChange<0)
MLL1_DNArep_sub <- subset(MLL1_vol, MLL1_vol$Genesym %in% DNArep)
MLL1_G2M_sub <- subset(MLL1_vol, MLL1_vol$Genesym %in% G2M)
MLL1_TP53TargetGenes_sub <- subset(MLL1_vol, MLL1_vol$Genesym %in% tp53tg$V1)
MLL2_TP53TargetGenes_sub <- subset(MLL2_vol, MLL2_vol$Genesym %in% tp53tg$V1)
MLL4_TP53TargetGenes_sub <- subset(MLL4_vol, MLL4_vol$Genesym %in% tp53tg$V1)

MLL1_E2FTargetGenes_sub <- subset(MLL1_vol, MLL1_vol$Genesym %in% e2ftg$V1)
MLL2_E2FTargetGenes_sub <- subset(MLL2_vol, MLL2_vol$Genesym %in% e2ftg$V1)
MLL4_E2FTargetGenes_sub <- subset(MLL4_vol, MLL4_vol$Genesym %in% e2ftg$V1)

MLL1_NFKBTargetGenes_sub <- subset(MLL1_vol, MLL1_vol$Genesym %in% nfkbtg$V1)
MLL2_NFKBTargetGenes_sub <- subset(MLL2_vol, MLL2_vol$Genesym %in% nfkbtg$V1)
MLL4_NFKBTargetGenes_sub <- subset(MLL4_vol, MLL4_vol$Genesym %in% nfkbtg$V1)
#FEN1;MCM7;RFC2;POLD2;POLE3;MCM5;MCM6
DNArep <- c("FEN1", "MCM7", "RFC2", "POLD2", "POLE3", "MCM6", "MCM5")
G2M <- c("MYC", "E2F1", "CDC25A")
EnhancedVolcano(MLL4_NFKBTargetGenes_sub,
                lab = MLL4_NFKBTargetGenes_sub$Genesym,
                x = 'log2FoldChange',
                y = 'padj', 
                pCutoff = 0.01, 
                FCcutoff = 1, 
                title = "MLL4-KO NFKB target gene",
                #drawConnectors = TRUE, 
                labSize = 0.0,
                boxedLabels = FALSE
                )
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
