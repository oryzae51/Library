library(EnhancedVolcano)
MLL1_vol=read.table("/Users/hmkim/data/deg_data/MLL-KO_volcano/MLL1_volcano.csv", sep=",", header=TRUE)
MLL1_vol_sub <- subset(MLL1_vol, MLL1_vol$log2FoldChange<0)
MLL1_DNArep_sub <- subset(MLL1_vol, MLL1_vol$Genesym %in% DNArep)
MLL1_G2M_sub <- subset(MLL1_vol, MLL1_vol$Genesym %in% G2M)
MLL1_TP53TargetGenes_sub <- subset(MLL1_vol, MLL1_vol$Genesym %in% tp53tg)
#FEN1;MCM7;RFC2;POLD2;POLE3;MCM5;MCM6
DNArep <- c("FEN1", "MCM7", "RFC2", "POLD2", "POLE3", "MCM6", "MCM5")
G2M <- c("MYC", "E2F1", "CDC25A")
EnhancedVolcano(MLL1_TP53TargetGenes_sub,
                lab = MLL1_TP53TargetGenes_sub$Genesym,
                x = 'log2FoldChange',
                y = 'padj', 
                pCutoff = 0.01, 
                FCcutoff = 1, 
                title = "MLL1-KO DNA replication",
                drawConnectors = TRUE, 
                labSize = 5.0,
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
