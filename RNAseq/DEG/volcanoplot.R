library(EnhancedVolcano)
MLL1_vol=read.table("/Users/hmkim/data/deg_data/MLL-KO_volcano/MLL1_volcano.csv", sep=",", header=TRUE)
MLL1_vol_sub <- subset(MLL1_vol, MLL1_vol$log2FoldChange<0)
EnhancedVolcano(MLL1_vol,
                lab = MLL1_vol$Genesym,
                x = 'log2FoldChange',
                y = 'padj', 
                pCutoff = 0.01, 
                FCcutoff = 1, 
                selectLab = c("MCM1", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7"),
                title = "MLL1-KO p53-associated gene",
                drawConnectors = TRUE, 
                labSize = 6.0,
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
