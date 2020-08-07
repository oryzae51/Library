###Ensembl -> Gene symbol converter
#and add gene symbel to res(DESeq2)
library("biomaRt")
Geneid_con <- c(rownames(resLFC_MLL1))
Geneid_con_split <- strsplit(Geneid_con, split = "[.]")
Geneid_rmdot <- unlist(Geneid_con_split)[2*(1:length(Geneid_con))-1]
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters = "ensembl_gene_id", attributes=c("ensembl_gene_id", "hgnc_symbol"), values = Geneid_rmdot, mart=mart)
G_list <- G_list[c(order(G_list$ensembl_gene_id)),]
resLFC_MLL1 <- resLFC_MLL1[c(order(rownames(resLFC_MLL1))),]
resLFC_MLL1$Genesym <- G_list$hgnc_symbol

Geneid_con <- c(rownames(resLFC_MLL2))
Geneid_con_split <- strsplit(Geneid_con, split = "[.]")
Geneid_rmdot <- unlist(Geneid_con_split)[2*(1:length(Geneid_con))-1]
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters = "ensembl_gene_id", attributes=c("ensembl_gene_id", "hgnc_symbol"), values = Geneid_rmdot, mart=mart)
G_list <- G_list[c(order(G_list$ensembl_gene_id)),]
resLFC_MLL2 <- resLFC_MLL2[c(order(rownames(resLFC_MLL2))),]
resLFC_MLL2$Genesym <- G_list$hgnc_symbol

Geneid_con <- c(rownames(resLFC_MLL4))
Geneid_con_split <- strsplit(Geneid_con, split = "[.]")
Geneid_rmdot <- unlist(Geneid_con_split)[2*(1:length(Geneid_con))-1]
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters = "ensembl_gene_id", attributes=c("ensembl_gene_id", "hgnc_symbol"), values = Geneid_rmdot, mart=mart)
G_list <- G_list[c(order(G_list$ensembl_gene_id)),]
resLFC_MLL4 <- resLFC_MLL4[c(order(rownames(resLFC_MLL4))),]
resLFC_MLL4$Genesym <- G_list$hgnc_symbol



library("biomaRt")
Geneid_rmdot <- read.table("~/data/deg_data/MLL-KO_diff/MLL4_KO_diff_down2_rm_dot.csv", header = TRUE)
Geneid_con_split <- strsplit(Geneid_con, split = "[.]")
Geneid_rmdot <- unlist(Geneid_con_split)[2*(1:length(Geneid_con))-1]
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters = "ensembl_gene_id", attributes=c("ensembl_gene_id", "hgnc_symbol"), values = Geneid_rmdot, mart=mart)
G_list <- G_list[c(order(G_list$ensembl_gene_id)),]
resLFC_MLL1 <- resLFC_MLL1[c(order(rownames(resLFC_MLL1))),]
resLFC_MLL1$Genesym <- G_list$hgnc_symbol

write.csv(as.data.frame(G_list), file = "/Users/hmkim/data/dasdfasdf.csv")
write.csv(as.data.frame(resLFC_MLL4), 
          file = "~/data/deg_data/MLL-KO_down/MLL4_down.csv")



