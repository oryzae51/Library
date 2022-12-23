##assign gene symbol to result####################################################
library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl_id <- c(row.names(resLFC_MLL2))
ensembl_id <- substr(ensembl_id, 1, 15)
res_genesym <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id", 
                     values = ensembl_id, 
                     mart = ensembl)
#remove duplicate -> 여기 if 문 넣어서 0일때는 넘어가게 해야함
if (length(which(duplicated(res_genesym$ensembl_gene_id))) == 0){
  res_genesym <- res_genesym
}else{
  res_genesym <- res_genesym[-which(duplicated(res_genesym$ensembl_gene_id)),]
}

#add empty values to not assigned id in res_genesym
notassigned <- setdiff(ensembl_id, res_genesym$ensembl_gene_id)
empty_col <- rep(" ", time=length(notassigned))
notassigned_df <- data.frame(notassigned, empty_col)
names(notassigned_df) <- c("ensembl_gene_id", "hgnc_symbol")
res_genesym <- rbind(res_genesym, notassigned_df)
#add genesym column to res
resLFC_MLL2$ensembl_gene_id <- c(ensembl_id)
resLFC_MLL2 <- resLFC_MLL2[c(order(resLFC_MLL2$ensembl_gene_id)), ]
res_genesym <- res_genesym[c(order(res_genesym$ensembl_gene_id)), ]
resLFC_MLL2$gene_sym <- res_genesym$hgnc_symbol
resLFC_MLL2 <- resLFC_MLL2[!(resLFC_MLL2$ensembl_gene_id %in% notassigned), ]

length(resLFC_MLL2$gene_sym)