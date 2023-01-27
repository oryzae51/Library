##edgeR with no replicate
##Automate count data.frame 
#Load data directory
src_dir <- c("/Users/hmkim/data/quant_data/LJK")
setwd("/Users/hmkim/data/quant_data/LJK")
src_files <- list.files(src_dir)
src_files <- src_files[!src_files %in% "summary"]
src_files
dirVec <- c()
for (i in src_files){
  direc <- paste(getwd(), "/", i, sep = "")
  dirVec[length(dirVec) + 1] <- direc
}
tableList <- list()
for (i in dirVec){
  tableList[[length(tableList) + 1]] <- read.table(i, sep = "\t", header = TRUE)
}
#count dataframe
cnts <- data.frame()
for (i in 1:length(tableList)){
  if (i == 1){
    cnts <- data.frame(c(tableList[[i]][, 7]))
  }else{
    cnts <- cbind(cnts, c(tableList[[i]][, 7]))
  }
}
rownames(cnts) <- c(tableList[[1]]$Geneid)
cnts <- cnts[rowSums(cnts > 1) >=length(cnts),]#drop genes with low counts, this is necessary for rld quality
names(cnts) <- c("S1", "WT") #change col name

##conduct edgeR
library(edgeR)
group <- factor(c(2, 1))
bcv <- 0.1
y <- DGEList(counts = cnts, group = group)
et <- exactTest(y, dispersion = bcv^2)
FDR <- p.adjust(et$table$PValue, method = "BH")
sum(FDR<0.05)
et$table <- cbind(et$table, FDR)
edge_sig <- subset(et$table, FDR<0.05)
edge_sig <- subset(edge_sig, logFC>1 | logFC< -1)
sum(edge_sig$FDR<0.01)
