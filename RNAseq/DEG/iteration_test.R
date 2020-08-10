src_dir <- c("/Users/hmkim/data/quant_data/LJK")
setwd("/Users/hmkim/data/quant_data/LJK")
src_files <- list.files(src_dir)
src_files
dirVec <- c()
for (i in src_files){
  direc <- paste(getwd(), "/", i, sep = "")
  dirVec[length(dirVec) + 1] <- direc
}
dirVec

tableList <- list()
for (i in dirVec){
  tableList[[length(tableList) + 1]] <- read.table(i, sep = "\t", header = TRUE)
}

for (i in 1:length(tableList)){
  tableList[[i]][, 7]
}
