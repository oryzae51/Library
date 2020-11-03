##Automate count data.frame 
#Load data directory
src_dir <- c("/Users/hmkim/data/quant_data/LJK")
setwd("/Users/hmkim/data/quant_data/LJK")

#Take source file directory and save file names to src_files
src_files <- list.files(src_dir)
src_files <- src_files[!src_files %in% "summary"]
src_files

#Save source files' directory to dirVec
dirVec <- c()
for (i in src_files){
  direc <- paste(getwd(), "/", i, sep = "")
  dirVec[length(dirVec) + 1] <- direc
}
dirVec

#read.table each element of dirVec and save to tablelist as 3D list
tableList <- list()
for (i in dirVec){
  tableList[[length(tableList) + 1]] <- read.table(i, sep = "\t", header = TRUE)
}

#Make count dataframe from 3D list
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
