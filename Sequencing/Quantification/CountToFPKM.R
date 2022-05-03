#import data as datafram(df)
##Pull down data, data preprocessing
setwd("/Users/hmkim/data/112421_FOXK1-KO_pub/quant_data/cnt")

##Automate count data.frame 
#Load data directory
src_dir <- c()
src_dir <- c("/Users/hmkim/data/112421_FOXK1-KO_pub/quant_data/cnt")
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
names(cnts) <- c("WT_1", "WT_2", "FOXK1KO_1", "FOXK1KO_2") #change col name

##Add gene length to count data
#assign annotations from Biomart
library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
#biomart needs query. So I extract gene annotation from count file in to 
#query. substr() is for remove version (.XX)
query <- c(row.names(cnts))
query <- substr(query, 1, 15)
##아래 주석처리 된 코드에서 원하는 attribute 찾을 수 있음
filters <- listFilters(ensembl)
aTtributes <- listAttributes(ensembl)
getBM_query_output <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "start_position", "end_position"),
                     filters = "ensembl_gene_id", 
                     values = query, 
                     mart = ensembl)
#add gene length column to getBM_query_output
geneLength <- abs(getBM_query_output$end-getBM_query_output$start_position)
getBM_query_output[, "geneLength"] <- geneLength

#remove duplicate
getBM_query_output <- getBM_query_output[-which(duplicated(getBM_query_output$ensembl_gene_id)),]
#add empty values to not assigned id in getBM_query_output
notassigned <- setdiff(query, getBM_query_output$ensembl_gene_id)
empty_col <- rep(" ", time=length(notassigned))
notassigned_df <- data.frame(notassigned, empty_col, empty_col, empty_col, empty_col)
names(notassigned_df) <- c("ensembl_gene_id", "hgnc_symbol", "start_position", "end_position", "geneLength")
getBM_query_output <- rbind(getBM_query_output, notassigned_df)
#add genesym, start_position, end_position, geneLength column to res
cnts$ensembl_gene_id <- c(query)
cnts <- cnts[c(order(cnts$ensembl_gene_id)), ]
getBM_query_output <- getBM_query_output[c(order(getBM_query_output$ensembl_gene_id)), ]
cnts$gene_sym <- getBM_query_output$hgnc_symbol
cnts$start_position <- getBM_query_output$start_position
cnts$end_position <- getBM_query_output$end_position
cnts$geneLength <- getBM_query_output$geneLength
cnts <- cnts[!(cnts$ensembl_gene_id %in% notassigned), ]


##countToFPKM
library(countToFPKM)
counts <- as.matrix(cnts$WT_1)
featureLength <- as.integer(cnts$geneLength)
meanLength <- as.numeric(317.532556)
fpkm_matrix <- fpkm (count, featureLength, meanLength)



file.readcounts <- system.file("extdata", "RNA-seq.read.counts.csv", package="countToFPKM")
file.annotations <- system.file("extdata", "Biomart.annotations.hg38.txt", package="countToFPKM")
file.sample.metrics <- system.file("extdata", "RNA-seq.samples.metrics.txt", package="countToFPKM")

# Import the read count matrix data into R.
counts <- as.matrix(read.csv(file.readcounts))

# Import feature annotations. 
# Assign feature lenght into a numeric vector.
gene.annotations <- read.table(file.annotations, sep="\t", header=TRUE)
featureLength <- gene.annotations$length

# Import sample metrics. 
# Assign mean fragment length into a numeric vector.
samples.metrics <- read.table(file.sample.metrics, sep="\t", header=TRUE)
meanFragmentLength <- samples.metrics$meanFragmentLength

fpkm_matrix <- fpkm (counts, featureLength, meanFragmentLength)

