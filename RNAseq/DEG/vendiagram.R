setwd("/Users/hmkim/data/buffer/MYC_KD")

# Load library
library(VennDiagram)

# Generate 3 sets of 200 words
set1 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set2 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set3 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")

set1 <- read.csv("EP400NL_up.txt", header = FALSE)
set2 <- read.csv("MYC_down.txt", header = FALSE)
set1 <- read.csv("/Users/hmkim/data/buffer/EP400NL/EP400NL RNAseq data/RNA-seq raw data/csv and txt files/EP400NL_f1-5_up_list.txt", header = FALSE)
set2 <- read.csv("/Users/hmkim/data/buffer/MYC_KD/MYC_down_f1.5_padj05.txt", header = FALSE)
set1 <- read.csv("/Users/hmkim/data/buffer/MYC_KD/EP400NL_up_f1.3.txt", header = FALSE)
set2 <- res_f1.3_padj05$gene_sym
  

# Chart
venn.diagram(
  x = list(set1$V1, set2),
  category.names = c("EP400NL induction\nUp-regulated genes", "MYC conditional KD\nDown-regulated gene"),
  filename = '#16_venn_diagramm.png',
  height = 480 ,
  width = 480 ,
  resolution = 300,
  output=TRUE, 
  #col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c("#440154ff", '#21908dff'),
  cex = 0.6,
  fontfamily = "sans",
  cat.cex = 0.35,
  cat.default.pos = "outer",
  cat.pos = c(+0, 0),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff'),
  lty = 'blank'
)

z<- c(set1$V1, set2$V1)
z<- intersect(set1$V1, set2$V1)
write.csv(z,file="EP_MYC_intersect.csv")
