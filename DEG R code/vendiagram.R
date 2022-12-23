results_dir <- "/Users/kimhyoungmin/Library/Mobile Documents/com~apple~CloudDocs/ETL/Cowork/060922_PKH_ZWC/"

setwd("/Users/kimhyoungmin/Library/Mobile Documents/com~apple~CloudDocs/ETL/Experiment\ data/2022\ data/082222_KHM_HeLaS3_KO_results/")
# Load library
library(VennDiagram)

#Load files
set1 <- read.csv("res_f15_padj05_FK1.csv", header = TRUE)
set2 <- read.csv("res_f15_padj05_FK2.csv", header = TRUE)
set3 <- read.csv("res_f15_padj05_FKD.csv", header = TRUE)
set4 <- read.csv("res_f15_padj05_MLL2.csv", header = TRUE)

set1 <- read.csv("res_f2_padj05_FK1.csv", header = TRUE)
set2 <- read.csv("res_f2_padj05_FK2.csv", header = TRUE)
set3 <- read.csv("res_f2_padj05_FKD.csv", header = TRUE)
set4 <- read.csv("res_f2_padj05_MLL2.csv", header = TRUE)


#Extract foldchange lager than 1 or less than -1
set1_down <- set1[set1$log2FoldChange>=0, ]
set2_down <- set2[set2$log2FoldChange>=0, ]
set3_down <- set3[set3$log2FoldChange>=0, ]
set4_down <- set4[set4$log2FoldChange>=0, ]


#Two categories
venn.diagram(
  x = list(set1_down$ensembl_gene_id, set2_down$ensembl_gene_id, set3_down$ensembl_gene_id, set4_down$ensembl_gene_id),
  category.names = c("FOXK1-KO", "FOXK2-KO", "FOXK-DKO", "MLL2-KO"),
  filename = 'Venn_up_f15.png',
  height = 480 ,
  width = 480 ,
  resolution = 300,  
  output=TRUE, 
  #col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c("#440154ff", '#21908dff', '#3B528BFF', '#FDE725FF'),
  cex = 0.6,
  fontfamily = "sans",
  cat.cex = 0.35,
  cat.default.pos = "outer",
  cat.pos = c(+0, 0, 0, 0), #sample 갯수에 따라 조절
  cat.dist = c(0.055, 0.055, 0.055, 0.055), #sample 갯수에 따라 조절
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#3B528BFF', '#FDE725FF'),
  lty = 'blank'
)

# Chart - Three categories
venn.diagram(
  x = list(set1_down$ensembl_gene_id, set2_down$ensembl_gene_id, set3_down$ensembl_gene_id),
  category.names = c("FOXK1-KO", "FOXK2-KO", "FOXK-DKO"),
  filename = 'Venn_down_FK1_FK2_FKD_f15.png',
  height = 480 ,
  width = 480 ,
  resolution = 300,
  output=TRUE, 
  #col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c("#440154ff", '#21908dff', '#fde725ff'),
  cex = 0.6,
  fontfamily = "sans",
  cat.cex = 0.35,
  cat.default.pos = "outer",
  cat.pos = c(+0, 0, 0),
  cat.dist = c(0.055, 0.055, -0.455),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
  lty = 'blank'
)

venn.diagram(
  x = list(set1_down$gene_sym, set2_down$gene_sym, set3_down$gene_sym),
  category.names = c("Klf2KO", "Klf2KORing1bKD", "Ring1bKD"),
  filename = 'Venn_down.png',
  height = 480 ,
  width = 480 ,
  resolution = 300,
  output=TRUE, 
  #col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c("#440154ff", '#21908dff', '#fde725ff'),
  cex = 0.6,
  fontfamily = "sans",
  cat.cex = 0.35,
  cat.default.pos = "outer",
  cat.pos = c(+0, 0, 0),
  cat.dist = c(0.055, 0.055, -0.455),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
  lty = 'blank'
)

#교집합 출력
intersect_up <- intersect(set1_down$X, set2$X)
intersect_up <- intersect(intersect_up, set3_up$gene_sym)
intersect_down <- intersect(set1_down$X, set2_down$X)
intersect_down <- intersect(intersect_down, set3_down$gene_sym)
resvsDOXWT_intersect <- set1_down[(set1_down$X %in% intersect_up), ]
write.csv(resvsDOXWT_intersect,file="vsDOX_up_vs_WTdox_f15.csv")
write.csv(intersect_down,file="2D_vs_3D_intersect_down.csv")

intersect_down <- intersect(set1_down$X, set3_down$X)
intersect_down <- intersect(intersect_down, set3_down$X)
intersect_down <- intersect(intersect_down, set4_down$X)
intersect_down <- set1_down[set1_down$X %in% intersect_down, ]
write.csv(intersect_down,file="intersect_down_FK1D_f15.csv")

intersect_up <- intersect(set1_down$X, set2_down$X)
intersect_up <- intersect(intersect_up, set3_down$X)
intersect_up <- intersect(intersect_up, set4_down$X)
intersect_up <- set1_down[set1_down$X %in% intersect_up, ]
write.csv(intersect_up,file="intersect_up_FK12D_f15.csv")


#차집합 출력 
setdiff_Mutnon_up <- setdiff(set1_up$gene_sym, set2_up$gene_sym)
setdiff_Mutnon_down <- setdiff(set1_down$gene_sym, set2_down$gene_sym)
setdiff_MutDox_up <- setdiff(set2_up$gene_sym, set1_up$gene_sym)
setdiff_MutDox_down <- setdiff(set2_down$gene_sym, set1_down$gene_sym)
write.csv(setdiff_Mutnon_up,file="MutDox_vs_Mutnon_ctrl_WTnon_setdiff_Mutnon_up.csv")
write.csv(setdiff_Mutnon_down,file="MutDox_vs_Mutnon_ctrl_WTnon_setdiff_Mutnon_down.csv")
write.csv(setdiff_MutDox_up,file="MutDox_vs_Mutnon_ctrl_WTnon_setdiff_MutDox_up.csv")
write.csv(setdiff_MutDox_down,file="MutDox_vs_Mutnon_ctrl_WTnon_setdiff_MutDox_down.csv")

setdiff_up <- setdiff(set1_up$gene_sym, set2_up$gene_sym)
setdiff_down <- setdiff(set1_down$gene_sym, set2_down$gene_sym)
write.csv(setdiff_up,file="2D_vs_3D_setdiff_up.csv")
write.csv(setdiff_down,file="2D_vs_3D_setdiff_down.csv")
