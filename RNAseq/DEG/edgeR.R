library(edgeR)

group <- factor(c(2, 2, 3, 3, 4, 4, 1, 1))
y <- DGEList(counts = cnts, group = group)
keep <- filterByExpr(y)
filtered_y <- y[keep, , keep.lib.sizes = FALSE]
norm_y <- calcNormFactors(filtered_y)
design <- model.matrix(~group)
dispers_y <- estimateDisp(norm_y, design)

#Perform quasi-likelihood F-test
fit <- glmQLFit(dispers_y, design)
qlf.2vs1 <- glmQLFTest(fit, coef = 2)
topTags(qlf.2vs1)
FDR <- p.adjust(qlf.2vs1$table$PValue)
sum(FDR<0.01)
edge_sig <- subset(edge_sig, logFC>1 | logFC< -1)
edge_sig <- subset(edge_sig, logFC>1)
FDR <- p.adjust(edge_sig$PValue)
edge_sig <- cbind(edge_sig, FDR)
sum(FDR<0.01)
edge_sig <- subset(edge_sig, FDR<0.01)
edge_sig <- subset(qlf.2vs1$table, FDR < 0.01)
write.csv(as.data.frame(edge_sig), 
          file="~/data/deg_confirm/edge_mine_up.csv")

library(venn)
library(gplots)
prot1 <- read.csv("~/data/deg_confirm/edge_mine_up.csv", header = FALSE)
prot2 <- read.csv("~/data/deg_confirm/deseq2_mine_up.csv", header = FALSE)
data <- list(A = prot1$V1, B = prot2$V1)
venn(data)

library(venn)
library(gplots)
prot1 <- read.csv("~/data/deg_confirm/edge_mine_up_down_rmdot.csv", header = FALSE)
prot2 <- read.csv("~/data/deg_confirm/deseq2_mine_up_down_rmdot.csv", header = FALSE)
data <- list(KHM_pipeline_edge_up_and_down_reg = prot1$V1, KHM_pipeline_deseq2_up_and_down_reg = prot2$V1)
venn(data)
