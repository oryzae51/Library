E2F4 <- read.csv("/Users/hmkim/data/deg_data/target_gene_list/E2F4.csv", header = FALSE, sep = ",")
E2F6 <- read.csv("/Users/hmkim/data/deg_data/target_gene_list/E2F6.csv", header = FALSE, sep = ",")
MLL1 <- read.csv("/Users/hmkim/data/deg_data/MLL-KO_down/mll1down_genesym.csv", header = FALSE, sep = ",")
MLL2 <- read.csv("/Users/hmkim/data/deg_data/MLL-KO_down/mll2down_genesym.csv", header = FALSE, sep = ",")
MLL4 <- read.csv("/Users/hmkim/data/deg_data/MLL-KO_down/mll4down_genesym.csv", header = FALSE, sep = ",")

#intersect
E2F4_MLL1_intersect <- length(intersect((E2F4$V1), (MLL1$V1)))
E2F4_MLL2_intersect <- length(intersect((E2F4$V1), (MLL2$V1)))
E2F4_MLL4_intersect <- length(intersect((E2F4$V1), (MLL4$V1)))
E2F6_MLL1_intersect <- length(intersect((E2F6$V1), (MLL1$V1)))
E2F6_MLL2_intersect <- length(intersect((E2F6$V1), (MLL2$V1)))
E2F6_MLL4_intersect <- length(intersect((E2F6$V1), (MLL4$V1)))

MLLintersect <- intersect(c(resLFC_MLL1df$Geneid), c(resLFC_MLL2df$Geneid))