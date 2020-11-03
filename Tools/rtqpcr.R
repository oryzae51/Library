#Load qpcr table data
qpcr <- read.table("/Users/hmkim/data/rapid.tsv", header = TRUE, sep = "\t")

mfc <- c(qpcr$Mean.fold.change)
mfc <- mfc[!is.na(mfc)]

sample <- c("WT", "MLL1-KO", "MLL2-KO", "MLL4-KO")

target <- unique(c(qpcr$Target))
target <- target[!target %in% "GAPDH"]

sem <- c(qpcr$SEM)
sem <- sem[!is.na(sem)]

CHK1 <- data.frame(Sample=sample, MeanFoldChange=c(mfc[1:4]), SEM=c(sem[1:4]))
MCM5 <- data.frame(Sample=sample, MeanFoldChange=c(mfc[5:8]), SEM=c(sem[5:8]))
MCM6 <- data.frame(Sample=sample, MeanFoldChange=c(mfc[9:12]), SEM=c(sem[9:12]))
MCM7 <- data.frame(Sample=sample, MeanFoldChange=c(mfc[13:16]), SEM=c(sem[13:16]))
ORC1 <- data.frame(Sample=sample, MeanFoldChange=c(mfc[17:20]), SEM=c(sem[17:20]))
RFC1 <- data.frame(Sample=sample, MeanFoldChange=c(mfc[21:24]), SEM=c(sem[21:24]))
RRM2 <- data.frame(Sample=sample, MeanFoldChange=c(mfc[25:28]), SEM=c(sem[25:28]))

write.csv(as.data.frame(CHK1), file = "~/data/CHK1.csv", row.names = FALSE)
write.csv(as.data.frame(MCM5), file = "~/data/MCM5.csv", row.names = FALSE)
write.csv(as.data.frame(MCM6), file = "~/data/MCM6.csv", row.names = FALSE)
write.csv(as.data.frame(MCM7), file = "~/data/MCM7.csv", row.names = FALSE)
write.csv(as.data.frame(ORC1), file = "~/data/ORC1.csv", row.names = FALSE)
write.csv(as.data.frame(RFC1), file = "~/data/RFC1.csv", row.names = FALSE)
write.csv(as.data.frame(RRM2), file = "~/data/RRM2.csv", row.names = FALSE)
