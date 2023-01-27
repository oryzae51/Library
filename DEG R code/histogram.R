#histogram with frequency table
foldchange <- c(MLL1_vol_sub$log2FoldChange)
hist(foldchange,
     breaks = seq(-7,0,by=0.5),
     axes = FALSE, 
     main = "Histogram of MLL1-KO down-regulation gene", 
     xlab = "log2FoldChange")
x_axis_tick=seq(-7,0,by=0.5)
axis(side=1,at=x_axis_tick)
y_axis_tick=seq(0,max(h$counts),by=1000)
axis(side=2,at=y_axis_tick)
text(h$mids,h$counts+100,labels=h$counts)

library(ggplot2)
binsize <- diff(range(MLL1_vol_sub$log2FoldChange))/30
ggplot(MLL1_vol_sub, aes(x=log2FoldChange)) +
  geom_histogram(
  binwidth=0.5, 
  fill="white", 
  colour="black") +
  scale_x_continuous(breaks = seq(-8,7.5,by=0.5))


