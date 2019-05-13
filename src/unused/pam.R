setEPS()
library(gplots)
library(RColorBrewer)

t = read.table("~/Desktop/t.csv", sep="\t")
c = readLines("~/Desktop/c.txt")
r = readLines("~/Desktop/r.txt")

postscript("~/Dropbox/pam.eps", width=16.5)
colnames(t) <- r
rownames(t) <- c
t <- t(t)
my_pal <- colorRampPalette(c("#FF0000", "#00B050"))(n=2)
heatmap.2(as.matrix(t), col=my_pal, Rowv=FALSE, Colv = FALSE, dendrogram="none", trace="none",  key = FALSE, lhei = c(1,8), lwid = c(1,8))
dev.off()