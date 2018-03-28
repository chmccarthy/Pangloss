library(UpSetR)
setEPS()

clusters<-read.table("noncore_matchtable.txt", header = F)
postscript("./noncore_upset.eps", height = 8.5, width = 8.5)
colnames(clusters) <- readLines("panoct_tags.txt")
upset(clusters, order.by = "freq", sets=rev(readLines("panoct_tags.txt")), keep.order = TRUE, mainbar.y.label = "Number of clusters", sets.x.label = "Number of proteins")

dev.off()
q()
