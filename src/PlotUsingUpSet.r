library(UpSetR)

clusters<-read.table("~/Documents/GitHub/PanGLOSS/work/Yeast/acc_pam.txt", header = F)
sets <- readLines("~/Documents/GitHub/PanGLOSS/work/Yeast/panoct_tags.txt")
colnames(clusters) <- sets
orders <- sets
upset(clusters, order.by = "freq", sets = rev(orders), keep.order = TRUE,
      mainbar.y.label = "Number of clusters", sets.x.label = "Number of proteins",
      mb.ratio = c(1, 0), nintersects=100)

#readLines("panoct_tags.txt")


