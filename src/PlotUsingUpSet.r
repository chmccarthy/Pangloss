library(UpSetR)
args = commandArgs(trailingOnly=TRUE)
setEPS()

clusters<-read.table(args[1], header = F)
postscript(args[2], height = 8.5, width = 8.5)
colnames(clusters) <- readLines("panoct_tags.txt")
upset(clusters, order.by = "freq", sets=rev(c("A1163", "AF10", "ISSFT021", "HMRAF706", "IF1SWF4", "AF293", "LMB35AA", "RP2014", "JCM10253", "Z5", "AF210", "HMRAF270")), keep.order = TRUE, mainbar.y.label = "Number of clusters", sets.x.label = "Number of proteins")

dev.off()
q()

#readLines("panoct_tags.txt")