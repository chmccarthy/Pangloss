library(karyoploteR)
setEPS()

tags = readLines("~/Dropbox/Genes_Paper/panoct_tags.txt")
karyotypes = read.table("~/Dropbox/Genes_Paper/karyotypes.txt")
lengths = read.table("~/Dropbox/Genes_Paper/genomes/lengths.txt")

for (tag in tags)
{
  postscript("test.eps", height=8.5, width=8.5)
  genes = karyotypes[karyotypes$V5 == tag, ]
  contigs = lengths[lengths$V4 == tag, ]
  
  names(genes) <- c("chr", "start", "end", "component", "tag", "orthologs")
  names(contigs) <- c("chr", "start", "end", "tag")
  contigs = within(contigs, rm("tag"))
  genes = within(genes, rm("tag"))
  pg <- c(core="00ff00", acc="ff0000")
  
  custom.genome <- toGRanges(contigs)
  custom.cytobands <- toGRanges(genes)
  
  pp <- getDefaultPlotParams(plot.type=1)
  pp$leftmargin = 0.15
  kp <- plotKaryotype(genome = custom.genome, plot.type = 1, plot.params = pp)
  kpAddMainTitle(kp, main=tag, col="black")
  kpPlotRegions(kp, data = custom.cytobands[custom.cytobands$component == "core"], avoid.overlapping = TRUE, data.panel = "ideogram", col = "darkgreen")
  kpPlotRegions(kp, data = custom.cytobands[custom.cytobands$component == "acc"], avoid.overlapping = TRUE, data.panel = "ideogram", col = "red")
  kpBars(kp, y0=start(custom.cytobands[custom.cytobands$orthologs]), y1=end(custom.cytobands[custom.cytobands$orthologs]), data.panel = 1)
  dev.off()
  }

