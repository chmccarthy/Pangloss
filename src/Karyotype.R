if (!require(karyoploteR))
{
  if (!require(BiocManager))
  {
    install.packages("BiocManager")
  }
  BiocManager::install("karyoploteR", version = "3.8")
}

library(karyoploteR)
setEPS()

args = commandArgs(trailingOnly = TRUE)

tags <- readLines(args[1])
karyotypes <- read.table(args[2])
lengths <- read.table(args[3])


for (tag in tags)
{
  postscript(stringr::str_interp("${name}.eps", list(name = tag)), height=8.5, width=8.5)
  genes = karyotypes[karyotypes$V5 == tag, ]
  contigs = lengths[lengths$V4 == tag, ]
  
  names(genes) <- c("chr", "start", "end", "component", "tag")#, "orthologs")
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
  dev.off()
  }

