## Install BiocManager and karyploteR if not already available.
#if (!require(karyoploteR))
#{
#  if (!require(BiocManager))
#  {
#    install.packages("BiocManager")
#  }
#  BiocManager::install("karyoploteR", version = "3.9")
#}

# Function to plot color bar, modified from https://gist.github.com/johncolby/993830.
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut))/(max)
  plot(c(0,10), c(min,max+1), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
  title("\n\n# of syntenic orthologs", cex.main = 1, font.main= 1)
  axis(2, ticks, las=1, at=(c(min+0.5, max+0.5)), labels=(c(min, max)), mgp=c(1.5,0.5,0))
  for (i in 1:(length(lut))) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }	
}

## Import and settings statements (we set individual karyotype plots within the loop).
library(karyoploteR)
library(Hmisc)
setEPS()
args = commandArgs(trailingOnly = TRUE)

## This script is run from within Pangloss as "Karyotype.R [tags] [karyotypes] [lengths]",
## so the cluster size file is parsed in as sizes.
if (length(args)==3)
{
  tags <- readLines(args[1])
  karyotypes <- read.table(args[2])
  lengths <- read.table(args[3])
}
  
## Loop through every strain in a dataset and run karyoplotting for that strain.
for (tag in tags)
{
  ## Open postscript file and extract strain-specific data from karyotype dataframe.
  genes = karyotypes[karyotypes$V6 == tag, ]
  contigs = lengths[lengths$V4 == tag, ]

  ## Extract relevant data from strain dataframes into new dataframes.
  names(genes) <- c("chr", "name", "start", "end", "component", "tag", "orthologs")
  names(contigs) <- c("chr", "start", "end", "tag")
  contigs = within(contigs, rm("tag"))
  comps = within(genes, rm("tag", "name", "orthologs"))
  ortho = within(genes, rm("tag", "name", "component"))

  ## Convert karyotype and length data to GRanges.
  custom.genome <- toGRanges(contigs)
  custom.compbands <- toGRanges(comps)
  custom.orthobands <- toGRanges(ortho)
  
  ## Generate continuous colour range for ortholog number plots.
  sets <- sort(unique(ortho$orthologs))
  colfunc <- colorRampPalette(c("red", "yellow", "darkgreen"))
  colscale <- colfunc(length(sets))
  
  ## Generate and write component karyotype plot and then close file.
  postscript(stringr::str_interp("${name}_components.eps", list(name = tag)), height=8.5, width=8.5)
  pp <- getDefaultPlotParams(plot.type=1)
  pp$leftmargin = 0.15
  kp <- plotKaryotype(genome = custom.genome, plot.type = 1, plot.params = pp)
  kpAddMainTitle(kp, main=tag, col="black")
  kpPlotRegions(kp, data = custom.compbands[custom.compbands$component == "core"], avoid.overlapping = TRUE, data.panel = "ideogram", col = "darkgreen", layer.margin = 0.15)
  kpPlotRegions(kp, data = custom.compbands[custom.compbands$component == "acc"], avoid.overlapping = TRUE, data.panel = "ideogram", col = "red", layer.margin = 0.15)
  dev.off()
  
  ## Generate and write ortholog number karyotype plot and then close file.
  postscript(stringr::str_interp("${name}_orthologs.eps", list(name = tag)), height=8.5, width=8.5)
  pp <- getDefaultPlotParams(plot.type=1)
  pp$leftmargin = 0.15
  kp <- plotKaryotype(genome = custom.genome, plot.type = 1, plot.params = pp)
  kpAddMainTitle(kp, main=tag, col="black")
  for (set in rev(sets))
  {
    kpPlotRegions(kp, data = custom.orthobands[custom.orthobands$orthologs == set], avoid.overlapping = TRUE, data.panel = "ideogram", col = colscale[set], layer.margin = 0.15)
  }
  subplot(color.bar(colscale, 1, length(tags), nticks=length(tags), title="Syntenic orthologs"), 0.9, 0.9)
  dev.off()
}

