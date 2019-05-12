## Install BiocManager and karyploteR if not already available.
if (!require(karyoploteR))
{
  if (!require(BiocManager))
  {
    install.packages("BiocManager")
  }
  BiocManager::install("karyoploteR", version = "3.8")
}

## Import and settings statements (we set individual karyotype plots within the loop).
library(karyoploteR)
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
  postscript(stringr::str_interp("${name}.eps", list(name = tag)), height=8.5, width=8.5)
  genes = karyotypes[karyotypes$V6 == tag, ]
  contigs = lengths[lengths$V4 == tag, ]

  ## Extract relevant data from strain dataframes into new dataframes.
  names(genes) <- c("chr", "name", "start", "end", "component", "tag")#, "orthologs")
  names(contigs) <- c("chr", "start", "end", "tag")
  contigs = within(contigs, rm("tag"))
  genes = within(genes, rm("tag", "name"))
  pg <- c(core="00ff00", acc="ff0000")

  ## Convert karyotype and length data to GRanges.
  custom.genome <- toGRanges(contigs)
  custom.cytobands <- toGRanges(genes)

  ## Generate and write karyotype plot and then close file.
  pp <- getDefaultPlotParams(plot.type=1)
  pp$leftmargin = 0.15
  kp <- plotKaryotype(genome = custom.genome, plot.type = 1, plot.params = pp)
  kpAddMainTitle(kp, main=tag, col="black")
  kpPlotRegions(kp, data = custom.cytobands[custom.cytobands$component == "core"], avoid.overlapping = TRUE, data.panel = "ideogram", col = "darkgreen")
  kpPlotRegions(kp, data = custom.cytobands[custom.cytobands$component == "acc"], avoid.overlapping = TRUE, data.panel = "ideogram", col = "red")
  dev.off()
  }

