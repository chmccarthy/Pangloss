## Install UpSetR if not already available.
#if (!require(UpSetR))
#{
#  install.packages("UpSetR")
#}

## Import and settings statements.
library(Cairo)
library(UpSetR)
x11(type="cairo")
args = commandArgs(trailingOnly=TRUE)
setEPS()
cairo_ps("UpSet.eps", height=8.5, width=8.5)

## This script is run from within Pangloss as "UpSet.R [matchtable] [tags]",
## so we parse these here.
if (length(args)==2) {
  matchtable <- args[1]
  tags <- args[2]
}

## Make dataframe for matchtable, assign the tags as column names.
## Note: have to read in matchtable without coercion to factors for the PAM transformation to work.
clusters <- read.table(matchtable, header = F, stringsAsFactors = F)
tags <- readLines(tags)
strains <- tags
colnames(clusters) <- strains

## Turn dataframe into PAM.
for (strain in strains)
{
  clusters[strain][clusters[strain] != "----------"] <- 1
  clusters[strain][clusters[strain] == "----------"] <- 0
}

## Convert the data within the dataframe to numerics from characters.
clusters <- as.data.frame(lapply(clusters, function(x) as.numeric(as.character(x))))

## Remove all core syntenic clusters from the dataframe, we're not interested in them for this plot.
## Note: there's a chance that the conversion step above may change symbols in column names,
##       so to avoid bugs from here on we use colnames(clusters) instead of strains as the
##       latter may not reflect the former.
clusters <- clusters[rowSums(clusters) < length(colnames(clusters)),]

## Generate UpSet plot for distribution of accessory orthologs in a pangenome.
upset(clusters, order.by = "freq", sets = rev(colnames(clusters)), keep.order = TRUE,
      mainbar.y.label = "Syntenic clusters", sets.x.label = "Accessory gene models")

## Close plot.
dev.off()

