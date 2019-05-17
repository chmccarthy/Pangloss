## Install UpSetR if not already available.
if (!require(ggplot2))
{
  install.packages("ggplot2")
}

if (!require(reshape2))
{
  install.packages("reshape2")
}

if (!require(dyplr))
{
  install.packages("dplyr")
}

## Import and settings statements.
library(ggplot2)
library(reshape2)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
setEPS()
postscript("Clusters.eps", height=8.5, width=8.5)

## This script is run from within Pangloss as "UpSet.R [matchtable] [tags]",
## so we parse these here.
if (length(args)==2) {
  matchtable <- args[1]
  tags <- args[2]
  karyotypes <- read.table(args[3])
}

matchtable <- "matchtable.txt"
tags <- "panoct_tags.txt"
karyotypes <- "karyotypes.txt"

## Make dataframe for matchtable, assign the tags as column names.
## Note: have to read in matchtable without coercion to factors for the PAM transformation to work.
clusters <- read.table(matchtable, header = F, stringsAsFactors = F)
k <- read.table(karyotypes)
tags <- readLines(tags)
strains <- tags
colnames(clusters) <- strains

## Turn dataframe into PAM.
for (strain in strains)
{
  clusters[strain][clusters[strain] != "----------"] <- 1
  clusters[strain][clusters[strain] == "----------"] <- 0
}

melted <- melt(t(clusters), value.name = "Presence", 
               varnames = c("Strain", "Cluster"))
p <- ggplot(melted, aes(x = Cluster, y = Strain,
                        fill = Presence)) + geom_tile(colour="white",size=0.25)

cluster_ids <- unique(melted$Cluster)
