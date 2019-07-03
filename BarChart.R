## Install ggplot2 if not already available.
#if (!require(ggplot2))
#{
#  install.packages("ggplot2")
#}

## Function for Chao lower bound estimate of pangenome size.
## Adapted from micropan's chao function by Snipen & Liland (see https://rdrr.io/cran/micropan/src/R/powerlaw.R).
chao <- function(df)
{
  y <- sum(df$Count)
  y1 <- df$Count[1]
  y2 <- df$Count[2]
  if (y2 == 0)
  {
    stop()
  }
  else
  {
    chao <- round(y + y1^2 / (2*y2))
  }
  return(chao)
}

## Import and settings statements.
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
setEPS()
postscript("BarChart.eps", height=8.5, width=8.5)

## This script is run from within Pangloss as "BarChart.R [name_of_cluster_size_file]",
## so the cluster size file is parsed in as sizes.
if (length(args)==1)
{
  sizes <- args[1]
}

## Read in cluster size file as dataframe.
df = read.table(sizes, header=TRUE, sep="\t")

## Get scale from size of dataframe (i.e. number of cluster sizes, number of strains in your dataset).
scale <- df$Size

## Run chao function and generate label for observed pangenome size and estimated pangenome size.
chao_label = sprintf("Observed N = %s\nChao N = %s", sum(df$Count), chao(df))

## Plot bar chart with pangenome sizes label and 1-to-n gradient scaling.
p <-ggplot(data = df, aes(x = Size, y = Count, fill = scale)) +
  geom_bar(stat = "identity") + scale_fill_gradient(low = "#FF8888",high = "#98FB98") +
  annotate("text", x = length(df$Count) / 4, y = Inf, label = chao_label, 
           hjust = 0, vjust = 2.5) +
  geom_text(aes(label = as.character(Count)), vjust=-0.5) +
  scale_x_continuous(breaks = df$Size) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, sum(df$Count))) +
  labs(x = "Cluster size", y = "# of clusters", fill = "Size") +
  theme_classic()

## Write bar chart plot to file and close.
p
dev.off()
