## Install ggplot2 if not already available.
#if (!require(ggplot2))
#{
#  install.packages("ggplot2")
#}

## Function for Chao lower bound estimate of pangenome size based on a similar function in micropan.
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

## Function for potenital Bohning lower bound estimate of pangenome size if I can figure out how to implement it?
bohning <- function(df)
{
  y <- sum(df$count)
  y2 <- df$count[2]
  y3 <- df$count[3]
  if (y3 == 0)
  {
    stop()
  }
  else
  {
    bohning <- round(y + y1^3 / (y2^2))
  }
  return(bohning)
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
  scale_x_continuous(breaks = df$Size)

## Write bar chart plot to file and close.
p
dev.off()
