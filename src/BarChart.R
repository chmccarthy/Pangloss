if (!require(ggplot2))
{
  install.packages("ggplot2")
}

library(ggplot2)
setEPS()

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

postscript("BarChart.eps", height=8.5, width=8.5)

df = read.table("/Users/charley/Dropbox/Genes_Paper/cluster_sizes.txt", header=TRUE, sep="\t")

scale <- df$Size

chao_label = sprintf("Observed N = %s\nChao N = %s", sum(df$Count), chao(df))

p <-ggplot(data = df, aes(x = Size, y = Count, fill = scale)) +
  geom_bar(stat = "identity") + scale_fill_gradient(low = "#FF8888",high = "#98FB98") +
  geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=-0.5) +
  annotate("text", x = length(df$Count) / 4, y = Inf, label = chao_label, 
           hjust = 0, vjust = 2.5) +
  scale_x_continuous(breaks = df$Size)

p
dev.off()
