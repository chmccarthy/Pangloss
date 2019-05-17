# install packages
install.packages(pkgs = c("ggplot2","dplyr"),dependencies = T)

library(ggplot2)
library(dplyr)

tags <- readLines("panoct_tags.txt")

df <- data.frame(BUSCO=character(), State=character(), Genome=character())

for (tag in tags)
{
  result <- read.table(stringr::str_interp("busco/run_${name}.busco/full_table_${name}.busco.tsv", list(name = tag)), skip = 4, fill=TRUE)
  result <- select(result, V2, V3)
  colnames(result) <- c("State", "Genome")
  result$Genome <- tag
  df <- rbind(df, result)
}



p <- ggplot(df, aes(x = BUSCO, y = Genome, fill = State)) + geom_tile(colour = "white", size = 0.05) +
     theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
     scale_fill_manual(values=rev(c("#f1eef6", "#bdc9e1", "#74a9cf", "#0570b0")))

p