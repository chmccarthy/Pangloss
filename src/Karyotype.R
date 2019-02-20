library(karyoploteR)

att = read.table("~/Dropbox/Genes_Paper/allatt.db")
tags = readLines("~/Dropbox/Genes_Paper/panoct_tags.txt")
lengths = read.table("~/Dropbox/Genes_Paper/genomes/lengths.txt")
matchtable = read.table("~/Dropbox/Genes_Paper/panoct/matchtable.txt")
names(matchtable) <- tags

for (tag in tags)
{
  genes = att[att$V6 == tag, ]
  contigs = lengths[lengths$V4 == tag, ]
  
  names(genes) <- c("chr", "name", "start", "end", "info", "tag")
  genes = within(genes, rm("info", "tag"))
  genes = genes[c("chr", "start", "end", "name")]
  for (gene in genes)
  {
    print(gene[3])
  }
  
  names(contigs) <- c("chr", "start", "end", "tag")
  contigs = within(contigs, rm("tag"))
  
  custom.genome <- toGRanges(contigs)
  custom.cytobands <- toGRanges(genes)
  
  kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, plot.type=1)
  kpAddBaseNumbers(kp, tick.dist = 1000000)
  break()
}