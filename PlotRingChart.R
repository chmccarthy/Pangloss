library(ggplot2)
library(ggrepel)
args = commandArgs(trailingOnly=TRUE)
setEPS()

if (length(args)==5) {
  core <- as.numeric(args[1])
  softcore <- as.numeric(args[2])
  accessory <- as.numeric(args[3])
  sizes <- strsplit(args[4], ",")
  counts <- as.numeric(unlist(strsplit(args[5], ",")))
}

postscript("pangenome.eps", height = 8.5, width = 8.5)
dat = data.frame(count=c(accessory, softcore, core), category=c("Accessory", "Softcore", "Core"))

dat$fraction = dat$count / sum(dat$count)
dat = dat[order(dat$fraction), ]
dat$ymax = cumsum(dat$fraction)
dat$ymin = c(0, head(dat$ymax, n=-1))

p1 = ggplot(dat, aes(fill=category, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
     geom_rect(stat="identity", color='black') +
     coord_polar(theta="y") +
     xlim(c(0, 6)) +
     theme(panel.grid=element_blank()) +
     theme(axis.text=element_blank()) +
     theme(axis.ticks=element_blank()) +
     annotate("text", x = 0, y = 0, label = "Teh") +
     labs(title="") +
     geom_label_repel(aes(label=paste(round(fraction*100, digits = 2),"%"), x=4 ,y=(ymin+ymax)/2), nudge_x = 2, arrow = arrow(angle = 0, length = unit(0.01, 'npc')), inherit.aes = FALSE, show.legend = FALSE, direction="both", max.iter = 3e3, ylim=c(6,8)) +
     guides(fill=guide_legend(override.aes=list(colour=NA)))

p1
dev.off()

postscript("cluster_sizes.eps", height = 8.5, width = 8.5)
dat2 = data.frame(count=c(counts), category=unlist(sizes))

dat2$fraction = dat2$count / sum(dat2$count)
dat2 = dat2[order(dat2$fraction), ]
dat2$ymax = cumsum(dat2$fraction)
dat2$ymin = c(0, head(dat2$ymax, n=-1))

p2 = ggplot(dat2, aes(fill=category, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
     geom_rect(stat="identity", color='black') +
     coord_polar(theta="y") +
     xlim(c(0, 6)) +
     theme(panel.grid=element_blank()) +
     theme(axis.text=element_blank()) +
     theme(axis.ticks=element_blank()) +
     annotate("text", x = 0, y = 0, label = "Teh") +
     labs(title="") +
     geom_label_repel(aes(label=paste(round(fraction*100, digits = 2),"%"), x=4 ,y=(ymin+ymax)/2), nudge_x = 2, arrow = arrow(angle = 0, length = unit(0.01, 'npc')), inherit.aes = FALSE, show.legend = FALSE, direction="both", max.iter = 3e3, ylim=c(6,8)) +
     guides(fill=guide_legend(override.aes=list(colour=NA)))

p2
dev.off()

q()
