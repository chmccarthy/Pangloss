library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
setEPS()

if (length(args)==2) {
  core <- max(as.numeric(args[1]), as.numeric(args[2]))
  accessory <- min(as.numeric(args[1]), as.numeric(args[2]))
}


postscript("pangenome.eps", height = 8.5, width = 8.5)
dat = data.frame(count=c(accessory, core), category=c("Accessory", "Core"))

dat$fraction = dat$count / sum(dat$count)
dat = dat[order(dat$fraction), ]
dat$ymax = cumsum(dat$fraction)
dat$ymin = c(0, head(dat$ymax, n=-1))


p1 = ggplot(dat, aes(fill=category, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
     geom_rect() +
     coord_polar(theta="y") +
     xlim(c(0, 4)) +
     theme(panel.grid=element_blank()) +
     theme(axis.text=element_blank()) +
     theme(axis.ticks=element_blank()) +
     annotate("text", x = 0, y = 0, label = "Teh") +
     labs(title="") +
     geom_label(aes(label=paste(round(fraction*100, digits = 2),"%"),x=3.5,y=(ymin+ymax)/2),inherit.aes = FALSE, show.legend = FALSE)

p1     
dev.off()
q()