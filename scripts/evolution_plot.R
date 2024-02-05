
library(ggplot2)
library(cowplot)


###################### SUPP FIG 6A #####################


muestras = c("85","107","161","185","207","243","261")
lista = list()
myplots <- vector('list', length(muestras))
i = 1
for ( mu in muestras) {
  data2 <- read.csv2(paste0("/home/jmartin/Documents/articulo/scripts/inputs/Supp_FIG6/counts_", mu, ".txt"), header =T, dec='.', sep=',')
  data2$names = factor(data2$names, levels = rev(c("tissue", "plasma_bl","tissue:plasma_bl",  "tissue:plasma_po" , "tissue:plasma_bl:plasma_po",
                                                   "plasma_bl:plasma_po", "tissue:plasma_re",  "plasma_bl:plasma_re", 
                                                   "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re",
                                                   "tissue:plasma_po:plasma_re",
                                                   "plasma_bl:plasma_po:plasma_re", 
                                                   "plasma_po", 
                                                   "plasma_po:plasma_re", 
                                                   "plasma_re")))
  
  data2$tipo = factor(data2$tipo, levels = c("Tissue", "Plasma-BL", "Plasma-PO","Plasma-RE"))
  
  col <- c("red4", "lightpink2",  "firebrick1",   "plum3", "slategray3",   "slategray1", "steelblue4", "turquoise3",  
           "deepskyblue3",  "aquamarine4",  "#91D1C2FF",  "darkseagreen2", "#E64B3533", "#B09C8599","#B09C85FF")
  g = ggplot(data2, aes(fill=names, y=len, x=tipo)) + 
    geom_bar(position="stack", stat="identity") + scale_fill_manual(values = col) + xlab(mu) + ylim(0, 1500)
  
  myplots[[i]] <- g
  i = i + 1
  
}

grobs <- ggplotGrob(g)$grobs 
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
a = plot_grid(myplots[[1]]  + theme(legend.position="none"), myplots[[2]] + theme(legend.position="none"),
              myplots[[3]] + theme(legend.position="none"), myplots[[4]] + theme(legend.position="none"), 
              myplots[[5]]+ theme(legend.position="none"),
              myplots[[6]]+ theme(legend.position="none"), myplots[[7]]+ theme(legend.position="none")+ theme(legend.position="none"),  align = "v", 
              axis = "lr",
              nrow = 4)


svg(file="/home/jmartin/Documents/articulo/scripts/figuras/Supp_Fig6a.svg", width = 20, height = 12)
plot_grid(a , legend, align = "h", 
          axis = "lr",
          ncol = 2,  rel_widths = c(6/8, 2/8))
dev.off()



###################### SUPP FIG 6B #####################


muestras = c("M1","M2","M3","M5","M6","M7","M8","M9","M10","M11","M12","M13","M14","M15")
lista = list()
myplots <- vector('list', length(muestras))
i = 1
for ( mu in muestras) {
  
  data2 <- read.csv2(paste0("/home/jmartin/Documents/articulo/scripts/inputs/Supp_FIG6/counts_", mu, ".txt"), header =T, dec='.', sep=',')
  data2$names = factor(data2$names, levels = rev(c("tissue", "plasma_bl","tissue:plasma_bl",  "tissue:plasma_po" , "tissue:plasma_bl:plasma_po",
                                                   "plasma_bl:plasma_po", "tissue:plasma_re",  "plasma_bl:plasma_re", 
                                                   "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re",
                                                   "tissue:plasma_po:plasma_re",
                                                   "plasma_bl:plasma_po:plasma_re", 
                                                   "plasma_po", 
                                                   "plasma_po:plasma_re", 
                                                   "plasma_re")))
  
  data2$tipo = factor(data2$tipo, levels = c("Tissue", "Plasma-BL", "Plasma-PO","Plasma-RE"))
  
  col <- c("red4", "lightpink2",  "firebrick1",   "plum3", "slategray3",   "slategray1", "steelblue4", "turquoise3",  
           "deepskyblue3",  "aquamarine4",  "#91D1C2FF",  "darkseagreen2", "#E64B3533", "#B09C8599","#B09C85FF")
  g = ggplot(data2, aes(fill=names, y=len, x=tipo)) + 
    geom_bar(position="stack", stat="identity") + scale_fill_manual(values = col) + xlab(mu) + ylim(0,2500)
  
  myplots[[i]] <- g
  i = i + 1
  
}

grobs <- ggplotGrob(g)$grobs 
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
a = plot_grid(myplots[[1]]  + theme(legend.position="none"), myplots[[2]] + theme(legend.position="none"),
              myplots[[3]] + theme(legend.position="none"), myplots[[4]] + theme(legend.position="none"), 
              myplots[[5]]+ theme(legend.position="none"),
              myplots[[6]]+ theme(legend.position="none"), myplots[[7]]+ theme(legend.position="none"),
              myplots[[8]]+ theme(legend.position="none"), myplots[[9]]+ theme(legend.position="none"),
              myplots[[10]]+ theme(legend.position="none"), myplots[[11]]+ theme(legend.position="none"),
              myplots[[12]]+ theme(legend.position="none"), myplots[[13]]+ theme(legend.position="none"),
              myplots[[14]]+ theme(legend.position="none"),  align = "v", 
          axis = "lr",
          nrow = 4)

svg(file="/home/jmartin/Documents/articulo/scripts/figuras/Supp_Fig6B.svg", width = 20, height = 12)
plot_grid(a , legend, align = "h", 
          axis = "lr",
          ncol = 2,  rel_widths = c(6/8, 2/8))

dev.off()
