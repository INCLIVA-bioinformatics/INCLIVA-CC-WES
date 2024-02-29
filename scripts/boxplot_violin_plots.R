library(ggplot2)
library(ggpubr)
# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


####### FIGURE 2B ###########

vafs <- read.csv2("/home/jmartin/Documents/articulo/scripts/inputs/input_FIG2B.csv", header =T, dec='.', sep='\t')
#vafs$

vafs$PORC = vafs$PORC*100

svg(file="/home/jmartin/Documents/articulo/scripts/figuras/FIG2B.svg", width = 15, height = 10)

p <- ggplot(vafs, aes(x=TIPO, y=PORC, color=TIPO, )) + 
  geom_violin() + stat_summary(fun.data=data_summary) + theme_bw() + ylim(0,100) + 
  stat_compare_means( aes(label = ..p.signif..), method = "wilcox.test", paired = TRUE, label.x = 1.5, label.y = 100) +
  geom_line(aes(group = ID), color = "grey") + scale_color_manual(values=c("#3C5488FF", "#E64B35FF")) + ylab("% of concordance") + 
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.position = "none",
        axis.text=element_text(size=25),
        axis.title=element_text(size=25,face="bold"),
        axis.title.x=element_blank())
p
dev.off()


################ FIGURE 2C #################


vafs <- read.csv2("/home/jmartin/Documents/articulo/scripts/inputs/input_FIG2C.csv", header =T, dec='.', sep='\t')
#vafs$

vafs$PORC = vafs$PORC*100

svg(file="/home/jmartin/Documents/articulo/scripts/figuras/FIG2C.svg", width = 15, height = 10)

p <- ggplot(vafs, aes(x=TIPO, y=PORC, color=TIPO, )) + 
  geom_violin() + stat_summary(fun.data=data_summary) + theme_bw() + ylim(0,100) + 
  stat_compare_means( aes(label = ..p.signif..), method = "wilcox.test", paired = TRUE, label.x = 1.5, label.y = 100) +
  geom_line(aes(group = ID), color = "grey") + scale_color_manual(values=c("#3C5488FF", "#E64B35FF")) +
  ylab("% of concordance") + 
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.position = "none",
        axis.text=element_text(size=25),
        axis.title=element_text(size=25,face="bold"),
        axis.title.x=element_blank())
p
dev.off()



########### Supp FIG A1 #########


vafs <- read.csv2("/home/jmartin/Documents/articulo/scripts/inputs/input_Supp_Fig1A.txt", header =T, dec='.', sep='\t')

vafs$PORC = vafs$PORC*100

svg(file="/home/jmartin/Documents/articulo/scripts/figuras/Supp_FIG1A.svg", width = 14, height = 10)
p <- ggplot(vafs, aes(x=TIPO, y=PORC, color=TIPO)) + 
  geom_violin() + stat_summary(fun.data=data_summary) + theme_bw() + ylim(0, 100) + scale_colour_manual(values =c("steelblue","salmon","darkolivegreen3")) +
  ylab("% of concordance") + 
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.position = "none",
        axis.text=element_text(size=25), axis.title.x=element_blank(),
        axis.title=element_text(size=25,face="bold"))
p
dev.off()

############### SUPPL FIGURE 5B ############## 


vafs <- read.csv2("/home/jmartin/Documents/articulo/scripts/inputs/input_Supp_FIG5B.csv", header =T, dec='.', sep='\t')
#vafs$

vafs$PORC = vafs$PORC*100

vafs$TIPO = factor(vafs$TIPO, levels = c("PRIMARY TISSUE-PLASMA RELAPSE", "PLASMA BASELINE-PLASMA RELAPSE"))

svg(file="/home/jmartin/Documents/articulo/scripts/figuras/Supp_FIG5B.svg", width = 16, height = 10)
p <- ggplot(vafs, aes(x=TIPO, y=PORC, color=TIPO, )) + 
  geom_violin() + stat_summary(fun.data=data_summary) + theme_bw() + ylim(0,100) + 
  stat_compare_means( aes(label = ..p.signif..), method = "wilcox.test", paired = TRUE, label.x = 1.5, label.y = 100) +
  geom_line(aes(group = ID), color = "grey") + scale_color_manual(values=c("#3C5488FF", "#E64B35FF")) + ylab("% of concordance") + 
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.position = "none",
        axis.text=element_text(size=25),
        axis.title=element_text(size=25,face="bold"),
        axis.title.x=element_blank())
p
dev.off()



############# Supp FIG3 ##########

vafs <- read.csv2("/home/jmartin/Documents/articulo/scripts/inputs/input_Supp_FIG3.csv", header =T, dec='.', sep='\t')


svg(file="/home/jmartin/Documents/articulo/scripts/figuras/Supp_FIG3.svg", width = 11, height = 10)
p <- ggplot(vafs, aes(x=TIPO, y=CONCORDANCE, fill = TIPO)) + 
  geom_boxplot() + theme_bw() +
  ylab("Median Concordance") + scale_fill_manual(values =c("#3C5488FF","#E64B35FF")) + ylim(0,1) +
  stat_compare_means(aes(group = TIPO), label = "p.signif", label.x = 1.5, label.y = 0.95) +
  theme(plot.title = element_text(size = 20, face = "bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.title.x=element_blank())
p

dev.off()




###### Supp FIG 10A ##########

vafs <- read.csv2("/home/jmartin/Documents/articulo/scripts/inputs/input_Supp_FIG10A.txt", header =T, dec='.', sep='\t')



vafs$TIPO <- factor( as.character(vafs$TIPO), levels=c("WBCs", "TISSUE", "PLASMA-BL", "PLASMA-PO", "PLASMA") )



svg(file="/home/jmartin/Documents/articulo/scripts/figuras/Supp_FIG10A.svg", width = 13, height = 10)
p <- ggplot(vafs, aes(x=TIPO, y=COV, color=INS, fill=INS)) + 
  geom_boxplot() + theme_bw() + scale_colour_manual(values =c("grey48","grey50")) +
  scale_fill_manual(values =c("#3C5488FF","#E64B35FF")) +
  ylab("Median Coverage") + 
  stat_compare_means(aes(group = INS), label = "p.signif") +
  theme(plot.title = element_text(size = 20, face = "bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.title.x=element_blank())
p

dev.off()



###### Supp FIG 10B ##########

vafs <- read.csv2("/home/jmartin/Documents/articulo/scripts/inputs/input_Supp_FIG10B.txt", header =T, dec='.', sep='\t')


vafs$TIPO <- factor( as.character(vafs$TIPO), levels=c("TISSUE", "PLASMA-BL", "PLASMA-PO", "PLASMA-RE") )

vafs$T_FRACTION = vafs$T_FRACTION*100

svg(file="/home/jmartin/Documents/articulo/scripts/figuras/Supp_FIG10B.svg", width = 13, height = 10)
p <- ggplot(vafs, aes(x=TIPO, y=T_FRACTION, color=INS, fill=INS)) + 
  geom_boxplot() + theme_bw() + scale_colour_manual(values =c("grey48","grey50")) + ylim(0, 100) +
  scale_fill_manual(values =c("#3C5488FF","#E64B35FF")) +
  ylab("Tumor Fraction") + 
  stat_compare_means(aes(group = INS), label = "p.signif") +
  theme(plot.title = element_text(size = 20, face = "bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=20))
p

dev.off()




###### Supp FIG 11C ##########

vafs <- read.csv2("/home/jmartin/Documents/articulo/scripts/inputs/input_Supp_FIG10C.txt", header =T, dec='.', sep=';')



vafs$TIPO <- factor( as.character(vafs$TIPO), levels=c("Tissue-Baseline", "Plasma-Baseline", "Plasma-PO", "Plasma-Relapse") )



svg(file="/home/jmartin/Documents/articulo/scripts/figuras/Supp_FIG10C.svg", width = 14, height = 10)
p <- ggplot(vafs, aes(x=TIPO, y=TMB, color=INS, fill=INS)) + 
  geom_boxplot() + theme_bw() + scale_colour_manual(values =c("grey48","grey50")) +
  scale_fill_manual(values =c("#3C5488FF","#E64B35FF")) +
  ylab("TMB") + 
  stat_compare_means(aes(group = INS), label = "p.signif") +
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=20))
p

dev.off()



