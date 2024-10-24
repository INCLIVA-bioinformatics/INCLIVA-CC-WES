
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



################ SUPP FIG 6C ############

data = read.csv("/home/jmartin/Documents/articulo/scripts/inputs/input_Supp_FIG6C.txt", header = TRUE, sep = "\t", encoding = "ASCII")

# Cargar las librerías necesarias
library(reshape2)


# Convertir ID a carácter para asegurarse de que se trata como categórico
data$ID <- as.character(data$ID)

# Transformar los datos a un formato largo para ggplot
data_long <- melt(data, id.vars = "ID", variable.name = "Group", value.name = "Y")

# Añadir una columna para clasificar los valores en "Epithelial" y "Mesenchymal"
data_long$Classification <- ifelse(data_long$Y > 0, "Epithelial", "Mesenchymal")

# Crear el gráfico
svg(file="//home/jmartin/Documents/articulo/scripts/repo_actualizado/INCLIVA-CC-WES/figures/Supp_FIG6C.svg", width = 12, height = 6)
ggplot(data_long, aes(x = ID, y = Y, color = Group)) +
  geom_point() +  # Dibujar puntos
  geom_hline(yintercept = 0, linetype = "dashed") +  # Añadir la línea horizontal en Y=0
  labs(x = "ID", y = "Value") +  # Añadir títulos y etiquetas
  scale_color_manual(values = c("Primary" = "#3C5488FF", "Metastasis" = "#E64B35FF")) +  # Especificar los colores
  theme_minimal() + ylab("Primary vs Metastasis EMT score") + # Usar un tema minimalista 
  theme(legend.title=element_blank(), legend.position = "top",
        legend.text=element_text(size=15),
        axis.text=element_text(size=12),
        axis.text.x=element_text(size=12, angle = 90, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=16)) +  # Rotar etiquetas del eje X
  annotate("text", x = Inf, y = 0.03, label = "Epithelial status", hjust = -0.1, vjust = -0.1, angle = 90, color = "black", size = 4.5) +  # Añadir etiqueta "Epithelial status"
  annotate("text", x = Inf, y = -0.33, label = "Mesenchymal status", hjust = -0.1, vjust = -0.1, angle = 90, color = "black", size = 4.5)  # Añadir etiqueta "Mesenchymal status"
dev.off()
