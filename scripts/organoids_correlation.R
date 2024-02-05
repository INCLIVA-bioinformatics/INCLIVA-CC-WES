library(cluster) 
library(ggplot2)
library(ggdendro)
library(tidyverse)
library(ComplexHeatmap)
library(reshape2)

############ SUPP FIG 10A #########

data <- read.delim("/home/jmartin/Documents/articulo/scripts/inputs/input_Supp_FIG10A.tsv")
data <- data %>% remove_rownames %>% column_to_rownames(var="X")

data2 = data[!(row.names(data) %in% "CTO51"),]

dd <- dist(data2, method = "binary")
hc <- hclust(dd, method = "ward.D2")

svg(file="/home/jmartin/Documents/articulo/scripts/figuras/Supp_FIG10A.svg", width = 15, height = 10)
ggdendrogram(hc, rotate = TRUE, theme_dendro = T)
dev.off()


############ SUPP FIG 10B #########

data = read.csv("/home/jmartin/Documents/articulo/scripts/inputs/input_Supp_FIG10B.txt", header = TRUE, sep = "\t", encoding = "ASCII", stringsAsFactors = FALSE)

data$MUESTRA <- factor( as.character(data$MUESTRA), levels=c("161-CTO65","49-CTO147","185-CTO119") )

m = acast(data, GEN~MUESTRA, value.var="TYPE")
m[is.na(m)] = ""

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  B = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#E64B35FF", col = NA))
  },
  RB = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "darkolivegreen4", col = NA))
  },
  R = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#3C5488FF", col = NA))
  }
)


col = c("R" = "#3C5488FF", "RB" = "darkolivegreen4", "B" = "#E64B35FF")


heatmap_legend_param = list(title = "Mutation stage", at = c("R", "RB", "B"), 
                            labels = c("ONLY IN THE PATIENT", "IN THE PATIENT AND IN THE ORGANOID", "ONLY IN THE ORGANOID"), direction = "horizontal")
column_title = "TARGETEABLE MUTATIONS CC PATIENTS - ORGANOID MODELS"

svg(file="/home/jmartin/Documents/articulo/scripts/figuras/Supp_FIG10B.svg", width = 8, height = 9)

plot_t =  oncoPrint(m,
                    alter_fun = alter_fun, col = col, 
                    remove_empty_columns = FALSE, remove_empty_rows = TRUE, column_title_gp = gpar(fontsize = 22),
                    pct_side = "right", row_names_side = "left", show_column_names = TRUE, column_names_gp = gpar(fontsize = 10), 
                    row_names_gp = gpar(fontsize = 10), column_order = c("161-CTO65","49-CTO147","185-CTO119"),
                    top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot()),
                    heatmap_legend_param = heatmap_legend_param) 
draw(plot_t)

dev.off()

