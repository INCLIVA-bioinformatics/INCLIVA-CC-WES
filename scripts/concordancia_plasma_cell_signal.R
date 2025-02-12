library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggrepel)

#juntar el df infiltration por la columna Mixture con el df porcentaje_exlusivo_PLASMA.BL_PLASMA por la columna ID.

porcentaje_exclusivo_PLASMA.BL_PLASMA = read.csv("/home/jmartin/Documents/articulo/scripts/inputs/Ext_Data_FIG6/input_Ext_Data_Fig6C1.txt", header = TRUE, sep = "\t", encoding = "ASCII", stringsAsFactors = FALSE)
infiltrarion = read.csv("/home/jmartin/Documents/articulo/scripts/inputs/Ext_Data_FIG6/input_Ext_Data_Fig6C2.csv", header = TRUE, sep = ";", encoding = "ASCII", stringsAsFactors = FALSE)

porcentaje_exclusivo_PLASMA.BL_PLASMA$ID <- paste0("EP", porcentaje_exclusivo_PLASMA.BL_PLASMA$ID)

merged_df <- merge(infiltrarion, porcentaje_exclusivo_PLASMA.BL_PLASMA, by.x = "Mixture", by.y = "ID")
colnames(merged_df)
rownames(merged_df) <- merged_df$Mixture
merged_df$Mixture <- NULL

correlation_matrix <- cor(merged_df)
correlation_long <- reshape2::melt(correlation_matrix)


cor_result <- cor.test(merged_df$CONCORDANCIA, merged_df$Plasma.cells)


svg(file="/home/jmartin/Documents/articulo/scripts/figuras/Ext_Data_Fig6C.svg", width = 10, height = 8)
ggplot(merged_df, aes(x = CONCORDANCIA, y = Plasma.cells
)) +
  geom_point(size = 3, color = "grey20", shape = 20, alpha = 0.7, aes(color = cor_result$p.value < 0.05)) +
  labs(x = "Mutational Concordance (Baseline Vs Relapse)", y = "Activated B cells") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(size = 20, face = "bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=20),
        axis.text=element_text(size=12),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20)) +  
  geom_smooth(method = "lm", se = TRUE, color = "black", alpha = 0.2) + 
  geom_text(x = max(merged_df$CONCORDANCIA), y = min(merged_df$Plasma.cells),
            label = "p-value = 0.0366",
            hjust = 1, vjust = 2, color = "black", size = 4)
dev.off()



