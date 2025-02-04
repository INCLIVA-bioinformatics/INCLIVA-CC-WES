
################# INCLIVA ###############3


de <- read.csv2("/home/jmartin/Documents/articulo/scripts/inputs/input_FIG5B1.csv", header =T, dec='.', sep='\t')


# add a column of NAs
de$diffexpressed <- "-"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$LOG2 > 0.6 & de$pvalue < 0.05] <- "Enriched_Relapse"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$LOG2 < -0.6 & de$pvalue < 0.05] <- "Enriched_Baseline"


# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=de, aes(x=LOG2, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("Enriched_Baseline", "Enriched_Relapse", "-")
p3 <- p2 + scale_colour_manual(values = mycolors)


de$delabel <- NA
de$delabel[de$diffexpressed != "-"] <- de$GENE[de$diffexpressed != "-"]

# Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
# load library
library(ggrepel)
# plot adding up all layers we have seen so far
svg(file="/home/jmartin/Documents/articulo/scripts/figuras/FIG4B1.svg", width = 18, height = 13)
ggplot(data=de, aes(x=LOG2, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = 37) + xlab(expression(log[2]*"FC")) + ylab(expression(-log[10]*"(pvalue)")) + ylim(0, 3.1) + xlim(-6,6) +
  scale_color_manual(values=c("black", "#3C5488FF", "#E64B35FF")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.title = element_blank(),
        legend.text=element_text(size=25),
        axis.text=element_text(size=15),
        axis.title=element_text(size=25))
dev.off()



################# MOMA ##################

de <- read.csv2("/home/jmartin/Documents/articulo/scripts/inputs/input_FIG5B2.csv", header =T, dec='.', sep='\t')


# add a column of NAs
de$diffexpressed <- "-"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$LOG2 > 0.6 & de$pvalue < 0.05] <- "Enriched_Relapse"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$LOG2 < -0.6 & de$pvalue < 0.05] <- "Enriched_Baseline"


# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=de, aes(x=LOG2, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("black", "steelblue", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("black", "steelblue", "red")
names(mycolors) <- c("Enriched_Baseline", "Enriched_Relapse", "-")
p3 <- p2 + scale_colour_manual(values = mycolors)


de$delabel <- NA
de$delabel[de$diffexpressed != "-"] <- de$GENE[de$diffexpressed != "-"]

# Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
# load library
library(ggrepel)
# plot adding up all layers we have seen so far
svg(file="/home/jmartin/Documents/articulo/scripts/figuras/FIG5B2.svg", width = 18, height = 13)
ggplot(data=de, aes(x=LOG2, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = 37) + xlab(expression(log[2]*"FC")) + ylab(expression(-log[10]*"(pvalue)")) + ylim(0, 3.1) + xlim(-6,6) +
  scale_color_manual(values=c("black", "#3C5488FF", "#E64B35FF")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.title = element_blank(),
        legend.text=element_text(size=25),
        axis.text=element_text(size=15),
        axis.title=element_text(size=25))
dev.off()
