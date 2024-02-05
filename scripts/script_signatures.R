library(deconstructSigs)
library("BSgenome.Hsapiens.UCSC.hg38")
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(patchwork)
library(gridExtra)
library(cowplot)

#################### SUPP FIGURE 7A ###############


###################### INPUT DATA   ######################

signatures = read.csv("/home/jmartin/Documents/articulo/scripts/inputs/input_Supp_FIG7A.csv", header = TRUE, sep = "\t", encoding = "ASCII", stringsAsFactors = FALSE)

clinica = read.csv("/home/jmartin/Documents/articulo/scripts/inputs/data_clinica_25.csv", header = TRUE, sep = "{", encoding = "ASCII", stringsAsFactors = FALSE)

signatures$Samples = factor(signatures$Samples, levels = c("EP107", "EP13", "EP136", "EP138", "EP150", "EP24", "EP45", "EP48", "EP49", "EP62", "EP63", "EP7", "EP8", "EP104",  "EP161", "EP185", "EP189", "EP204", "EP207", "EP242", "EP243", "EP158", "EP85", "EP259","EP261"))
clinica$Samples = factor(clinica$Samples, levels = c("EP107", "EP13", "EP136", "EP138", "EP150", "EP24", "EP45", "EP48", "EP49", "EP62", "EP63", "EP7", "EP8", "EP104",  "EP161", "EP185", "EP189", "EP204", "EP207", "EP242", "EP243","EP158", "EP85", "EP259","EP261"))


########## PLOT ##########


# Stacked + percent
#tiff(file="/home/jmartin/Documentos/exomas/nuevos_exomas/TANDA3/signatures/SigProfiler/COSMIC_SBS96_percentages.tiff", width = 50, height = 15,units = 'cm', res = 600)
p = ggplot(signatures, aes(fill=Signature, y=Percentage, x=Samples)) +
  scale_fill_manual(values = c("palegreen",  "khaki2", "steelblue",  "lemonchiffon3",  "wheat3", "seagreen", "mediumorchid1"), labels= c("SBS1 (DNA mismatch)", "SBS15 (MSI)","SBS23 (Unknown)", "SBS43 (Sequencing artifacts)", "SBS46 (Sequencing artifacts)", "SBS5 (ERCC2 mutations)",  "SBS96B")) +
  geom_bar(position="fill", stat="identity") + theme_bw() + ggtitle("Plasma at Relapse Signatures Distribution (SBS)") +
  xlab("Sample") + ylab("Percentage") +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14),
        axis.title=element_text(size=20,face="bold"),
        axis.title.x = element_blank())
p
#dev.off()



age = ggplot(clinica, aes(Samples, AGE)) + geom_point() + theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ,axis.title.x = element_blank())

tto = ggplot(clinica, aes(fill=TTO.Y.N., y=Num, x=Samples))  +
  scale_fill_manual(values = c("snow3", "cadetblue")) +
  geom_bar(position="fill", stat="identity") + theme_void() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())

msi = ggplot(clinica, aes(fill=MSI, y=Num, x=Samples))  +
  scale_fill_manual(values = c( "snow3", "mediumorchid")) +
  geom_bar(position="fill", stat="identity") + theme_void() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())

gender_p = ggplot(clinica, aes(fill=Gender, y=Num, x=Samples))  +
  geom_bar(position="fill", stat="identity") + theme_void() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())

stage = ggplot(clinica, aes(fill=StStage, y=Num, x=Samples))  +
  geom_bar(position="fill", stat="identity") + theme_void() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())

site = ggplot(clinica, aes(fill=Location, y=Num, x=Samples))  +
  scale_fill_manual(values = c("darkolivegreen1",  "steelblue1")) +
  geom_bar(position="fill", stat="identity") + theme_void() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())

batch = ggplot(clinica, aes(fill=BATCH, y=Num, x=Samples))  +
  scale_fill_manual(values = c("black",  "grey", "grey95")) +
  geom_bar(position="fill", stat="identity") + theme_void() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())


svg(file="/home/jmartin/Documents/articulo/scripts/figuras/Supp_Fig7A.svg", width = 28, height = 15)
plot_grid(age + theme(legend.position = "left", legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=8)),
          tto+ theme( legend.text = element_text(size=8), legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=8)),
          batch+ theme(legend.position = "left", legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=8)), 
          msi+ theme(legend.text = element_text(size=8), legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=8)), 
          site+ theme(legend.position = "left", legend.text = element_text(size=8), legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=8)), 
          stage+ theme( legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=8), legend.direction = "horizontal"),
          gender_p+ theme( legend.text = element_text(size=8), legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=8), legend.position = "left"),
          p, 
          align = "v", 
          axis = "lr",
          nrow = 8, 
          rel_heights = c(3/80, 3/80, 3/80, 3/80,3/80,3/80,3/80, 7/10))
dev.off()



#################### SUPP FIGURE 7b ###############


###################### INPUT DATA ######################

signatures = read.csv("/home/jmartin/Documents/articulo/scripts/inputs/input_Supp_FIG7B1.csv", header = TRUE, sep = "\t", encoding = "ASCII", stringsAsFactors = FALSE)

signatures2 = read.csv("/home/jmartin/Documents/articulo/scripts/inputs/input_Supp_FIG7A.csv", header = TRUE, sep = "\t", encoding = "ASCII", stringsAsFactors = FALSE)

signatures3 = signatures2[signatures2$Samples %in% c("EP107", "EP138", "EP8",  "EP161", "EP185", "EP204", "EP207", "EP242", "EP243", "EP158","EP85","EP261"),]

clinica = read.csv("/home/jmartin/Documents/articulo/scripts/inputs/input_Supp_FIG672.csv", header = TRUE, sep = "{", encoding = "ASCII", stringsAsFactors = FALSE)

signatures$Samples = factor(signatures$Samples, levels = c("EPDx107", "EPDx138", "EPDx8",  "EPDx161", "EPDx185", "EPDx204", "EPDx207", "EPDx242", "EPDx243", "EPDx158","EPDx85","EPDx261"))
signatures3$Samples = factor(signatures3$Samples, levels = c("EP107", "EP138", "EP8",  "EP161", "EP185", "EP204", "EP207", "EP242", "EP243", "EP158","EP85","EP261"))
clinica$Samples = factor(clinica$Samples, levels = c("EP107", "EP138", "EP8",  "EP161", "EP185", "EP204", "EP207", "EP242", "EP243", "EP158","EP85","EP261"))


########## PLOT ##########


# Stacked + percent
#tiff(file="/home/jmartin/Documentos/exomas/nuevos_exomas/TANDA3/signatures/SigProfiler/COSMIC_SBS96_percentages.tiff", width = 50, height = 15,units = 'cm', res = 600)
p = ggplot(signatures, aes(fill=Signature, y=Percentage, x=Samples)) +
  scale_fill_manual(values = c("palegreen",  "khaki2",  "lemonchiffon3",  "wheat3", "seagreen", "mediumorchid1"), labels= c("SBS1 (DNA mismatch)", "SBS15 (MSI)", "SBS43 (Sequencing artifacts)", "SBS46 (Sequencing artifacts)", "SBS5 (ERCC2 mutations)")) +
  geom_bar(position="fill", stat="identity") + theme_bw() + ggtitle("Plasma at Baseline Signatures Distribution (SBS)") + ylab("Percentage") +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14),
        axis.title=element_text(size=20,face="bold"),
        axis.title.x = element_blank())

p2 = ggplot(signatures3, aes(fill=Signature, y=Percentage, x=Samples)) +
  scale_fill_manual(values = c("palegreen",  "khaki2", "steelblue",  "lemonchiffon3",  "wheat3", "seagreen", "mediumorchid1"), labels= c("SBS1 (DNA mismatch)", "SBS15 (MSI)","SBS23 (Unknown)", "SBS43 (Sequencing artifacts)", "SBS46 (Sequencing artifacts)", "SBS5 (ERCC2 mutations)",  "SBS96B")) +
  geom_bar(position="fill", stat="identity") + theme_bw() + ggtitle("Plasma at Relapse Signatures Distribution (SBS)") +
  xlab("Sample") + ylab("Percentage") +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14),
        axis.title=element_text(size=20,face="bold"),
        axis.title.x = element_blank())
#dev.off()



age = ggplot(clinica, aes(Samples, AGE)) + geom_point() + theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ,axis.title.x = element_blank())

tto = ggplot(clinica, aes(fill=TTO.Y.N., y=Num, x=Samples))  +
  scale_fill_manual(values = c("snow3", "cadetblue")) +
  geom_bar(position="fill", stat="identity") + theme_void() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())

msi = ggplot(clinica, aes(fill=MSI, y=Num, x=Samples))  +
  scale_fill_manual(values = c( "snow3", "mediumorchid")) +
  geom_bar(position="fill", stat="identity") + theme_void() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())

gender_p = ggplot(clinica, aes(fill=Gender, y=Num, x=Samples))  +
  geom_bar(position="fill", stat="identity") + theme_void() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())

stage = ggplot(clinica, aes(fill=StStage, y=Num, x=Samples))  +
  geom_bar(position="fill", stat="identity") + theme_void() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())

site = ggplot(clinica, aes(fill=Location, y=Num, x=Samples))  +
  scale_fill_manual(values = c("darkolivegreen1",  "steelblue1")) +
  geom_bar(position="fill", stat="identity") + theme_void() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())

batch = ggplot(clinica, aes(fill=BATCH, y=Num, x=Samples))  +
  scale_fill_manual(values = c("black",  "grey","grey95")) +
  geom_bar(position="fill", stat="identity") + theme_void() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())


svg(file="/home/jmartin/Documents/articulo/scripts/figuras/Supp_Fig7B.svg", width = 18, height = 15)
plot_grid(age + theme(legend.position = "left", legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=8)),
          tto+ theme( legend.text = element_text(size=8), legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=8)),
          batch+ theme(legend.position = "left", legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=8)), 
          msi+ theme(legend.text = element_text(size=8), legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=8)), 
          site+ theme(legend.position = "left", legend.text = element_text(size=8), legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=8)), 
          stage+ theme( legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=8), legend.direction = "horizontal"),
          gender_p+ theme( legend.text = element_text(size=8), legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=8), legend.position = "left"),
          p, p2,
          align = "v", 
          axis = "lr",
          nrow = 9, 
          rel_heights = c(3/80, 3/80, 3/80, 3/80,3/80,3/80,3/80, 7/20, 7/20))
dev.off()