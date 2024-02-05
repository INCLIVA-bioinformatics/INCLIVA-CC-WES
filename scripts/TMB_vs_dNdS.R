library(ggplot2)
library(vcd)
library(openintro)
library(tidyverse)
library(effsize)
library(ggplot2)
library(broom)
library(reshape2)

####### FIGURE 4A ##########

######## DISCOVERY COHORT #####

data = read.csv("/home/jmartin/Documents/articulo/scripts/inputs/input_FIG4A1.csv", header = TRUE, sep = "\t", encoding = "ASCII")

svg(file="/home/jmartin/Documents/articulo/scripts/figuras/FIG4A1.svg", width = 15, height = 10)
ggplot(data, aes(x=TMB, y=dN.dS, shape=SAMPLE_TYPE, color = SAMPLE_TYPE))  +
  geom_point() +
  theme(legend.position="top")  + geom_smooth(method = lm, fill = "grey80") + 
  scale_color_manual(values=c('#3C5488FF','#E64B35FF')) + annotate("text", x = 40, y = 1.04, color = "#3C5488FF",
                                                                     label = "paste(\"rho = 0.07\", \" p-value = 0.834\")", parse = TRUE) +  theme_bw() + ylab("dN/dS") +
  annotate("text", x = 40, y = 1.06, color = "#E64B35FF", label = "paste(\"rho = 0.454\", \" p-value = 0.0228\")", parse = TRUE) + 
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"), legend.position = "top")
dev.off()


######## VALIDATION COHORT #####

data = read.csv("/home/jmartin/Documents/articulo/scripts/inputs/input_FIG4A2.csv", header = TRUE, sep = ";", encoding = "ASCII")

plasma = cor.test(data$TMB_PLASMA, data$dnds_PLASMA, method = "spearman")
plasma1 = cor.test(data$TMB_PLASMA1, data$dnds_PLASMA1, method = "spearman")

pvalue_p = tidy(plasma)$p.value
pvalue_p1 = tidy(plasma1)$p.value

rhop = tidy(plasma)$estimate
rhop1 = tidy(plasma1)$estimate

value = data$TMB_PLASMA
value2 = data$dnds_PLASMA
type = rep("PLASMA-BL", times=nrow(data))
plasma = as.data.frame(cbind(type,value,value2))
colnames(plasma) = c("SAMPLE_TYPE", "TMB", "dN.dS")


value = data$TMB_PLASMA1
value2 = data$dnds_PLASMA1
type = rep("PLASMA-RE", times=nrow(data))
plasma1 = as.data.frame(cbind(type,value,value2))
colnames(plasma1) = c("SAMPLE_TYPE", "TMB", "dN.dS")

data2 = rbind(plasma, plasma1)
data2$TMB = as.numeric(as.character(data2$TMB))
data2$dN.dS = as.numeric(as.character(data2$dN.dS))

##### PLOTS ####

svg(file="/home/jmartin/Documents/articulo/scripts/figuras/FIG4A2.svg", width = 15, height = 10)
ggplot(data2, aes(x=TMB, y=dN.dS, shape=SAMPLE_TYPE, color = SAMPLE_TYPE))  +
  geom_point() +
  theme(legend.position="top")  + geom_smooth(method = lm, fill = "grey80") + 
  scale_color_manual(values=c('#3C5488FF','#E64B35FF')) + annotate("text", x = 11, y = 1.06, color = "orangered3",
                                                                     label = "paste(\"rho = 0.5541\", \" p-value = 0.0321\")", parse = TRUE) +  theme_bw() + ylab("dN/dS") +
  annotate("text", x = 11, y = 1.05, color = "#3C5488FF", label = "paste(\"rho = 0.3107\", \" p-value = 0.2592\")", parse = TRUE) + 
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"), legend.position = "top")
dev.off()

######### SUPP FIGURE 8 #########

#################### SUPP FIGURE 8A ########

tmb_plasma <- read_delim("/home/jmartin/Documents/articulo/scripts/inputs/input_Supp_FIG8A.csv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

colnames(tmb_tejido) = c("NHC", "Tissue Relapse", "Primary Tissue")

wilcox.test(tmb_tejido$`Tissue Relapse`,tmb_tejido$`Primary Tissue`, paired = TRUE)

tmb_tejido_long <- gather(tmb_tejido, condition, value,`Tissue Relapse`:`Primary Tissue`, factor_key=TRUE)
tmb_tejido_long$condition <- car::recode(tmb_tejido_long$condition,
                                         "'TPB-PORCENTAJE'='TPB';
                                         'MP-PORCENTAGE'='MP'")

svg(file="/home/jmartin/Documents/articulo/scripts/figuras/Supp_FIG8B.svg", width = 10, height = 8)
p <- ggpaired(tmb_tejido_long, x = "condition", y = "value",
              color = "condition", line.color = "gray", line.size = 0.4)+
  scale_color_manual(values=c("#3C5488FF", "#E64B35FF")) +
  stat_compare_means(paired = TRUE, method = "wilcox.test",
                     label.x = 1.4,
                     label.y = 60) + ylab("TMB") +
  theme(plot.title = element_text(size = 15, face = "bold"),
        legend.position = "none",
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))
p
dev.off()


#################### SUPP FIGURE 8B ############



tmb_tejido <- read_delim("/home/jmartin/Documents/articulo/scripts/inputs/input_Supp_FIG8B.csv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

colnames(tmb_plasma) = c("NHC", "Plasma Relapse", "Plasma Baseline")
wilcox.test(tmb_plasma$`Plasma Relapse`,tmb_plasma$`Plasma Baseline`, paired = TRUE)

tmb_plasma_long <- gather(tmb_plasma, condition, value,`Plasma Relapse`:`Plasma Baseline`, factor_key=TRUE)
tmb_plasma_long$condition <- car::recode(tmb_plasma_long$condition,
                                         "'TPB-PORCENTAJE'='TPB';
                                         'MP-PORCENTAGE'='MP'")

svg(file="/home/jmartin/Documents/articulo/scripts/figuras/Supp_FIG8A.svg", width = 10, height = 8)
p <- ggpaired(tmb_plasma_long, x = "condition", y = "value",
              color = "condition", line.color = "gray", line.size = 0.4)+
  scale_color_manual(values=c("#3C5488FF", "#E64B35FF")) +
  stat_compare_means(paired = TRUE, method = "wilcox.test",
                     label.x = 1.4,
                     label.y = 50) + ylab("TMB") +
  theme(plot.title = element_text(size = 15, face = "bold"),
        legend.position = "none",
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))
p
dev.off()


#################### SUPP FIGURE 8C ########

tmb_plasma <- read_delim("/home/jmartin/Documents/articulo/scripts/inputs/input_Supp_FIG8C.csv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

colnames(tmb_plasma) = c("NHC", "Plasma Relapse", "Plasma Baseline")

wilcox.test(tmb_plasma$`Plasma Relapse`,tmb_plasma$`Plasma Baseline`, paired = TRUE)

tmb_plasma_long <- gather(tmb_plasma, condition, value,`Plasma Relapse`:`Plasma Baseline`, factor_key=TRUE)
tmb_plasma_long$condition <- car::recode(tmb_plasma_long$condition,
                                         "'TPB-PORCENTAJE'='TPB';
                                         'MP-PORCENTAGE'='MP'")

svg(file="/home/jmartin/Documents/articulo/scripts/figuras/Supp_FIG8C.svg", width = 10, height = 8)
p <- ggpaired(tmb_plasma_long, x = "condition", y = "value",
              color = "condition", line.color = "gray", line.size = 0.4)+
  scale_color_manual(values=c("#3C5488FF", "#E64B35FF")) +
  stat_compare_means(paired = TRUE, method = "wilcox.test",
                     label.x = 1.4,
                     label.y = 12) + ylab("TMB") +
  theme(plot.title = element_text(size = 15, face = "bold"),
        legend.position = "none",
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))
p
dev.off()



######## SUPP FIGURE 9 ######

data = read.csv("/home/jmartin/Documents/articulo/scripts/inputs/input_Supp_FIG9.csv", header = TRUE, sep = "\t", encoding = "ASCII")

svg(file="/home/jmartin/Documents/articulo/scripts/figuras/Supp_FIG9.svg", width = 15, height = 10)
ggplot(data, aes(x=TMB, y=dN.dS, shape=SAMPLE_TYPE, color = SAMPLE_TYPE))  +
  geom_point() +
  theme(legend.position="top")  + geom_smooth(method = lm, fill = "grey80") + 
  scale_color_manual(values=c('#3C5488FF','#E64B35FF')) + annotate("text", x = 50, y = 1.25, color = "#3C5488FF",
                                                                     label = "paste(\"rho = 0.393\", \" p-value = 0.0785\")", parse = TRUE) +
  annotate("text", x = 50, y = 1.22, color = "#E64B35FF", label = "paste(\"rho = 0.4818\", \" p-value = 0.0199\")", parse = TRUE) + theme_bw() + ylab("dN/dS") +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"), legend.position = "top")
dev.off()
