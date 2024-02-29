### Script para generar los diagramas de Venn de un listado de VCFs ###
# jfcatala@incliva.es #
#   20230315 - 1355   #

# Notas: Por ahora asumo que los VCFs son monomuestra

library(vcfR)
library(tidyverse)
library(RColorBrewer)
library(UpSetR)
library(ggVennDiagram)


extract_variants_from_vcf <- function(vcf_file, ID, type = "all", filter = "PASS"){
  #Falta filtro de separación del indel-SNV
  #ID = "full"  <- chr_pos_ref_alt
  #ID = "small" <- chr_pos

  vcf <- read.vcfR( vcf_file, verbose = FALSE ) #Leo el fichero VCF

  vcf_fix <- as_tibble(getFIX(vcf)) %>% #Obtenfo la parte "fija" del VCF
    tidyr::separate_rows(ALT, sep = ",") %>% #Separo variantes multialélicas
    mutate(type = if_else(nchar(REF) == nchar(ALT), "SNV", "Indel")) #Defino SNVs e Indels

  vcf_fix <- vcf_fix %>% filter(FILTER %in% c(filter,NA))

  if(type == "SNV"){
    vcf_fix <- vcf_fix %>% filter(type %in% c("SNV"))
  } else if (type == "Indel"){
    vcf_fix <- vcf_fix %>% filter(type %in% c("Indel"))
  } else {
    vcf_fix <- vcf_fix
  }

  if(ID == "full"){
    variants <- str_c(vcf_fix$CHROM, vcf_fix$POS, vcf_fix$REF, vcf_fix$ALT, sep = "_")
  } else if (ID == "small"){
    variants <- str_c(vcf_fix$CHROM, vcf_fix$POS, sep = "_")
  } else {
    stop()
  }
  return(variants)
  
}

intersect_variants <- function(list, pl_type = "ggvenn", title = ""){
  if (pl_type == "upset"){
    # No he encontrado forma fácil de poner título a los Upset
    g <- upset(fromList(list), order.by = "degree", nsets = length(list))
  } else if (pl_type == "ggvenndiagram"){
    g <- ggVennDiagram(list) + 
      scale_fill_distiller(palette = "Reds", trans = "reverse") +
      ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, vjust = 4.5, size = 25))
      
  }
  return(g)
}

# List of VCFs | Futurible función para argparse
#list_of_vcf = [name1 = path1, name2 = path2]
list_of_vcf = list(tejido = "/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/ET-M2_FFPE_somatic.vcf.gz",
                plasma_bl = "/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/EPDx-M2_plasma_somatic.vcf.gz",
                plasma_po = "/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/EP-PO-M2_plasma_somatic.vcf.gz",
                plasma_re = "/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/EP-M2_plasma_somatic.vcf.gz" )

list_of_variants <- lapply(list_of_vcf, function(x){(extract_variants_from_vcf(x, ID = "full", type = "all"))})

intersect_variants(list_of_variants, pl_type = "upset")


#############3 MOMA ################



varst <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/tissue_vars.txt", header =F, dec='.', sep='\t')
varsvt = as.character(varst$V1)


varspb <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/plasma-bl_vars.txt", header =F, dec='.', sep='\t')
varsvpb = as.character(varspb$V1)

varspo <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/plasma-po_vars.txt", header =F, dec='.', sep='\t')
varsvpo = as.character(varspo$V1)

varsp <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/plasma_vars.txt", header =F, dec='.', sep='\t')
varsvp = as.character(varsp$V1)

list_of_variants = list(tissue = varsvt, plasma_bl = varsvpb, plasma_po = varsvpo, plasma_re = varsvp)

upset(fromList(list_of_variants), order.by = "degree", nsets = length(list_of_variants))






varst <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/tissue_genes.txt", header =F, dec='.', sep='\t')
varsvt = as.character(varst$V1)


varspb <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/plasma-bl_genes.txt", header =F, dec='.', sep='\t')
varsvpb = as.character(varspb$V1)

varspo <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/plasma-po_genes.txt", header =F, dec='.', sep='\t')
varsvpo = as.character(varspo$V1)

varsp <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/plasma_genes.txt", header =F, dec='.', sep='\t')
varsvp = as.character(varsp$V1)

list_of_variants = list(tissue = varsvt, plasma_bl = varsvpb, plasma_po = varsvpo, plasma_re = varsvp)

upset(fromList(list_of_variants), order.by = "degree", nsets = length(list_of_variants))



pre = c(varsvt, varsvpb, varsvpo)
diff_plasma = setdiff(varsvp, pre)

write.csv(diff_plasma, "/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/plasma-exclusive_genes.txt",col.names=F, row.names = FALSE, quote = FALSE)




############ JOIN ##########

overlapGroups <- function (listInput, sort = TRUE) {
  # listInput could look like this:
  # $one
  # [1] "a" "b" "c" "e" "g" "h" "k" "l" "m"
  # $two
  # [1] "a" "b" "d" "e" "j"
  # $three
  # [1] "a" "e" "f" "g" "h" "i" "j" "l" "m"
  listInputmat    <- fromList(listInput) == 1
  #     one   two three
  # a  TRUE  TRUE  TRUE
  # b  TRUE  TRUE FALSE
  #...
  # condensing matrix to unique combinations elements
  listInputunique <- unique(listInputmat)
  grouplist <- list()
  # going through all unique combinations and collect elements for each in a list
  for (i in 1:nrow(listInputunique)) {
    currentRow <- listInputunique[i,]
    myelements <- which(apply(listInputmat,1,function(x) all(x == currentRow)))
    attr(myelements, "groups") <- currentRow
    grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <- myelements
    myelements
    # attr(,"groups")
    #   one   two three 
    # FALSE FALSE  TRUE 
    #  f  i 
    # 12 13 
  }
  if (sort) {
    grouplist <- grouplist[order(sapply(grouplist, function(x) length(x)), decreasing = TRUE)]
  }
  attr(grouplist, "elements") <- unique(unlist(listInput))
  return(grouplist)
  # save element list to facilitate access using an index in case rownames are not named
}

li = overlapGroups(list_of_variants)

ta = attr(li, "elements")

li[1]
tissue = attr(li, "elements")[li[[1]]]



d = as.data.frame(lengths(li))
library(tibble)
d <- tibble::rownames_to_column(d, "names")
dt= cbind(d, rep("Tissue", length(rownames(d))))
colnames(dt) = c("names", "len", "tipo")
dpb= cbind(d, rep("Plasma-BL", length(rownames(d))))
colnames(dpb) = c("names", "len", "tipo")
dpo= cbind(d, rep("Plasma-PO", length(rownames(d))))
colnames(dpo) = c("names", "len", "tipo")
dp= cbind(d, rep("Plasma-RE", length(rownames(d))))
colnames(dp) = c("names", "len", "tipo")

data = rbind(dt, dpb, dpo, dp)






dt$names = factor(dt$names, levels = rev(c("tissue", "tissue:plasma_bl", "tissue:plasma_bl:plasma_po",
                                           "tissue:plasma_po" , "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re", "tissue:plasma_po:plasma_re",
                                           "tissue:plasma_re","plasma_bl","plasma_bl:plasma_po", "plasma_bl:plasma_po:plasma_re",
                                           "plasma_bl:plasma_re", "plasma_po", "plasma_po:plasma_re", "plasma_re")))

dpb$names = factor(dpb$names, levels = rev(c("tissue", "tissue:plasma_bl", "tissue:plasma_bl:plasma_po",
                                             "tissue:plasma_po" , "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re", "tissue:plasma_po:plasma_re",
                                             "tissue:plasma_re","plasma_bl","plasma_bl:plasma_po", "plasma_bl:plasma_po:plasma_re",
                                             "plasma_bl:plasma_re", "plasma_po", "plasma_po:plasma_re", "plasma_re")))

dpo$names = factor(dpo$names, levels = rev(c("tissue", "tissue:plasma_bl", "tissue:plasma_bl:plasma_po",
                                             "tissue:plasma_po" , "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re", "tissue:plasma_po:plasma_re",
                                             "tissue:plasma_re","plasma_bl","plasma_bl:plasma_po", "plasma_bl:plasma_po:plasma_re",
                                             "plasma_bl:plasma_re", "plasma_po", "plasma_po:plasma_re", "plasma_re")))

dp$names = factor(dp$names, levels = rev(c("tissue", "tissue:plasma_bl", "tissue:plasma_bl:plasma_po",
                                           "tissue:plasma_po" , "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re", "tissue:plasma_po:plasma_re",
                                           "tissue:plasma_re","plasma_bl","plasma_bl:plasma_po", "plasma_bl:plasma_po:plasma_re",
                                           "plasma_bl:plasma_re", "plasma_po", "plasma_po:plasma_re", "plasma_re")))

#write.csv(data, "/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/data.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#data <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/data.txt", header =T, dec='.', sep=',')


###########3 PLOT ##############

library(ggpubr)
t = ggplot(dt,aes(x = tipo, y = len, fill = names))  +
  geom_col() +
  scale_fill_manual(values = ifelse(grepl("tissue", levels(dt$names), fixed = TRUE), "grey50", "white"))  + coord_flip() +
  theme(plot.title = element_blank(), legend.position = "none",
        axis.title=element_blank(),
        axis.text.x = element_blank())

pb =  ggplot(dpb,aes(x = tipo, y = len, fill = names))  +
  geom_col() +
  scale_fill_manual(values = ifelse(grepl("plasma_bl", levels(dpb$names), fixed = TRUE), "grey50", "white")) +coord_flip() +
  theme(plot.title = element_blank(), legend.position = "none",
        axis.title=element_blank(),
        axis.text.x = element_blank())

po = ggplot(dpo,aes(x = tipo, y = len, fill = names))  +
  geom_col() +
  scale_fill_manual(values = ifelse(grepl("plasma_po", levels(dpo$names), fixed = TRUE), "grey50", "white"))+coord_flip() +
  theme(plot.title = element_blank(), legend.position = "none",
        axis.title=element_blank(),
        axis.text.x = element_blank())
p =  ggplot(dp,aes(x = tipo, y = len, fill = names))  +
  geom_col() +
  scale_fill_manual(values = ifelse(grepl("plasma_re", levels(dp$names), fixed = TRUE), "grey50", "white"))+coord_flip() +
  theme(plot.title = element_blank(), legend.position = "none",
        axis.title=element_blank(),
        axis.text.x = element_blank())


library(dplyr)
library(RColorBrewer)
library(patchwork)
library(gridExtra)
library(cowplot)

tiff(file="/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots//enriquecimiento_MOMA_genes_total.tiff", width = 18, height = 8,units = 'cm', res = 400)
plot_grid(t, pb, po, p,  align = "v", 
          axis = "lr",
          nrow = 4)

dev.off()





############### FISH PLOT ###################
library("ggsci")
library("gridExtra")



write.csv(data, "/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/counts.txt",col.names=F, row.names = T, quote = FALSE)

data2 <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/counts_mod.txt", header =T, dec='.', sep=',')
data2$names = factor(data2$names, levels = rev(c("tissue", "plasma_bl","tissue:plasma_bl",  "tissue:plasma_po" , "tissue:plasma_bl:plasma_po",
                                                 "plasma_bl:plasma_po", "tissue:plasma_re",  "plasma_bl:plasma_re", 
                                                 "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re",
                                                 "tissue:plasma_po:plasma_re",
                                                 "plasma_bl:plasma_po:plasma_re", 
                                                 "plasma_po", 
                                                 "plasma_po:plasma_re", 
                                                 "plasma_re")))

data2$tipo = factor(data2$tipo, levels = c("Tissue", "Plasma-BL", "Plasma-PO","Plasma-RE"))

col <- c("#DC0000FF", "#E64B35FF",  "#DC0000B2",   "#3C5488B2", "#8491B4B2",   "#3C5488D8", "#3C5488FF", "#4DBBD5FF",  
         "#4DBBD5B2",  "#00A087FF",  "#91D1C2FF",  "#91D1C2B2",  "#E64B3533", "#B09C85B2", "#B09C85FF")


svg(file="/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/fish_plots/FISH_plot_MOMA1.svg", width = 10, height = 8)
ggplot(data2, aes(fill=names, y=len, x=tipo)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = col)
dev.off()



#############


data2 <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/counts_mod2.txt", header =T, dec='.', sep=',')
data2$names = factor(data2$names, levels = rev(c("tissue", "no_tissue", "plasma_bl","no_plasma_bl","tissue:plasma_bl", "no_tissue:plasma_bl", "tissue:plasma_po" , "no_tissue:plasma_po", "tissue:plasma_bl:plasma_po", "no_tissue:plasma_bl:plasma_po",
                                                 "plasma_bl:plasma_po", "no_plasma_bl:plasma_po","tissue:plasma_re", "no_tissue:plasma_re", "plasma_bl:plasma_re", "no_plasma_bl:plasma_re",
                                                 "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re",
                                                 "no_tissue:plasma_bl:plasma_re", "tissue:plasma_po:plasma_re", "no_tissue:plasma_po:plasma_re",
                                                 "plasma_bl:plasma_po:plasma_re", "no_plasma_bl:plasma_po:plasma_re",
                                                "plasma_po", "no_plasma_po",
                                                 "plasma_po:plasma_re", "no_plasma_po:plasma_re", "plasma_re", "no_plasma_re")))

data2$tipo = factor(data2$tipo, levels = c("Tissue", "Plasma-BL", "Plasma-PO","Plasma-RE"))

availableColours <- sample(col_vector, n)
legendColours <- ifelse(str_detect(data2$names, fixed("no")), "white", availableColours)
ifelse(str_detect(levels(data2$names), fixed("no")), "white", availableColours)

data3 = data2[order(data2$names),]

col <- c("white",  "#DC0000FF", "white",   "#E64B35FF", "white",   "#DC0000B2",  "white",   "#3C5488B2",
         "white",   "#8491B4B2", "white",   "#3C5488D8", "#3C5488FF","white",   "#4DBBD5FF","white",   
         "#4DBBD5B2", "white",   "#00A087FF",  "white",   "#91D1C2FF","white",   "#91D1C2B2", "white",   "#E64B3533","white","#B09C85B2", "white",   "#B09C85FF")

#names(col) <- as.character(data2$color)

svg(file="/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/fish_plots/FISH_plot_MOMA2.svg", width = 12, height = 8)
ggplot(data3, aes(fill=names, y=len, x=tipo)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = col)
dev.off()


write.csv(a, "/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/a.txt",col.names=F, row.names = F, quote = FALSE)







############# COLORES ################

library("scales")
pal_npg("nrc")(10)
pal_npg("nrc", alpha = 0.3)(10)










##########3 INCLIVA ###################



varst <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/tissue_vars.txt", header =F, dec='.', sep='\t')
varsvt = as.character(varst$V1)


varspb <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/plasma-bl_vars.txt", header =F, dec='.', sep='\t')
varsvpb = as.character(varspb$V1)

varspo <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/plasma-po_vars.txt", header =F, dec='.', sep='\t')
varsvpo = as.character(varspo$V1)

varsp <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/plasma_vars.txt", header =F, dec='.', sep='\t')
varsvp = as.character(varsp$V1)

list_of_variants = list(tissue = varsvt, plasma_bl = varsvpb, plasma_po = varsvpo, plasma_re = varsvp)

upset(fromList(list_of_variants), order.by = "degree", nsets = length(list_of_variants))






varst <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/tissue_genes.txt", header =F, dec='.', sep='\t')
varsvt = as.character(varst$V1)


varspb <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/plasma-bl_genes.txt", header =F, dec='.', sep='\t')
varsvpb = as.character(varspb$V1)

varspo <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/plasma-po_genes.txt", header =F, dec='.', sep='\t')
varsvpo = as.character(varspo$V1)

varsp <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/plasma_genes.txt", header =F, dec='.', sep='\t')
varsvp = as.character(varsp$V1)

list_of_variants = list(tissue = varsvt, plasma_bl = varsvpb, plasma_po = varsvpo, plasma_re = varsvp)

upset(fromList(list_of_variants), order.by = "degree", nsets = length(list_of_variants))



pre = c(varsvt, varsvpb, varsvpo)
diff_plasma = setdiff(varsvp, pre)

write.csv(diff_plasma, "/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/plasma-exclusive_genes.txt",col.names=F, row.names = FALSE, quote = FALSE)




############ JOIN ##########

overlapGroups <- function (listInput, sort = TRUE) {
  # listInput could look like this:
  # $one
  # [1] "a" "b" "c" "e" "g" "h" "k" "l" "m"
  # $two
  # [1] "a" "b" "d" "e" "j"
  # $three
  # [1] "a" "e" "f" "g" "h" "i" "j" "l" "m"
  listInputmat    <- fromList(listInput) == 1
  #     one   two three
  # a  TRUE  TRUE  TRUE
  # b  TRUE  TRUE FALSE
  #...
  # condensing matrix to unique combinations elements
  listInputunique <- unique(listInputmat)
  grouplist <- list()
  # going through all unique combinations and collect elements for each in a list
  for (i in 1:nrow(listInputunique)) {
    currentRow <- listInputunique[i,]
    myelements <- which(apply(listInputmat,1,function(x) all(x == currentRow)))
    attr(myelements, "groups") <- currentRow
    grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <- myelements
    myelements
    # attr(,"groups")
    #   one   two three 
    # FALSE FALSE  TRUE 
    #  f  i 
    # 12 13 
  }
  if (sort) {
    grouplist <- grouplist[order(sapply(grouplist, function(x) length(x)), decreasing = TRUE)]
  }
  attr(grouplist, "elements") <- unique(unlist(listInput))
  return(grouplist)
  # save element list to facilitate access using an index in case rownames are not named
}

li = overlapGroups(list_of_variants)

ta = attr(li, "elements")

li[1]
tissue = attr(li, "elements")[li[[1]]]



d = as.data.frame(lengths(li))
library(tibble)
d <- tibble::rownames_to_column(d, "names")
dt= cbind(d, rep("Tissue", length(rownames(d))))
colnames(dt) = c("names", "len", "tipo")
dpb= cbind(d, rep("Plasma-BL", length(rownames(d))))
colnames(dpb) = c("names", "len", "tipo")
dpo= cbind(d, rep("Plasma-PO", length(rownames(d))))
colnames(dpo) = c("names", "len", "tipo")
dp= cbind(d, rep("Plasma", length(rownames(d))))
colnames(dp) = c("names", "len", "tipo")

data = rbind(dt, dpb, dpo, dp)






dt$names = factor(dt$names, levels = rev(c("tissue", "tissue:plasma_bl", "tissue:plasma_bl:plasma_po",
                                           "tissue:plasma_po" , "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re", "tissue:plasma_po:plasma_re",
                                           "tissue:plasma_re","plasma_bl","plasma_bl:plasma_po", "plasma_bl:plasma_po:plasma_re",
                                           "plasma_bl:plasma_re", "plasma_po", "plasma_po:plasma_re", "plasma_re")))

dpb$names = factor(dpb$names, levels = rev(c("tissue", "tissue:plasma_bl", "tissue:plasma_bl:plasma_po",
                                             "tissue:plasma_po" , "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re", "tissue:plasma_po:plasma_re",
                                             "tissue:plasma_re","plasma_bl","plasma_bl:plasma_po", "plasma_bl:plasma_po:plasma_re",
                                             "plasma_bl:plasma_re", "plasma_po", "plasma_po:plasma_re", "plasma_re")))

dpo$names = factor(dpo$names, levels = rev(c("tissue", "tissue:plasma_bl", "tissue:plasma_bl:plasma_po",
                                             "tissue:plasma_po" , "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re", "tissue:plasma_po:plasma_re",
                                             "tissue:plasma_re","plasma_bl","plasma_bl:plasma_po", "plasma_bl:plasma_po:plasma_re",
                                             "plasma_bl:plasma_re", "plasma_po", "plasma_po:plasma_re", "plasma_re")))

dp$names = factor(dp$names, levels = rev(c("tissue", "tissue:plasma_bl", "tissue:plasma_bl:plasma_po",
                                           "tissue:plasma_po" , "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re", "tissue:plasma_po:plasma_re",
                                           "tissue:plasma_re","plasma_bl","plasma_bl:plasma_po", "plasma_bl:plasma_po:plasma_re",
                                           "plasma_bl:plasma_re", "plasma_po", "plasma_po:plasma_re", "plasma_re")))

#write.csv(data, "/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/data.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#data <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/data.txt", header =T, dec='.', sep=',')


###########3 PLOT ##############

library(ggpubr)
t = ggplot(dt,aes(x = tipo, y = len, fill = names))  +
  geom_col() +
  scale_fill_manual(values = ifelse(grepl("tissue", levels(dt$names), fixed = TRUE), "grey50", "white"))  + coord_flip() +
  theme(plot.title = element_blank(), legend.position = "none",
        axis.title=element_blank(),
        axis.text.x = element_blank())

pb =  ggplot(dpb,aes(x = tipo, y = len, fill = names))  +
  geom_col() +
  scale_fill_manual(values = ifelse(grepl("plasma_bl", levels(dpb$names), fixed = TRUE), "grey50", "white")) +coord_flip() +
  theme(plot.title = element_blank(), legend.position = "none",
        axis.title=element_blank(),
        axis.text.x = element_blank())

po = ggplot(dpo,aes(x = tipo, y = len, fill = names))  +
  geom_col() +
  scale_fill_manual(values = ifelse(grepl("plasma_po", levels(dpo$names), fixed = TRUE), "grey50", "white"))+coord_flip() +
  theme(plot.title = element_blank(), legend.position = "none",
        axis.title=element_blank(),
        axis.text.x = element_blank())
p =  ggplot(dp,aes(x = tipo, y = len, fill = names))  +
  geom_col() +
  scale_fill_manual(values = ifelse(grepl("plasma_re", levels(dp$names), fixed = TRUE), "grey50", "white"))+coord_flip() +
  theme(plot.title = element_blank(), legend.position = "none",
        axis.title=element_blank(),
        axis.text.x = element_blank())



tiff(file="/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots//enriquecimiento_INCLIVA_vars.tiff", width = 18, height = 8,units = 'cm', res = 400)
plot_grid(t, pb, po, p,  align = "v", 
          axis = "lr",
          nrow = 4)

dev.off()

###########################33

write.csv(data, "/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/counts_INCLIVA.txt",col.names=F, row.names = T, quote = FALSE)

data2 <- read.csv2("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/counts_INCLIVA.txt", header =T, dec='.', sep=',')
data2$names = factor(data2$names, levels = rev(c("tissue", "plasma_bl","tissue:plasma_bl",  "tissue:plasma_po" , "tissue:plasma_bl:plasma_po",
                                                 "plasma_bl:plasma_po", "tissue:plasma_re",  "plasma_bl:plasma_re", 
                                                 "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re",
                                                 "tissue:plasma_po:plasma_re",
                                                 "plasma_bl:plasma_po:plasma_re", 
                                                 "plasma_po", 
                                                 "plasma_po:plasma_re", 
                                                 "plasma_re")))

data2$tipo = factor(data2$tipo, levels = c("Tissue", "Plasma-BL", "Plasma-PO","Plasma-RE"))

col <- c("#DC0000FF", "#E64B35FF",  "#DC0000B2",   "#3C5488B2", "#8491B4B2",   "#3C5488D8", "#3C5488FF", "#4DBBD5FF",  
         "#4DBBD5B2",  "#00A087FF",  "#91D1C2FF",  "#91D1C2B2",  "#E64B3533", "#B09C85B2", "#B09C85FF")


svg(file="/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/fish_plots/FISH_plot_INCLIVA_1.svg", width = 10, height = 8)
ggplot(data2, aes(fill=names, y=len, x=tipo)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = col) + ylim(0, 8000) + scale_y_continuous(limits = c(0, 8000), breaks = c(0, 2500, 5000, 7500)) 

dev.off()







######################3 INDIVIDUALES MOMA#####################


muestras = c("M1","M2","M3","M5","M6","M7","M8","M9","M10","M11","M12","M13","M14","M15")
lista = list()
myplots <- vector('list', length(muestras))
i = 1
for ( mu in muestras) {
  varst <- read.csv2(paste0("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/ind/", mu, "_tissue_genes.txt"), header =F, dec='.', sep='\t')
  varsvt = as.character(varst$V1)
  
  
  varspb <- read.csv2(paste0("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/ind/", mu, "_plasma-bl_genes.txt"), header =F, dec='.', sep='\t')
  varsvpb = as.character(varspb$V1)
  
  varspo <- read.csv2(paste0("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/ind/", mu, "_plasma-po_genes.txt"), header =F, dec='.', sep='\t')
  varsvpo = as.character(varspo$V1)
  
  varsp <- read.csv2(paste0("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/ind/", mu, "_plasma_genes.txt"), header =F, dec='.', sep='\t')
  varsvp = as.character(varsp$V1)
  
  list_of_variants = list(tissue = varsvt, plasma_bl = varsvpb, plasma_po = varsvpo, plasma_re = varsvp)
  
  li = overlapGroups(list_of_variants)
  
  ta = attr(li, "elements")
  
  li[1]
  tissue = attr(li, "elements")[li[[1]]]
  
  
  
  d = as.data.frame(lengths(li))
  library(tibble)
  d <- tibble::rownames_to_column(d, "names")
  dt= cbind(d, rep("Tissue", length(rownames(d))))
  colnames(dt) = c("names", "len", "tipo")
  dpb= cbind(d, rep("Plasma-BL", length(rownames(d))))
  colnames(dpb) = c("names", "len", "tipo")
  dpo= cbind(d, rep("Plasma-PO", length(rownames(d))))
  colnames(dpo) = c("names", "len", "tipo")
  dp= cbind(d, rep("Plasma-RE", length(rownames(d))))
  colnames(dp) = c("names", "len", "tipo")
  
  data = rbind(dt, dpb, dpo, dp)
  
  
  
  
  
  
  dt$names = factor(dt$names, levels = rev(c("tissue", "tissue:plasma_bl", "tissue:plasma_bl:plasma_po",
                                             "tissue:plasma_po" , "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re", "tissue:plasma_po:plasma_re",
                                             "tissue:plasma_re","plasma_bl","plasma_bl:plasma_po", "plasma_bl:plasma_po:plasma_re",
                                             "plasma_bl:plasma_re", "plasma_po", "plasma_po:plasma_re", "plasma_re")))
  
  dpb$names = factor(dpb$names, levels = rev(c("tissue", "tissue:plasma_bl", "tissue:plasma_bl:plasma_po",
                                               "tissue:plasma_po" , "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re", "tissue:plasma_po:plasma_re",
                                               "tissue:plasma_re","plasma_bl","plasma_bl:plasma_po", "plasma_bl:plasma_po:plasma_re",
                                               "plasma_bl:plasma_re", "plasma_po", "plasma_po:plasma_re", "plasma_re")))
  
  dpo$names = factor(dpo$names, levels = rev(c("tissue", "tissue:plasma_bl", "tissue:plasma_bl:plasma_po",
                                               "tissue:plasma_po" , "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re", "tissue:plasma_po:plasma_re",
                                               "tissue:plasma_re","plasma_bl","plasma_bl:plasma_po", "plasma_bl:plasma_po:plasma_re",
                                               "plasma_bl:plasma_re", "plasma_po", "plasma_po:plasma_re", "plasma_re")))
  
  dp$names = factor(dp$names, levels = rev(c("tissue", "tissue:plasma_bl", "tissue:plasma_bl:plasma_po",
                                             "tissue:plasma_po" , "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re", "tissue:plasma_po:plasma_re",
                                             "tissue:plasma_re","plasma_bl","plasma_bl:plasma_po", "plasma_bl:plasma_po:plasma_re",
                                             "plasma_bl:plasma_re", "plasma_po", "plasma_po:plasma_re", "plasma_re")))
  
  #write.csv(data, "/home/jmartin/Documentos/exomas/MOMA/INCLIVA/metastasis/upset_plots/data.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  #data <- read.csv2("/home/jmartin/Documentos/exomas/MOMA/INCLIVA/metastasis/upset_plots/data.txt", header =T, dec='.', sep=',')
  
  
  write.csv(data, paste0("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/ind/counts_", mu, ".txt"),col.names=F, row.names = T, quote = FALSE)
  
  system(paste0("python /home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/ind/reorder_counts.py /home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/ind/counts_",
                mu, ".txt /home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/ind/counts2_", mu, ".txt"), wait=TRUE)
  
  data2 <- read.csv2(paste0("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/ind/counts2_", mu, ".txt"), header =T, dec='.', sep=',')
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
  svg(file=paste0("/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/fish_plots/ind/", mu, "_fish_plot.svg"), width = 10, height = 8)
  g = ggplot(data2, aes(fill=names, y=len, x=tipo)) + 
    geom_bar(position="stack", stat="identity") + scale_fill_manual(values = col) + xlab(mu) + ylim(0,2500)

  print(g)
  dev.off()
  
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

svg(file="/home/jmartin/Documents/articulo/scripts/figuras/fish/Supp_Fig7A.svg", width = 20, height = 12)
plot_grid(a , legend, align = "h", 
          axis = "lr",
          ncol = 2,  rel_widths = c(6/8, 2/8))

dev.off()


svg(file="/home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/fish_plots/FISH_plot_MOMA_ind.svg", width = 20, height = 12)
plot_grid(a , legend, align = "h", 
          axis = "lr",
          ncol = 2,  rel_widths = c(6/8, 2/8))

dev.off()

######################3 INDIVIDUALES INCLIVA #####################


muestras = c("85","107","161","185","207","243","261")
lista = list()
myplots <- vector('list', length(muestras))
i = 1
for ( mu in muestras) {
  varst <- read.csv2(paste0("/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/ind/", mu, "_tissue_genes.txt"), header =F, dec='.', sep='\t')
  varsvt = as.character(varst$V1)
  
  
  varspb <- read.csv2(paste0("/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/ind/", mu, "_plasma-bl_genes.txt"), header =F, dec='.', sep='\t')
  varsvpb = as.character(varspb$V1)
  
  varspo <- read.csv2(paste0("/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/ind/", mu, "_plasma-po_genes.txt"), header =F, dec='.', sep='\t')
  varsvpo = as.character(varspo$V1)
  
  varsp <- read.csv2(paste0("/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/ind/", mu, "_plasma_genes.txt"), header =F, dec='.', sep='\t')
  varsvp = as.character(varsp$V1)
  
  list_of_variants = list(tissue = varsvt, plasma_bl = varsvpb, plasma_po = varsvpo, plasma_re = varsvp)
  
  li = overlapGroups(list_of_variants)
  
  ta = attr(li, "elements")
  
  li[1]
  tissue = attr(li, "elements")[li[[1]]]
  
  
  
  d = as.data.frame(lengths(li))
  library(tibble)
  d <- tibble::rownames_to_column(d, "names")
  dt= cbind(d, rep("Tissue", length(rownames(d))))
  colnames(dt) = c("names", "len", "tipo")
  dpb= cbind(d, rep("Plasma-BL", length(rownames(d))))
  colnames(dpb) = c("names", "len", "tipo")
  dpo= cbind(d, rep("Plasma-PO", length(rownames(d))))
  colnames(dpo) = c("names", "len", "tipo")
  dp= cbind(d, rep("Plasma-RE", length(rownames(d))))
  colnames(dp) = c("names", "len", "tipo")
  
  data = rbind(dt, dpb, dpo, dp)
  
  
  
  
  
  
  dt$names = factor(dt$names, levels = rev(c("tissue", "tissue:plasma_bl", "tissue:plasma_bl:plasma_po",
                                             "tissue:plasma_po" , "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re", "tissue:plasma_po:plasma_re",
                                             "tissue:plasma_re","plasma_bl","plasma_bl:plasma_po", "plasma_bl:plasma_po:plasma_re",
                                             "plasma_bl:plasma_re", "plasma_po", "plasma_po:plasma_re", "plasma_re")))
  
  dpb$names = factor(dpb$names, levels = rev(c("tissue", "tissue:plasma_bl", "tissue:plasma_bl:plasma_po",
                                               "tissue:plasma_po" , "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re", "tissue:plasma_po:plasma_re",
                                               "tissue:plasma_re","plasma_bl","plasma_bl:plasma_po", "plasma_bl:plasma_po:plasma_re",
                                               "plasma_bl:plasma_re", "plasma_po", "plasma_po:plasma_re", "plasma_re")))
  
  dpo$names = factor(dpo$names, levels = rev(c("tissue", "tissue:plasma_bl", "tissue:plasma_bl:plasma_po",
                                               "tissue:plasma_po" , "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re", "tissue:plasma_po:plasma_re",
                                               "tissue:plasma_re","plasma_bl","plasma_bl:plasma_po", "plasma_bl:plasma_po:plasma_re",
                                               "plasma_bl:plasma_re", "plasma_po", "plasma_po:plasma_re", "plasma_re")))
  
  dp$names = factor(dp$names, levels = rev(c("tissue", "tissue:plasma_bl", "tissue:plasma_bl:plasma_po",
                                             "tissue:plasma_po" , "tissue:plasma_bl:plasma_po:plasma_re","tissue:plasma_bl:plasma_re", "tissue:plasma_po:plasma_re",
                                             "tissue:plasma_re","plasma_bl","plasma_bl:plasma_po", "plasma_bl:plasma_po:plasma_re",
                                             "plasma_bl:plasma_re", "plasma_po", "plasma_po:plasma_re", "plasma_re")))
  
  #write.csv(data, "/home/jmartin/Documentos/exomas/MOMA/INCLIVA/metastasis/upset_plots/data.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  #data <- read.csv2("/home/jmartin/Documentos/exomas/MOMA/INCLIVA/metastasis/upset_plots/data.txt", header =T, dec='.', sep=',')
  
  
  write.csv(data, paste0("/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/ind/counts_", mu, ".txt"),col.names=F, row.names = T, quote = FALSE)
  
  system(paste0("python /home/jmartin/Documents/EXOMAS/MOMA/resultados/oncogenic/upset_plots/ind/reorder_counts.py /home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/ind/counts_",
                mu, ".txt /home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/ind/counts2_", mu, ".txt"), wait=TRUE)
  
  data2 <- read.csv2(paste0("/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/ind/counts2_", mu, ".txt"), header =T, dec='.', sep=',')
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
  svg(file=paste0("/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/ind/", mu, "_fish_plot.svg"), width = 10, height = 8)
  g = ggplot(data2, aes(fill=names, y=len, x=tipo)) + 
    geom_bar(position="stack", stat="identity") + scale_fill_manual(values = col) + xlab(mu) + ylim(0, 1500)
  
  print(g)
  dev.off()
  
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


svg(file="/home/jmartin/Documents/articulo/scripts/figuras/fish/Supp_Fig7B.svg", width = 20, height = 12)
plot_grid(a , legend, align = "h", 
          axis = "lr",
          ncol = 2,  rel_widths = c(6/8, 2/8))
dev.off()


svg(file="/home/jmartin/Documents/EXOMAS/MOMA/INCLIVA/MSI_REMOVING/upset_plots/FISH_plot_INCLIVA_ind.svg", width = 20, height = 12)
plot_grid(a , legend, align = "h", 
          axis = "lr",
          ncol = 2,  rel_widths = c(6/8, 2/8))
dev.off()



