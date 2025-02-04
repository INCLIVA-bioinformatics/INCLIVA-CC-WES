library(MSstats)
library(MSstatsConvert)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

names = as.list(c("107", "136", "259", "63", "185", "13", "242", "243", "189", "204", "49", "261", "104"))

all_abundance = as.data.frame(c("Suffix_Patient" = c(), "Patient" = c(), "metastasis" = c(), "primary" = c()))


for (name in names) {

  evidence <- read.csv2(paste0("/home/jmartin/Documents/articulo/scripts/repo_actualizado/INCLIVA-CC-WES/inputs/FIG5D/", name, "_evidence.txt"), header =T, dec='.', sep='\t')
  
  neo_data <- read.csv2(paste0("/home/jmartin/Documents/articulo/scripts/repo_actualizado/INCLIVA-CC-WES/inputs/FIG5D/", name, "_sequences.txt"), header =T, dec='.', sep='\t')
  
  protein_groups <- read.csv2(paste0("/home/jmartin/Documents/articulo/scripts/repo_actualizado/INCLIVA-CC-WES/inputs/FIG5D/", name, "_proteinGroups.txt"), header =T, dec='.', sep='\t')
  
  neoepitopes = unique(neo_data$Sequence)
  
  if ('Experiment' %in% names(evidence)) {
    evidence %>%
      distinct(`Raw.file`) %>%
      mutate(Condition = evidence %>%distinct(`Experiment`)) %>%
      arrange(Raw.file) %>%
      mutate(BioReplicate = 1:n()) %>%
      ungroup() %>%
      add_column(IsotypeLabelType="L")  -> annotation
    
  } else {
    
    evidence %>%
      distinct(`Raw.file`) %>%
      mutate(Condition = str_replace(Raw.file,"^.*-","")) %>%
      arrange(Raw.file) %>%
      group_by(Condition) %>%
      mutate(BioReplicate = 1:n()) %>%
      ungroup() %>%
      add_column(IsotypeLabelType="L")  -> annotation
    
    
    if (annotation$Raw.file[1] != annotation$Condition[1]) {
      annotation$Condition[[1]] = paste0("ET", name)
      annotation$Condition[[2]] = paste0("EM", name)
    }
  } 
  
  annotation
  
  MaxQtoMSstatsFormat(
    evidence = evidence,
    annotation = annotation,
    proteinGroups = protein_groups,
    removeFewMeasurements = FALSE,
    useUniquePeptide = FALSE,
  ) -> raw_data
  

  dataProcess(
    raw_data,
    logTrans = 2,
    n_top_feature = 1
  ) -> quantified_data
  

  # Extract the base peptide sequence (before the underscore) for matching
  quantified_data$FeatureLevelData$PEPTIDE_BASE <- sub("_.*", "", quantified_data$FeatureLevelData$PEPTIDE)
  
  # Filter the dataframe to include only the rows with target peptides
  quantified_neoepitopes <- quantified_data$FeatureLevelData[quantified_data$FeatureLevelData$PEPTIDE_BASE %in% neoepitopes, ]
  
  abund = quantified_neoepitopes %>%
    arrange(GROUP,originalRUN) %>%
    mutate(originalRUN=factor(originalRUN, levels=unique(originalRUN)))
  

  df <- replace(abund, is.na(abund), 0)
  data_abund = df[, c("GROUP", "ABUNDANCE")]
  

  meta <- data_abund$ABUNDANCE[data_abund$GROUP == paste0("EM", name)]
  primary <- data_abund$ABUNDANCE[data_abund$GROUP == paste0("ET", name)]
  

  result_test = t.test(meta, primary, paired = "TRUE",alternative = "less") # Prueba t
  result_test


  patient_abund <- data_abund %>%
    mutate(
      Condition = ifelse(grepl("ET", GROUP), "primary", "metastasis"),
      Patient = as.numeric(gsub("[^0-9]", "", GROUP))
    )
  
  
  patient_abund$ABUNDANCE = as.numeric(patient_abund$ABUNDANCE)
  
  patient_abund2 <- patient_abund %>%
    group_by(Patient, Condition) %>%
    mutate(suffix = row_number()) %>%
    ungroup() 
  
  patient_abund2 <- patient_abund2 %>%
    mutate(Suffix_Patient = paste0(Patient, "_", suffix))
  
  df_wide <- patient_abund2 %>%
    select(Suffix_Patient, Patient, Condition, ABUNDANCE) %>%
    pivot_wider(names_from = Condition, values_from = ABUNDANCE)
  
  patient_abundance <- df_wide %>%
    filter(!is.na(primary) | !is.na(metastasis))
  
  all_abundance = rbind(all_abundance, patient_abundance)
  

}


result_global = wilcox.test(all_abundance$metastasis, all_abundance$primary, paired = TRUE, alternative = "less")


data_long <- pivot_longer(all_abundance, cols = c("metastasis", "primary"), 
                          names_to = "group", 
                          values_to = "abundance")

t_test_results <- data_long %>%
  group_by(Patient) %>%
  summarise(p_value = t.test(abundance[group == "metastasis"], abundance[group == "primary"], paired = TRUE, alternative = "less")$p.value)


t_test_results <- t_test_results %>%
  mutate(significance = ifelse(p_value < 0.05, "Significativo", "No Significativo"))


data_mean <- data_long %>%
  group_by(Patient, group) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE)) %>%
  ungroup()


data_mean$group <- factor(data_mean$group, levels = c("primary", "metastasis"), 
                          labels = c("Primary", "Metastasis"))


data_mean <- left_join(data_mean, t_test_results, by = "Patient")

data_mean$Patient <- factor(data_mean$Patient, levels = sort(unique(data_mean$Patient), decreasing = TRUE))

svg(file="/home/jmartin/Documents/articulo/scripts/repo_actualizado/INCLIVA-CC-WES/figures/FIG5D.svg", width = 4, height = 8)

p = ggplot(data_mean, aes(x = group, y = factor(Patient), fill = mean_abundance)) +
  geom_tile(color = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue", guide = "colourbar", name = "Neoepitope\nAbundance") + 
  labs(x = "Time Point", 
       y = "Patient") +
  theme_minimal() + scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(size = 15), 
        legend.title=element_text(size=20),
        legend.text=element_text(size=0),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 20, colour = "black")) + 
  geom_text(aes(x = 2, y = factor(Patient),label = ifelse(significance == "Significativo", "*", "")), 
            color = "black", vjust = 1, hjust = -4.5, size = 4) + 
  annotate("segment", x = 1, xend = 2, y = -0.01, yend = -0.01, color = "black", linewidth = 0.5) + 
  annotate("segment", x = 1, xend = 1, y = -0.01, yend = 0.5, color = "black", linewidth = 0.3) + 
  annotate("segment", x = 2, xend = 2, y = -0.01, yend = 0.5, color = "black", linewidth = 0.3) +
  annotate("text", x = 1.5, y = 0.5, label = "p-value < 0.001", size = 4.7, vjust = 1.5, colour = "black") 
p
dev.off()




############ FIGURE 5E ##################


#libraries
library(pheatmap)

summary_df <- read.csv2("/home/jmartin/Documents/articulo/scripts/repo_actualizado/INCLIVA-CC-WES/inputs/input_FIG5E.csv", header =T, dec='.', sep=',')


#plot
heatmap_data <- summary_df %>%
  select(variable, mean_ratio_wt_prim, mean_ratio_mut_prim, p_value) %>%
  column_to_rownames("variable")

rownames(heatmap_data) <- paste(rownames(heatmap_data), "p=", round(heatmap_data$p_value, 4), sep = " ")

pheatmap(
  heatmap_data[, c("mean_ratio_wt_prim", "mean_ratio_mut_prim")],
  color = colorRampPalette(c("white", "blue"))(100),              
  show_rownames = TRUE,                                          
  fontsize = 10                                                   
)
