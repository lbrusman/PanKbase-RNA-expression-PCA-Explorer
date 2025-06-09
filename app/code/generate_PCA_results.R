library(shiny)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(Hmisc)
library(janitor)
library(Matrix)

## This file is to generate results for PCA visualization
## That way you don't need to run PCA every time (it's slow)

#Wrangle metadata
# Might not even need metadata
metadata <- read.csv("/Users/lbrusman/Desktop/Gaulton_lab/shiny_apps/PanKbase-RNA-expression-PCA/app/data/metadata_for_ADA.csv") 


cell_types <- c("Acinar", "Alpha", "Beta", "Cycling Alpha", "Delta", "Ductal", "Endothelial", "Gamma + Epsilon", "Immune (Macrophages)", "MUC5B+ Ductal", "Quiescent Stellate")

for (c in cell_types) {
  # Get pseudobulk data for this cell type and filter for only untreated samples
  fname <- paste0("/Users/lbrusman/Desktop/Gaulton_lab/shiny_apps/PanKbase-RNA-expression-PCA/app/data/CPM_matrices/", c, "_pseudobulk_cpm.csv")
  pseudo_df <- read.csv(fname)
  
  # Get all gene names
  genes <- colnames(pseudo_df)[12:ncol(pseudo_df)]
  pseudo_df[,genes] <- lapply(pseudo_df[,genes], as.numeric)
  
  #filter metadata to samples in this cell type
  metadata_c <- metadata
  
  #get new samples column
  #get meta_all ready to merge...
  metadata_c$samples_split <- ifelse(metadata_c$srr %in% pseudo_df$samples_split, metadata_c$srr, metadata_c$rrid)
  metadata_c <- metadata_c %>% filter(samples_split %in% pseudo_df$samples)
  metadata_c$samples <- metadata_c$samples_split
  metadata_c <- metadata_c %>% distinct(samples, .keep_all = TRUE)
  
  # Merge with all metadata
  pseudo_df <- metadata_c %>% merge(pseudo_df, on = "samples")
  
  # Get gene expr. matrix to do PCA
  to_plot <- pseudo_df[,genes]
  
  # Do PCA
  res_pca <- prcomp(to_plot, scale = TRUE)
  
  setwd("/Users/lbrusman/Desktop/Gaulton_lab/shiny_apps/PanKbase-RNA-expression-PCA/app/code/PCA_results/")
  fname <- paste0(c, "_all_PCA_results.rds")
  saveRDS(res_pca, fname)
  
  #merge PCA results with metadata
  res <- res_pca$x %>% as.data.frame()
  res$samples <- pseudo_df$samples
  res <- res %>% merge(metadata_c, on = "samples")
  
  #save results
  setwd("/Users/lbrusman/Desktop/Gaulton_lab/shiny_apps/PanKbase-RNA-expression-PCA/app/code/PCA_results/")
  fname <- paste0(c, "_PCA_results.csv")
  write.csv(res, fname)
}
