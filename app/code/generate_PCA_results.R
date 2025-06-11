library(shiny)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(Hmisc)
library(janitor)
library(Matrix)

## This file is to generate results for PCA visualization
## That way you don't need to run PCA every time (it's slow)

set.seed(123)

cell_types <- c("Acinar", "Active Stellate", "Alpha", "Beta", "Cycling Alpha", "Delta", "Ductal", "Endothelial", "Gamma + Epsilon", "Immune (Macrophages)", "MUC5B+ Ductal", "Quiescent Stellate")

for (c in cell_types) {
  # Get pseudobulk data for this cell type
  fname <- paste0("../outputs/CPM_matrices/", c, "_pseudobulk_cpm.csv")
  pseudo_df <- read.csv(fname)

  # Get all gene names
  genes <- colnames(pseudo_df)[76:ncol(pseudo_df)]
  pseudo_df[,genes] <- lapply(pseudo_df[,genes], as.numeric) # Make sure columns are numeric
  
  # Get gene expr. matrix to do PCA
  to_plot <- pseudo_df[,genes]
  
  # Do PCA
  res_pca <- prcomp(to_plot, scale = TRUE)
  
  # Save RDS for factor contributions
  setwd("../outputs/PCA_results/")
  fname <- paste0(c, "_all_PCA_results.rds")
  saveRDS(res_pca, fname)
  
  # Merge PCA results with metadata
  res <- res_pca$x %>% as.data.frame()
  res$samples <- pseudo_df$samples
  res <- res %>% merge(pseudo_df[,1:75], on = "samples")
  
  # Save PC results merged with metadata
  fname <- paste0(c, "_PCA_results.csv")
  write.csv(res, fname)
}
