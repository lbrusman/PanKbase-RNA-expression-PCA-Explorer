suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(janitor))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(Matrix))


set.seed(123)
# Read in data
pankbase_data <- readRDS("/tscc/nfs/home/lebrusman/Gaulton_lab/data/ADA_object/250424_ADA_object_metadata_v3_3.rds")
colnames(pankbase_data@meta.data) <- colnames(pankbase_data@meta.data) %>% make_clean_names()

setwd("/tscc/nfs/home/lebrusman/Gaulton_lab/output/for_PCA_app/")

# Get celltype proportions
props <- prop.table(table(pankbase_data$cell_type, pankbase_data$samples), margin=2) %>% as.data.frame()
colnames(props) <- c("cell_type", "samples", "percent")

props_wide <- props %>% pivot_wider(id_cols = samples, names_from = cell_type, values_from = percent)
colnames(props_wide)[colnames(props_wide) != "samples"] <- paste0(colnames(props_wide)[colnames(props_wide) != "samples"], "_proportion")

# Get celltype counts in case you want to filter this way
all_ncells <- table(pankbase_data$cell_type, pankbase_data$samples) %>% as.data.frame()
colnames(all_ncells) <- c("cell_type", "samples", "n_cells")

all_ncells_wide <- all_ncells %>% pivot_wider(id_cols = samples, names_from = cell_type, values_from = n_cells)
colnames(all_ncells_wide)[colnames(all_ncells_wide) != "samples"] <- paste0(colnames(all_ncells_wide)[colnames(all_ncells_wide) != "samples"], "_ncells")


cell_types <- unique(pankbase_data@meta.data$cell_type)
Idents(pankbase_data) <- "cell_type"
for (c in cell_types) {
    # Subset Seurat object to only cell type of interest
    seurat_obj <- subset(pankbase_data, idents = c)
  	# If you want to subset genes
    rawcts <- LayerData(seurat_obj, assay = "RNA", layer = "counts")
    countcells <- rowSums(rawcts != 0)
    countcells <- countcells[countcells >= ncol(seurat_obj)*0.05] # Get genes expressed in min 5% of cells in cluster
    totalcts <- rawcts[which(rowSums(rawcts)>=100),] # Subset to min of 100 total counts for that gene
    keep_genes <- intersect(names(countcells), rownames(totalcts)) #keep only genes that meet these criteria

    # Subset seurat obj to those genes
    seurat_obj <- subset(seurat_obj, features = keep_genes)

    #get seurat obj metadata
    pankbase_cols <- colnames(seurat_obj@meta.data)[grep("pan_kbase", colnames(seurat_obj@meta.data))]
    metadata <- seurat_obj@meta.data[,c("samples", "treatments", "samples_split", "chemistry", "source", pankbase_cols)] %>% distinct(samples, treatments, samples_split, chemistry, source, .keep_all = TRUE) %>% as.data.frame()
    print(head(metadata))

    
    metadata <- metadata %>% merge(props_wide, on = "samples")
    metadata <- metadata %>% merge(all_ncells_wide, on = c("samples"))
    metadata <- metadata %>% filter(samples %in% unique(seurat_obj$samples))
    
    #for cpm
    pseudo_df <- AggregateExpression(seurat_obj, assay = "RNA", group.by = "samples", normalization.method = "RC", scale.factor = 1e6, return.seurat = TRUE)[["RNA"]]$data %>% as.data.frame()

    # Merge with metadata
    pseudo_df <- pseudo_df %>% t() %>% as.data.frame()
    pseudo_df$samples <- rownames(pseudo_df)

    all_df <- metadata %>% merge(pseudo_df, on = "samples")

    fname <- paste0(c, "_pseudobulk_cpm.csv")
    write.csv(all_df, fname)
    
}
