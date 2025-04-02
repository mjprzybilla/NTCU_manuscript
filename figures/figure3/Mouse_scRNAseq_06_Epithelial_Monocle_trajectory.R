################################################################################################################################################
##                                                                                                                      
##  MOUSE NTCU MODEL SINGLE-CELL GEX DATA - MONOCLE TRAJECTORY OF EPITHELIAL CELLS 
##                                                                                                                      
##  Date: 03 Feb 2024                                                                                                                   
##  
##  Author: Ahmed Alhendi                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

set.seed(123)


# Load all required packages
library(Seurat)
library(data.table)
library(monocle3)
library(monocle)
library(SeuratData)
library(SeuratWrappers)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(magrittr)
library(dplyr)
library(scales)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(viridis)
library(pheatmap)
library(reshape2)
library(RColorBrewer)

# Read in Seurat objects
seurat_integrated <- readRDS("data/scRNAseq/NTCU_Trachea_Epithelial_integrated_seurat.rds")


# Define cluster numbers and corresponding labels
clusters <- c(0:12)
labels <- c("Basal regenerative", "Basal proliferative", "Krt4/13+ cells", "Secretory", "Secretory Mecom+",
            "Basal Tgm2+", "Basal", "Basal Mecom+", "Deuterosomal cells", "Ciliated Cells",
            "lonocytes", "Tuft", "Neuroendocrine")

# Create dataframe for clusters and labels
df_clusters <- data.frame(clusters = clusters, labels = labels)

# Define cell type order for visualization
cells.type.order <- c('Basal proliferative', 'Basal', 'Basal Tgm2+', 'Basal Mecom+', 'Basal regenerative',
                      'Krt4/13+ cells', 'Secretory', 'Secretory Mecom+', 'Deuterosomal cells', 
                      'Ciliated Cells', 'Neuroendocrine', 'lonocytes', 'Tuft')

# Define custom colours for each cluster
colors <- c('#3333ff', '#0bb1da', '#ffcc00', '#ff527d', '#cc0099', '#32b34e', '#7a00cc', 
            '#ccaaff', '#99ffff', '#008080', '#ff66cc', '#ffff66', '#999999')

# Limit colours to match the number of clusters
colors <- colors[1:13]

# Create dataframe linking labels with their respective colours
df_colors <- data.frame(labels = cells.type.order, colors = colors)

# Merge labels and colours with clusters dataframe
df_color_clusters <- left_join(df_clusters, df_colors, by = "labels")
df_color_clusters


# Subset the data based on cell types of interest
cluster.included <- c('Basal proliferative', 'Basal', 'Basal Tgm2+', 'Basal Mecom+', 'Basal regenerative',
                      'Krt4/13+ cells', 'Secretory', 'Secretory Mecom+')

# Subset Seurat object for these specific cell types
seurat.obj <- subset(seurat_integrated, celltypes %in% cluster.included)

seurat.obj@meta.data$cell_id <- rownames(seurat.obj@meta.data)
cell_id <- seurat.obj@meta.data$cell_id

# Calculate the number of cells to keep (10% of the total cells)
num_cells <- length(cell_id)
num_cells_to_keep <- round(0.10 * num_cells)

# Randomly select the indices of cells to keep
indices_to_keep <- sample(1:num_cells, num_cells_to_keep)
cell_id_to_keep <- cell_id[indices_to_keep]


# Subset seurat 
seurat.obj <- subset(seurat.obj, cell_id %in% cell_id_to_keep)

# Filter cluster labels to match the selected cell types
df_color_clusters.filtered <- df_color_clusters %>% 
                            filter(labels %in% cluster.included)


# Convert Seurat object to a CellDataSet object for Monocle
cds <- as.CellDataSet(seurat.obj)

# Estimate size factors and dispersions for the dataset
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# Assign celltypes as a factor for trajectory analysis
pData(cds)$celltypes <- factor(pData(cds)$celltypes, levels = cells.type.order)

# Perform differential gene test to identify genes for clustering
clustering_DEG_genes <- differentialGeneTest(cds, fullModelFormulaStr = '~celltypes', cores = 8)

# Use top 400 genes from the differential expression results
clustering_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:400]

# Set the ordering genes for trajectory analysis
cds <- monocle::setOrderingFilter(cds, ordering_genes = clustering_ordering_genes)

# Reduce dimensionality using DDRTree method
cds <- monocle::reduceDimension(cds, method = 'DDRTree')

# Order the cells based on their trajectory
cds <- monocle::orderCells(cds)

# Plot the DDRTree, coloured by cell state
plot.cell.state <- monocle::plot_cell_trajectory(cds, color_by = "State") +
  theme_cowplot(font_size = 25) +
  theme(legend.position = "top")

# Save the plot as a PDF
ggsave(plot.cell.state, file = "Mouse_epithelial_cell_state_monocle.pdf", width = 6, height = 7)

# Save the trajectory object (cds) for later use
saveRDS(cds, "Mouse_epithelial_moncole2_ddrtree_cds.rds")

# Load the saved cds object to continue from where we left off
cds <- readRDS("Mouse_epithelial_moncole2_ddrtree_cds.rds")

# Function to assign the root state for the trajectory
GM_state <- function(cds) {
  if (length(unique(cds$State)) > 1) {
    T0_counts <- table(cds$State, cds$celltypes)[,"Basal proliferative"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return(1)
  }
}

# Order the cells based on root state
cds <- monocle::orderCells(cds, root_state = GM_state(cds))

# Plot DDRTree colored by cluster cell types
plot.cell.types <- monocle::plot_cell_trajectory(cds, color_by = "celltypes") +
  theme_cowplot(font_size = 25) +
  scale_color_manual(values = df_color_clusters.filtered$colors) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(legend.position = "none")

# Save the plot as a PDF
ggsave(plot.cell.types, file = "Mouse_epithelial_cell_types_monocle.pdf", width = 9, height = 7)

# Plot DDRTree colored by pseudotime
plot.pseudotime <- monocle::plot_cell_trajectory(cds, color_by = "Pseudotime") +
  theme_cowplot(font_size = 25) +
  theme(legend.position = "top") +
  scale_color_viridis(option = "plasma")

# Save the plot as a PDF
ggsave(plot.pseudotime, file = "Mouse_epithelial_pseudotime_monocle.pdf", width = 6, height = 7)

# Reorder celltypes in the cds object according to predefined celltype order
pData(cds)$celltypes <- factor(pData(cds)$celltypes, levels = cells.type.order)

# Plot the trajectory tree of the cells, grouped by phenotype
plot.tree <- plot_complex_cell_trajectory(cds[, ], color_by = 'celltypes', 
                                          show_branch_points = TRUE, 
                                          cell_size = 0.3, cell_link_size = 1) + 
  facet_wrap(~phenotype, nrow = 1) + 
  theme_cowplot(font_size = 25) +
  scale_color_manual(values = df_color_clusters.filtered$colors) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(legend.title = element_blank()) +
  labs(x = "Component 1", y = "Component 2")

# Save the plot as a PDF
ggsave(plot.tree, file = "Mouse_epithelial_trajectory_tree_monocle.pdf", width = 10, height = 6)

############################################################################
##           COMPARATIVE TRAJECTORY ANALYSIS BETWEEN CONTROL AND NTCU
############################################################################

# Split the cds into control and ntcu
control_cells <- row.names(pData(cds)[pData(cds)$phenotype == "Control",])
cds.control <- cds[,control_cells]
cds.control

ntcu_cells <- row.names(pData(cds)[pData(cds)$phenotype == "NTCU-treated",])
cds.ntcu <- cds[,ntcu_cells]
cds.ntcu


# NTCU treated cell abundance
state_cluster_stat <- table(pData(cds)[, c('State', 'celltypes')])

str(state_cluster_stat)
state_cluster_stat <- apply(state_cluster_stat, 2, function(x) x / sum(x))
state_cluster_stat_ordered <- t(state_cluster_stat)

pheatmap::pheatmap(state_cluster_stat_ordered, cluster_cols = F, cluster_rows = F, fontsize = 20,
                   color = colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Reds'))(10))

# CONTROL cell abundance
state_cluster_stat <- table(pData(cds.control)[, c('State', 'celltypes')])


state_cluster_stat.cnt <- apply(state_cluster_stat, 2, function(x) x )
      
state_cluster_stat.cnt <- t(state_cluster_stat.cnt)
state_cluster_stat.cnt


state_cluster_stat <- apply(state_cluster_stat, 2, function(x) x / sum(x))
state_cluster_stat_ordered <- t(state_cluster_stat)

state_cluster_stat_ordered.cnt <- state_cluster_stat_ordered                            
                            
pheatmap::pheatmap(state_cluster_stat_ordered, main="Control", cluster_cols = F, cluster_rows = F, fontsize = 20,
                   color = colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Blues'))(50))


pesudo_count=1e-5

log_fold_change <- log2((state_cluster_stat_ordered.ntcu + pesudo_count) / (state_cluster_stat_ordered.cnt + pesudo_count) )

log_fold_change_df <- as.data.frame(log_fold_change)

log_fold_change_melted <- reshape2::melt(log_fold_change_df, id.vars = "Cell_Type")

log_fold_change_melted$Cell_Type <- factor(log_fold_change_melted$Cell_Type, levels= rev(df_color_clusters$labels))
log_fold_change_melted


# Create heatmap of log2 changes in cell abundance over the trajactory
plot.cell.abundance <- ggplot(log_fold_change_melted, aes(x = variable, y = Cell_Type, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0, name = "Log2FC") +
  labs(title="Cell abundance",
       x = "Cell state",
       y = "") +
theme_cowplot(font_size = 22)

ggsave(plot.cell.abundance, file="Mouse_epithelial_trajectory_cell_abundance_monocle.pdf", width=8, height=6)                           

# Save the final trajectory object for future use
saveRDS(cds, file="final_monocle_trajectory_NTCU.rds")
