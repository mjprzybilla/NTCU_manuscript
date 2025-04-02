################################################################################################################################################
##                                                                                                                      
##  MOUSE NTCU MODEL SINGLE-CELL GEX DATA VISUALIZATION - UMAP, Dotplot, and Cell proportions
##                                                                                                                      
##  Date: 03 Jan 2024                                                                                                                   
##  
##  Author: Ahmed Alhendi                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

# Load all required packages
library(Seurat)
library(data.table)
library(SeuratData)
library(SeuratWrappers)
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

# Set Seurat options
options(Seurat.object.assay.version = "v3")


############################################################################
##                     Load Functions
############################################################################

# Load custom plotting functions
source("scripts/scRNAseq_functions.R")


############################################################################
##           Visualize UMAP of All Cells 
############################################################################

# Read in Seurat object
seurat_integrated_all <- readRDS("data/scRNAseq/NTCU_GEX_Trachea_processed_Celltype_seurat.rds")

# Define color palettes for all cells
colors.pheno <- c("#1f77b4", "#ff7f0e")
colors.celltype <- c('#1f77b4', '#ff7f0e', '#2ca02c')
colors.sample <- c('#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd', '#8c564b', '#e377c2')


# Custom theme for all plots
custom_theme <- theme_classic() +
  theme(
    text = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  )

# UMAP: Cells colored by Phenotype
Idents(seurat_integrated_all) <- "phenotype"
levels(seurat_integrated_all) <- c("Control", "NTCU-treated")

plot.pheno <- DimPlot(seurat_integrated_all, reduction = "X_umap", 
                      cols = colors.pheno, alpha = 0.75, pt.size = 0.1, shuffle = TRUE) +
  custom_theme +
  guides(color = guide_legend(override.aes = list(size = 6, alpha = 1)))

plot.pheno <- LabelClusters(plot = plot.pheno, id = "ident", repel = TRUE, size = 6)

ggsave(plot.pheno, filename = "Mouse_UMAP_All_cells_colored_by_phenotypes.pdf", width = 8, height = 7)

# UMAP: Cells colored by Cell Type
Idents(seurat_integrated_all) <- "Cell_types"

plot.celltype <- DimPlot(seurat_integrated_all, reduction = "X_umap", 
                         cols = colors.celltype, alpha = 0.75, pt.size = 0.1) +
  custom_theme+
  guides(color = guide_legend(override.aes = list(size = 6, alpha = 1)))

plot.celltype <- LabelClusters(plot = plot.celltype, id = "ident", repel = TRUE, size = 6)

ggsave(plot.celltype, filename = "Mouse_UMAP_All_cells_colored_by_cell_types.pdf", width = 8, height = 7)

# UMAP: Cells colored by Sample ID
Idents(seurat_integrated_all) <- "patientID"

plot.sampleID <- DimPlot(seurat_integrated_all, reduction = "X_umap", 
                         cols = colors.sample, alpha = 0.75, pt.size = 0.1) +
  custom_theme +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 6, alpha = 1)))

ggsave(plot.sampleID, filename = "Mouse_UMAP_All_cells_colored_by_SampleID.pdf", width = 8, height = 7)


############################################################################
##                    Visualize UMAP of Epithelial cells
############################################################################

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


# UMAP plot coloured by phenotype
colors.pheno <- c("#1f77b4", "#ff7f0e")
Idents(object = seurat_integrated) <- "phenotype"
levels(seurat_integrated) <- c("Control", "NTCU-treated")

plot1 <- DimPlot(seurat_integrated, reduction = "umap", cols = colors.pheno, 
alpha = 0.75, pt.size = 0.1, shuffle = TRUE) +
custom_theme +
  guides(color = guide_legend(override.aes = list(size = 6, alpha = 1)))

# Add cluster labels
plot.pheno <- LabelClusters(plot = plot1, id = "ident", repel = TRUE, size = 6)

# Save UMAP plot
ggsave(plot.pheno, file = "Mouse_UMAP_epithelial_colored_by_phenotypes.pdf", width = 8, height = 7)

# UMAP plot coloured by cluster cell types
Idents(seurat_integrated) <- "celltypes"
levels(seurat_integrated) <- df_color_clusters$labels

plot.cluster.umap <- DimPlot(seurat_integrated, reduction = "umap", 
label = TRUE, label.size = 5) +
  custom_theme +
  scale_color_manual(values = df_color_clusters$colors, labels = df_color_clusters$labels) + 
  guides(color = guide_legend(override.aes = list(size = 8)))

# Save UMAP plot by cell type clusters
ggsave(plot.cluster.umap, file = "Mouse_UMAP_colored_by_cluster_cell_types.pdf", width = 10, height = 7)


############################################################################
##           Visualize Marker Expression as Dotplots
############################################################################

# Visualize cell markers of all cells
all_markers <- fread("data/scRNAseq/Mouse_all_cell_selected_markers.tsv")
features <- all_markers$GeneID
Idents(seurat_integrated_all) <- "Cell_types"
plot.markers.all <- plotExp(object=seurat_integrated_all, features=features)
ggsave(plot.markers, file = "Mouse_all_cells_markers_dotplot.pdf", width=10, height=4)


# Visualize cell markers of epithelial cells
df.markers <- fread("data/scRNAseq/Mouse_epithelial_selected_markers.csv")

# Set cell type identity
Idents(seurat_integrated) <- "celltypes"
levels(seurat_integrated) <- df_color_clusters$labels
plot.markers.epithelial <- plotExp(object = seurat_integrated, features = df.markers$GeneID)
ggsave(plot.markers.epithelial, file = "Mouse_epithelial_markers_dotplot.pdf", width = 16, height = 7)


############################################################################
##           Visualize Cell Number and Proportions
############################################################################

# Calculate proportions

meta <- seurat_integrated@meta.data

dat <- meta %>% dplyr::group_by(phenotype, celltypes) %>% dplyr::count()
total.counts <- meta %>% dplyr::group_by(phenotype) %>% dplyr::count()

dat <- left_join(dat, total.counts, by="phenotype") %>%
mutate(prop=n.x/n.y)

# Define cell type order

dat <- dat %>% 
    mutate(celltypes = factor(celltypes, levels = cells.type.order))


# Plot proportions

alpha = 0.3
dat = dat
g <- "phenotype"
stat.by <- "celltypes"

dat_area <- dat[rep(seq_len(nrow(dat)), each = 2), , drop = FALSE]
dat_area[[g]] <- as.numeric(as.factor(dat_area[[g]]))
dat_area[seq(1, nrow(dat_area), 2), g] <- dat_area[seq(1, nrow(dat_area), 2), g] - 0.3
dat_area[seq(2, nrow(dat_area), 2), g] <- dat_area[seq(2, nrow(dat_area), 2), g] + 0.3

dat[[g]] <- as.numeric(as.factor(dat[[g]]))


p1 <- ggplot(dat, aes(x = .data[[g]], y = prop * 100, fill = .data[[stat.by]])) +
    geom_area(
    data = dat_area, mapping = aes(x = .data[[g]], fill = .data[[stat.by]]),
    alpha = alpha, color = "grey50") +

geom_col(aes(fill = .data[[stat.by]]), width = 0.6, color = "black")+
theme_cowplot(font_size = 22 )+
theme(axis.text.x =element_text(size=22, color="black", angle = 60, vjust = 1, hjust=1),
                              legend.position="right",
                               legend.title=element_blank())+
scale_x_continuous(breaks = 1:2, labels = c("Control", "NTCU-treated"))+
ylab("Proportion of cells")+xlab("")+scale_fill_manual(values = colors)


p2 <- ggplot(dat, aes(x = .data[[g]], y = n.x, fill = .data[[stat.by]])) +
    geom_area(
    data = dat_area, mapping = aes(x = .data[[g]], fill = .data[[stat.by]]),
    alpha = alpha, color = "grey50") +

geom_col(aes(fill = .data[[stat.by]]), width = 0.6, color = "black")+
theme_cowplot(font_size = 22 )+
theme(axis.text.x =element_text(size=22, color="black", angle = 60, vjust = 1, hjust=1),
                              legend.position="right",
                               legend.title=element_blank())+
scale_x_continuous(breaks = 1:2, labels = c("Control", "NTCU-treated"))+
ylab("Number of cells")+xlab("")+scale_fill_manual(values = colors)

ggsave(p1, file="Mouse_epithelial_Cell_Cell_Proportion.pdf", width=6, height=7)
ggsave(p2, file="Mouse_epithelial_Cell_Number.pdf", width=6, height=7)

