################################################################################################################################################
##                                                                                                                      
##  Human SINGLE-CELL GEX data - MILO DIFFERENTIAL CELL ABUNDANCE FOR EPITHELIAL CELLS
##                                                                                                                      
##  Date: 18 July 2024                                                                                                                   
##  
##  Author: Ahmed Alhendi                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################


# Load required libraries
library(Seurat)
library(SingleCellExperiment)
library(miloR)
library(dplyr)
library(ggplot2)
library(patchwork)
library(data.table)
library(ggbeeswarm)

# Load the Seurat object
seurat_integrated <- readRDS("Human_tracheal_epithelial_harmony_seurat.rds")

# Convert Seurat object to Assay class
seurat_integrated[["RNA"]] <- as(object = seurat_integrated[["RNA"]], Class = "Assay")

# Convert Seurat object to SingleCellExperiment object
sce <- as.SingleCellExperiment(seurat_integrated)

# Initialize Milo object for spatial analysis
epith_milo <- Milo(sce)

# Build the graph using k-nearest neighbours (k = 30) and a distance of 30
epith_milo <- buildGraph(epith_milo, k = 30, d = 30, reduced.dim = "HARMONY_RNA")

# Create neighbourhoods with a proportion of 0.1, refine the neighbourhoods
epith_milo <- makeNhoods(epith_milo, prop = 0.1, k = 30, d = 30, refined = TRUE, reduced_dims = "HARMONY_RNA")

# Plot the histogram of neighbourhood sizes
plotNhoodSizeHist(epith_milo)

# Count cells in the neighbourhoods, adding meta data for sample and origin
epith_milo <- countCells(epith_milo, meta.data = as.data.frame(colData(epith_milo)), sample = "orig.ident")

# Convert Smoking status and dataset to factors and define their levels
epith_milo@colData$Smoking_status <- factor(epith_milo@colData$Smoking_status, levels = c("Non-smoker", "Current-smoker"))
epith_milo@colData$dataset <- factor(epith_milo@colData$dataset, levels = unique(epith_milo@colData$dataset))

# Create the design matrix for analysis based on smoking status and dataset
epith_milo_design <- data.frame(colData(epith_milo))[, c("orig.ident", "Smoking_status", "dataset")]

epith_milo_design <- distinct(epith_milo_design)
rownames(epith_milo_design) <- epith_milo_design$orig.ident

# Print the design matrix
print(epith_milo_design)

# Calculate the distance between neighbourhoods
epith_milo <- calcNhoodDistance(epith_milo, d = 30, reduced.dim = "HARMONY_RNA")

# Perform differential analysis between neighbourhoods based on Smoking status
da_results <- testNhoods(epith_milo, design = ~ Smoking_status, design.df = epith_milo_design)

# Plot the differential analysis results as a scatter plot
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1)

# Build the neighbourhood graph based on differential analysis
epith_milo <- buildNhoodGraph(epith_milo)

# Plot single-cell UMAP coloured by Smoking status, with annotations
umap_pl <- plotReducedDim(epith_milo, dimred = "UMAPAFTERHARMONY_RNA", colour_by = "Smoking_status", 
                          text_by = "Cluster_Annotations", 
                          text_size = 3, point_size = 0.5) + 
  guides(fill = "none")

# Plot neighbourhood graph coloured by differential analysis results
nh_graph_pl <- plotNhoodGraphDA(epith_milo, da_results, layout = "UMAPAFTERHARMONY_RNA", alpha = 0.1) +
  scale_fill_gradientn(colors = c("darkblue", "white", "darkred"))

# Combine the UMAP and neighbourhood graph plots into one layout
combined.miloR.ump.plot <- umap_pl + nh_graph_pl + 
  plot_layout(guides = "collect")

# Save combined plot as a PDF
ggsave(combined.miloR.ump.plot, file = "MiloR.heatmap.nhood.plot.pdf", width = 14, height = 7)

# Annotate neighbourhoods with cluster annotations
da_results <- annotateNhoods(epith_milo, da_results, coldata_col = "Cluster_Annotations")

# Summarise the differential analysis results by cluster annotations
da_results_summary <- da_results %>%
  select(Cluster_Annotations, logFC, PValue, SpatialFDR) %>%
  group_by(Cluster_Annotations) %>%
  summarise(
    mean_logFC = median(logFC),
    mean_PValue = median(PValue),
    mean_SpatialFDR = median(SpatialFDR)
  ) %>%
  arrange(-mean_logFC, mean_SpatialFDR)

# Update Cluster_Annotations for mixed clusters with low fraction
da_results$Cluster_Annotations <- ifelse(da_results$Cluster_Annotations_fraction < 0.7, "Mixed", da_results$Cluster_Annotations)

# Filter out "Mixed" clusters and reorder annotations
da_results2 <- da_results %>% filter(Cluster_Annotations != "Mixed") %>%
  mutate(Cluster_Annotations = factor(Cluster_Annotations, 
                                      levels = rev(da_results_summary$Cluster_Annotations)))

# Plot the differential analysis results as a beeswarm plot
plot.final <- plotDAbeeswarm(da_results2, group.by = "Cluster_Annotations", alpha = 0.1) +
  scale_color_gradientn(colors = c("darkblue", "white", "darkred")) +
  xlab("") + 
  theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"))

# Save the beeswarm plot as a PDF
ggsave(plot.final, file = "MiloR.plotDAbeeswarm.human_trachea.pdf", width = 7, height = 7)
