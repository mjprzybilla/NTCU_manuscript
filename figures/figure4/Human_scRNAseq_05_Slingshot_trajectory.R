################################################################################################################################################
##                                                                                                                      
##  Human SINGLE-CELL GEX data - SLINGSHOT CELL LINEAGE INFERENCE FOR EPITHELIAL CELLS 
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
library(slingshot)
library(dplyr)
library(ggplot2)
library(cowplot)
library(condiments)
library(reshape2)
library(tidyr)
library(scales)
library(ggpubr)

set.seed(123)

# Load the Seurat object
seurat_integrated <- readRDS("Human_tracheal_epithelial_harmony_seurat.rds")


# Subset the data based on cell types of interest
cluster.to.exclude <- c('Ciliated', 'Deuterosomal', 'Ionocyte', 'Neuroendocrine')

# Subset Seurat object for these specific cell types
seurat.obj <- subset(seurat_integrated, Cluster_Annotations %in% cluster.to.exclude, invert = TRUE)

# Filter cluster labels to match the selected cell types
df_color_clusters.filtered <- df_color_clusters %>% 
  filter(!labels %in% cluster.to.exclude)

# Convert Seurat object to SingleCellExperiment object
sce <- as.SingleCellExperiment(seurat.obj, assay = "RNA")

# Bind UMAP results and metadata into a dataframe
df <- bind_cols(
  as.data.frame(reducedDims(sce)$UMAP),
  as.data.frame(colData(sce))
) %>% sample_frac(1)

rownames(df) <- df$cell_id

# Order dataframe by UMAP coordinates
df <- df[order(match(df$cell_id, rownames(dimred))),]

# Plot Slingshot trajectories on UMAP
ggplot(df, aes(x = umap_1, y = umap_2, col = Cluster_Annotations)) +
  geom_point() +
  scale_color_manual(values = df_color_clusters.filtered$colors) +
  geom_path(data = slingCurves(sds, as.df = TRUE) %>% arrange(Order) %>% filter(Lineage == 1),
            aes(group = Lineage), col = "black", size = 1, linetype = 1) +
  geom_path(data = slingCurves(sds, as.df = TRUE) %>% arrange(Order) %>% filter(Lineage == 2),
            aes(group = Lineage), col = "black", size = 1, linetype = 1) +
  geom_path(data = slingCurves(sds, as.df = TRUE) %>% arrange(Order) %>% filter(Lineage == 3),
            aes(group = Lineage), col = "black", size = 1, linetype = 1) +
  cowplot::theme_cowplot(font_size = 22)


# Dimensionality reduction (UMAP) for cells
dimred <- seurat.obj@reductions$umap@cell.embeddings
clustering <- seurat.obj@meta.data$Cluster_Annotations

# Perform Slingshot lineage inference starting from 'Basal cycling' cluster
sds <- slingshot(dimred, clustering, start.clus = 'Basal cycling')

# Perform topology test to compare lineage progression between smoking conditions
top_res <- topologyTest(sds = sds, conditions = df$Smoking_status)
print(paste0("Topology Test stats: ", top_res))

# Perform progression test to check if smoking status affects lineage progression
prog_res.prolif <- progressionTest(sds, conditions = df$Smoking_status, global = TRUE, lineages = TRUE)
prog_res.prolif$FDR <- p.adjust(prog_res.prolif$p.value, method = "BH")
print(paste0("Progression Test stats: ", prog_res.prolif))

# Perform cell fate selection test between smoking conditions
fate_sel_res.prof <- fateSelectionTest(sds, conditions = df$Smoking_status, pairwise = TRUE)
fate_sel_res.prof$FDR <- p.adjust(fate_sel_res.prof$p.value, method = "BH")
print(paste0("Cell fate stats: ", fate_sel_res.prof))

# Calculate and plot pseudotime distributions
psts.prolif <- slingPseudotime(sds) %>%
  as.data.frame() %>%
  mutate(cells = rownames(.),
         conditions = df$Smoking_status) %>%
  pivot_longer(starts_with("Lineage"), values_to = "pseudotime", names_to = "lineages")

psts.plotp <- ggplot(psts.prolif, aes(x = pseudotime, fill = conditions)) +
  geom_density(alpha = .5) +
  scale_fill_manual(values = c("blue", "darkorange")) +
  facet_wrap(~lineages, scales = "free") +
  cowplot::theme_cowplot(font_size = 22) +
  theme(legend.position = "top")

# Save pseudotime distribution plot
ggsave(psts.plotp, file = "Human_Pseudotime_distribution_Slingshot.pdf", width = 14, height = 5)

# Mean curve weight calculation and plotting
df_w <- condiments:::.sling_reassign(sds) %>% 
  as.data.frame() %>%
  mutate(cell_id = rownames(.))

df_w <- df_w %>% pivot_longer(starts_with("Lineage"), names_to = "Curve", values_to = "weights")
df_w <- df_w %>% left_join(df, by = "cell_id") %>% 
  select(cell_id, Curve, weights, Smoking_status)

# Normalising weights per cell
df_w <- df_w %>%
  group_by(cell_id) %>%
  mutate(weights = weights / sum(weights)) %>%
  ungroup() %>%
  group_by(Smoking_status, Curve) %>%
  summarise(weights = mean(weights), .groups = NULL)

# Plot mean curve weights by smoking status
weights.plot.prolif <- ggplot(df_w, aes(x = Curve, fill = Smoking_status, y = weights)) +
  geom_col(position = "dodge", color = "black") +
  scale_fill_manual(values = c("blue", "darkorange")) +
  labs(x = "", y = "Mean curve weight") +
  cowplot::theme_cowplot(font_size = 22) +
  theme(legend.position = "top")

# Save mean curve weight plot
ggsave(weights.plot.prolif, file = "Human_Mean_weight_Slingshot.pdf", width = 8, height = 6)

# Boxplot of pseudotime values by cluster annotations
lineages <- getLineages(dimred, clustering, start.clus = 'Basal cycling')
curves <- getCurves(lineages)

sp.data <- as.data.frame(slingPseudotime(curves))
sp.data$cell_id <- rownames(sp.data)

# Add cluster annotations and condition information
sp.data <- left_join(sp.data, Cluster_Annotations.ids, by = "cell_id")
sp.data$condition <- ifelse(sp.data$cell_id %in% group.ids$cell_id[group.ids$Smoking_status == "Non-smoker"], 
                            "Non-smoker", "Current-smoker")

# Reshape and summarise data for boxplot
sp.data <- reshape2::melt(sp.data)

sp.data.sum <- sp.data %>%
  filter(!is.na(value)) %>%
  group_by(Cluster_Annotations) %>%
  summarise(median_p = median(value)) %>%
  arrange(median_p)

# Reorder factors based on summary
sp.data <- sp.data %>% 
  mutate(Cluster_Annotations = factor(Cluster_Annotations, levels = sp.data.sum$Cluster_Annotations))

# Reorder colour labels for consistent plotting
df_color_clusters.filtered <- df_color_clusters.filtered[order(match(df_color_clusters.filtered$labels, 
                                                                     sp.data.sum$Cluster_Annotations)), ]

sp.data <- sp.data[sp.data$Cluster_Annotations %in% cells.type.order, ]

# Plot pseudotime as boxplot
plot.p <- ggplot(sp.data, aes(value, Cluster_Annotations, fill = Cluster_Annotations)) + 
  geom_boxplot() +
  scale_fill_manual(values = df_color_clusters.filtered$colors) +
  theme_cowplot(font_size = 22) +
  labs(x = "Pseudotime", y = "") + 
  theme(legend.position = "none")

# Save pseudotime boxplot
ggsave(plot.p, file = "Human_Peudotime_boxplot_Slingshot.pdf", width = 6, height = 6)


# Save the final trajectory object for future use
saveRDS(sds, file="slingshot_trajectory_Human_epithelail.rds")

