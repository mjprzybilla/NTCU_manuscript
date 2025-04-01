
################################################################################################################################################
##                                                                                                                      
##  Human TRACHEAL SINGLE-CELL GEX DATA VISULIZATION - FIGURE 4 & SUPP FIGURE 6
##                                                                                                                      
##  Date: 18 july 2024                                                                                                                   
##  
##  Author: Ahmed Alhendi                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

##############################################################################
##                         LOAD REQUIRED PACKAGES
##############################################################################

# List of required packages
required_packages <- c(
  # Core analysis
  "Seurat", "SeuratData", "SeuratWrappers", "data.table", "dplyr", 
  "tidyverse", "ggplot2", "reshape2", "scales", "ggpubr", "ggsignif",
  
  # Single-cell analysis
  "miloR", "SingleCellExperiment", "scater", 
  "slingshot", "tradeSeq", "condiments", "DESeq2",
  
  # Visualization
  "patchwork", "cowplot", "viridis", "pheatmap", "RColorBrewer", 
  "wesanderson", "plotrix", "grid", "gridExtra"
)

# Load all required packages
suppressMessages({
  lapply(required_packages, library, character.only = TRUE)
})

############################################################################
##                     Load Data & Functions
############################################################################


# Read in Seurat object
seurat_integrated <- readRDS("data/scRNAseq/Human_tracheal_epithelial_harmony_seurat.rds")

# Load custom plotting functions
source("scripts/scRNAseq_functions.R")


############################################################################
##                    Define Cluster Colors and Labels
############################################################################
                        
# Define cell type order for visualization
cells.type.order <- c('Basal cycling', 'Basal 1','Basal 2','Basal 3','Basal SMG','KRT13/KRT4',
                      'Suprabasal 1','Suprabasal 2','Suprabasal 3', 'Secretory','Serous',
                      'Deuterosomal','Ciliated','Ionocyte','Neuroendocrine'
                     )


# Define cluster numbers and corresponding labels
clusters <- c(0:14)

labels <-               c("Basal 1",
                         "Suprabasal 2",
                         "Basal 3", 
                         "Secretory",
                         "Suprabasal 1",
                         "Basal 2",
                         "KRT13/KRT4",
                         "Ciliated",
                         "Basal cycling",
                         "Serous",
                         "Suprabasal 3",
                         "Basal SMG",
                         "Deuterosomal",
                         "Ionocyte",
                         "Neuroendocrine"
                        )

# Create dataframe for clusters and labels
df_clusters <- data.frame(clusters = clusters, labels = labels)


# Define custom colours for each cluster
colors <- c('#1f77b4', '#aec7e8', '#AA336A', '#ff9896', '#ff7f0e', '#ffbb78', '#2ca02c', 
'#90EE90', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', 
'#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5')

# Limit colours to match the number of clusters
colors <- colors[1:15]

# Create dataframe linking labels with their respective colours
df_colors <- data.frame(labels = cells.type.order, colors = colors)

# Merge labels and colours with clusters dataframe
df_color_clusters <- left_join(df_clusters, df_colors, by = "labels")
df_color_clusters


############################################################################
##           Visualize UMAP of Epithelial Cells 
############################################################################


# UMAP plot coloured by dataset
Idents(seurat_integrated) <- "dataset"

plot.umap.dataset <- DimPlot(seurat_integrated, reduction = "umap", 
alpha = 0.75, pt.size = 0.1, shuffle = TRUE)+
theme(text=element_text(size=20, color="black"))

ggsave(plot.umap.dataset, file = "Human_UMAP_epithelial_colored_by_dataset.pdf", width = 8, height = 6)


# Plot UMAP with condition colourings
colors.condition <- c("#1f77b4", "#ff7f0e")
Idents(object = seurat_integrated) <- "Smoking_status"
levels(seurat_integrated) <- c("Non-smoker", "Current-smoker")

plot.umap.condition <- DimPlot(seurat_integrated, reduction = "umap", 
alpha = 0.75, pt.size = 0.1, shuffle = TRUE, cols=colors.condition)+
theme(text=element_text(size=20, color="black"))

ggsave(plot.umap.condition, file = "Human_UMAP_epithelial_colored_by_condition.pdf", width = 8, height = 6)


# Plot UMAP with SampleID (Donor)
Idents(object = seurat_integrated) <- "Donor"

plot.umap.condition <- DimPlot(seurat_integrated, reduction = "umap", 
alpha = 0.75, pt.size = 0.1, shuffle = TRUE)+
theme(text=element_text(size=20, color="black"))

ggsave(plot.umap.condition, file = "Human_UMAP_epithelial_colored_by_condition.pdf", width = 8, height = 6)


# UMAP plot coloured by cluster cell types
Idents(seurat_integrated) <- "Cluster_Annotations"
levels(seurat_integrated) <- df_color_clusters$labels

plot.cluster.umap <- DimPlot(seurat_integrated, reduction = "umap", cols=df_color_clusters$colors
alpha = 0.75, pt.size = 0.1, shuffle = TRUE)+
theme(text=element_text(size=20, color="black"))

# Save UMAP plot by cell type clusters
ggsave(plot.cluster.umap, file = "Human_UMAP_colored_by_cluster_cell_types.pdf", width = 10, height = 7)



############################################################################
##           Visualize Marker Expression as Dotplots
############################################################################

# Read in the marker genes
df.markers <- fread("data/scRNAseq/Human.selected_markers.csv")

# Set cell type identity
Idents(seurat_integrated) <- "Cluster_Annotations"
levels(seurat_integrated) <- df_color_clusters$labels

# Plot expression of selected markers
plot.markers <- plotExp(object = seurat_integrated, features = df.markers$GeneID)

# Save marker dotplot
ggsave(plot.markers, file = "Human_epithelial_markers_dotplot.pdf", width=14, height=8)


############################################################################
##           MILO DIFFERNIAL CELL ABUDNACE FOR EPITHELIAL CELLS 
############################################################################

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

# Ensure unique rows by selecting distinct values and setting row names
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

## Plot single-cell UMAP coloured by Smoking status, with annotations
umap_pl <- plotReducedDim(epith_milo, dimred = "UMAPAFTERHARMONY_RNA", colour_by = "Smoking_status", 
                          text_by = "Cluster_Annotations", 
                          text_size = 3, point_size = 0.5) + 
  guides(fill = "none")

## Plot neighbourhood graph coloured by differential analysis results
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

#############################################################################
##           SLINGSHOT CELL LINEAGE INFERENCE FOR EPITHELIAL CELLS 
############################################################################

set.seed(123)

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


#############################################################################
##     Signature Overlap with Human Epithelial from previous publications
#############################################################################

# Identfing Human epithelail cell markers from our current study 
#markers <- FindAllMarkers(object = seurat_integrated, 
#                          only.pos = TRUE,
#                          min.pct = 0.25,
#		     			            logfc.threshold = 0.25)

markers <- fread("data/scRNAseq/Human.finder.markers.by.celltype.tsv")

markers <- markers %>%
           arrange(p_val_adj)  # Sort by significance

# Select top 50 markers per cluster
top_n_genes <- 50
de_markers <- markers %>%
              group_by(cluster) %>%
              slice_head(n = top_n_genes) %>%
              ungroup()

# Read gene signatures from previous studies
dictionary <- fread("data/scRNAseq/Human.airway.cell.types.signatures.from.publications.tsv") %>%
              rename(ref = V1, disease = V2, term = V3, gene = V4)

# Format term names for uniqueness
dictionary <- dictionary %>%
              mutate(term = paste0(ref, "_", term))

# Filter for epithelial-related terms
key_disease <- "Epithelial"
genesets <- dictionary %>%
            filter(disease == key_disease) %>%
            group_by(term) %>%
            slice_head(n = top_n_genes) %>%
            ungroup() %>%
            split(.$term, .$gene)

# Prepare cluster-specific gene lists
genes_by_cluster <- split(de_markers$gene, de_markers$cluster)

# Overrepresentation analysis using Fisher's exact test
fisher_results <- list()
for (cluster in names(genes_by_cluster)) {
  de_genes <- genes_by_cluster[[cluster]]
  bg_genes <- markers$gene[markers$cluster == cluster]
  bg_genes <- setdiff(bg_genes, de_genes)

  cluster_results <- data.frame(
    Geneset = names(genesets),
    Cluster = cluster,
    p.value = NA,
    overlap = NA
  )

  for (geneset in names(genesets)) {
    A <- sum(de_genes %in% genesets[[geneset]])
    B <- sum(bg_genes %in% genesets[[geneset]])
    C <- sum(!de_genes %in% genesets[[geneset]])
    D <- sum(!bg_genes %in% genesets[[geneset]])

    contingency_table <- matrix(c(A, B, C, D), nrow = 2,
                                dimnames = list(c("DE", "Not.DE"),
                                                c("In.gene.set", "Not.in.gene.set")))

    cluster_results$p.value[cluster_results$Geneset == geneset] <- fisher.test(contingency_table, alternative = "greater")$p.value
    cluster_results$overlap[cluster_results$Geneset == geneset] <- round((A / (A + B + C)) * 100, 3)
  }

  cluster_results$p.adjust <- p.adjust(cluster_results$p.value, method = "BH")
  fisher_results[[cluster]] <- cluster_results
}

# Save results
fisher_data <- bind_rows(fisher_results)
fwrite(fisher_data, "Human.fisher.exact.test.results.tsv", sep="\t")

# Define cell type and gene set order
geneset_order <- c(
  'Travaglini et al_Proliferating Basal', 'Travaglini et al_Basal',
  'Travaglini et al_Differentiating Basal', 'Travaglini et al_Club',
  'Travaglini et al_Goblet', 'Travaglini et al_Serous',
  'Travaglini et al_Ciliated', 'Travaglini et al_Ionocyte',
  'Travaglini et al_Neuroendocrine', 'Goldfarbmuren et al_Proliferating.basal',
  'Goldfarbmuren et al_Proteasomal.basal', 'Goldfarbmuren et al_Differentiating.basal',
  'Goldfarbmuren et al_KRT8.high', 'Goldfarbmuren et al_SMG.basal',
  'Goldfarbmuren et al_SMG.secretory', 'Goldfarbmuren et al_Mucus.secretory',
  'Goldfarbmuren et al_Ciliated', 'Goldfarbmuren et al_Ionocyte.tuft',
  'Deprez et al_Cycling Basal', 'Deprez et al_Basal', 'Deprez et al_Suprabasal',
  'Deprez et al_Secretory', 'Deprez et al_Serous', 'Deprez et al_Multiciliated',
  'Sikkema et al_Hillock cells'
)

cell_type_order <- c(
  'Basal cycling', 'Basal 1', 'Basal 2', 'Basal 3', 'Basal SMG', 'KRT13/KRT4',
  'Suprabasal 1', 'Suprabasal 2', 'Suprabasal 3', 'Secretory', 'Serous',
  'Deuterosomal', 'Ciliated', 'Ionocyte', 'Neuroendocrine'
)

# Reorder data factors
fisher_data <- fisher_data %>%
  mutate(Cluster = factor(Cluster, levels = rev(cell_type_order)),
         Geneset = factor(Geneset, levels = geneset_order))

# Define colours for references
ref_colours <- c(
  "Deprez et al" = "navy",
  "Goldfarbmuren et al" = "orange",
  "Travaglini et al" = "purple",
  "Sikkema et al" = "red"
)

ref_color_map <- ref_colours[unique(dictionary$ref)]

# Plot Overlap Results
plot_data <- ggplot(fisher_data, aes(x = Geneset, y = Cluster, fill = -log10(p.adjust + 1e-9))) +
  geom_tile(color = "white", size = 1) +
  scale_fill_gradient2(low = "white", high = "darkred", name = "-log10(p.adjust)") +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
    axis.ticks.x = element_line(color = ref_color_map, size = 5, 
                                arrow = grid::arrow(length = unit(8, "pt"), angle = 90))
  ) +
  labs(x = NULL, y = "Current Study") +
  expand_limits(y = -0.2)

# Add legend
legend_data <- data.frame(
  Geneset = 1, Cluster = 1,
  label = factor(names(ref_colours), levels = names(ref_colours))
)

plot_with_legend <- plot_data +
  geom_point(data = legend_data, aes(x = Geneset, y = Cluster, color = label), inherit.aes = FALSE) +
  scale_colour_manual(values = ref_colours, name = "Studies") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  coord_cartesian(clip = "off")

# Save plot
ggsave(plot_with_legend, file = "Human.Airway.signatures.overlap.pdf", width = 16, height = 10)


#############################################################################
##          Visualize KRT4/KRT13 cell signature
#############################################################################

set.seed(123)

# Load the KRT4/KRT13 signature
table.path <- "data/scRNAseq/Human.finder.markers.by.celltype.tsv"
krt13.smokers <- fread(table.path) %>% filter(cluster == "KRT13/KRT4") %>% head(50)
signature <- unique(krt13.smokers$gene)

# Load Galon data from DOI: 10.1038/s41586-019-1330-0
load("data/scRNAseq/galon.data.RData")

pheno <- jg.pheno
data <- jg.data[, pheno$colname]
pheno <- pheno %>% select(colname, name, group) %>% 
mutate(name = factor(name), group = factor(group)) %>%
column_to_rownames('colname')

# Normalise Galon data using DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(data), colData = pheno, design = ~ group)
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
norm.counts <- counts(dds, normalized=TRUE)

# Transform normalized counts to long format
norm.count.df <- melt(norm.counts)
colnames(norm.count.df) <- c("Gene", "colname", "norm.exp")
norm.count.df <- left_join(norm.count.df, jg.pheno, by="colname")

# Define histological conditions
histos <- c("normal normofluorescent", "hyperplasia", "metaplasia", "mild dysplasia", "moderate dysplasia",
            "severe dysplasia", "carcinoma in situ", "squamous cell carcinoma")
histos.upletter <- c("Normal", "Hyperplasia", "Metaplasia", "Mild dysplasia", "Moderate dysplasia",
                     "Severe dysplasia", "Carcinoma in situ", "Squamous cell carcinoma")

# Filter and preprocess signature data
signature.df <- norm.count.df %>% filter(Gene %in% signature, norm.exp < 30, histology %in% histos) %>%
  mutate(histology = gsub(" normofluorescent", "", histology)) %>%
  mutate(histology = firstup(histology)) %>%
  arrange(match(histology, histos))
rownames(signature.df) <- NULL

# Define developmental stages
signature.df <- signature.df %>% mutate(devstages = case_when(
  histology %in% c("Normal") ~ "Normal",
  histology %in% c("Hyperplasia", "Metaplasia", "Mild dysplasia", "Moderate dysplasia") ~ "Low grade",
  histology %in% c("Severe dysplasia", "Carcinoma in situ") ~ "High grade",
  TRUE ~ "LUSC"
)) %>%
  mutate(devstages = factor(devstages), group="signature-SDC_NCL_ITGA_ITGB")

# Summarise data for visualization
signature.df2 <- signature.df %>% drop_na(Gene) %>% group_by(histology, ID) %>%
  summarise(mean = mean(norm.exp), sd = std.error(norm.exp)) %>%
  arrange(match(histology, histos.upletter))

# Perform statistical comparisons
anno_df <- compare_means(mean ~ histology, ref.group = "Normal", data = signature.df2)

# Define p-value annotations
max.exp <- 3
anno_df <- anno_df %>% mutate(comparison = paste0(group1, ":", group2)) %>%
  select(comparison, p.adj) %>%
  filter(comparison %in% c("Normal:Hyperplasia", "Normal:Metaplasia", "Normal:Mild dysplasia",
                           "Normal:Moderate dysplasia", "Normal:Severe dysplasia", "Normal:Carcinoma in situ",
                           "Normal:Squamous cell carcinoma")) %>%
  mutate(start = "Normal", 
         end = c("Hyperplasia", "Metaplasia", "Mild dysplasia", "Moderate dysplasia", 
                 "Severe dysplasia", "Carcinoma in situ", "Squamous cell carcinoma"),
         y = seq(max.exp + 0.5, max.exp + 3.5, by = 0.5))

# Generate plot
plot <- signature.df2 %>%
  mutate(histology = factor(histology, levels=histos.upletter)) %>%
  ggplot(aes(x=histology, y=mean)) +
  geom_boxplot() +
  geom_jitter(aes(fill=histology), width=0.2, size=3, color="black", pch=21) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  ylab(expression("KRT13+ cells signature")) + xlab("") +
  cowplot::theme_cowplot(font_size = 20) +
  theme(axis.text.x = element_text(size=20, color="black", angle=60, vjust=1, hjust=1),
        legend.position="none") +
  geom_signif(data = anno_df,
              aes(xmin = start, xmax = end, annotations = p.adj, y_position = y),
              textsize = 6, manual = TRUE, inherit.aes = FALSE) +
  ylim(NA, max.exp + 3.5)

# Save plot
ggsave(plot, file="KRT13.signature.score.normal.Top50genes.pdf", width=7, height=8)
