################################################################################################################################################
##                                                                                                                      
##  MOUSE NTCU MODEL SINGLE-CELL GEX DATA VISULIZATION - FIGURE 3 & SUPP FIGURE 4
##                                                                                                                      
##  Date: 03 Jan 2024                                                                                                                   
##  
##  Author: Ahmed Alhendi                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

##############################################################################
##                         LOAD REQUIRED PACKAGES
##############################################################################

# List of required packages
required_packages <- c("Seurat", "data.table", "monocle3", "monocle", "SeuratData", 
                       "SeuratWrappers", "ggplot2", "patchwork", "magrittr", 
                       "dplyr", "scales", "tidyverse", "ggpubr", "cowplot", 
                       "viridis", "pheatmap", "reshape2", "RColorBrewer")


# Load all required packages
suppressMessages({
  lapply(required_packages, library, character.only = TRUE)
})


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


############################################################################
##           MONOCLE TRAJECTORY OF EPITHELIAL CELLS 
############################################################################
set.seed(123)

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


#############################################################################
##    Signature Overlap with Mouse Epithelial from previous publications
#############################################################################

# Identfing Mouse epithelial cell markers from our current study 
#markers <- FindAllMarkers(object = seurat_integrated, 
#                          only.pos = TRUE,
#                          min.pct = 0.25,
#		     			            logfc.threshold = 0.25)

# Read the Mouse epithelial cell markers from our current studys
markers <- fread("data/scRNAseq/NTCU-treated.finder.markers.by.celltype.tsv")

markers <- markers %>%
           arrange(p_val_adj)  # Sort by significance

# Select top 50 markers per cluster
top_n_genes <- 50
de_markers <- markers %>%
              group_by(cluster) %>%
              slice_head(n = top_n_genes) %>%
              ungroup()

# Read gene signatures from previous studies
dictionary <- fread("data/scRNAseq/Mouse.epithelial.signatures.from.publications.tsv") %>%
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
fwrite(fisher_data, "Mouse.fisher.exact.test.results.tsv", sep="\t")


# Order factors
cells.type.order <- c('Basal proliferative','Basal','Basal Tgm2+','Basal Mecom+','Basal regenerative','Krt4/13+ cells','Secretory',
                   'Secretory Mecom+', 'Deuterosomal cells', 'Ciliated Cells', 'Neuroendocrine', 'lonocytes', 'Tuft')


term2ref <- unique(dict[, c("term", "ref")])$ref
names(term2ref) <- unique(dict[, c("term", "ref")])$term

geneset.orders <- c('Plasschaert et al_Cycling Basal',
            'Plasschaert et al_Basal',
            'Plasschaert et al_Krt4/13+',
            'Plasschaert et al_Secretory',
            'Plasschaert et al_Ciliated',
            'Plasschaert et al_PNEC',
            'Plasschaert et al_Ionocytes',
            'Plasschaert et al_Brush',
            'Montoro et al_Basal',
            'Montoro et al_Krt4/Krt13',
            'Montoro et al_Club',
            'Montoro et al_Ciliated',
            'Montoro et al_Neuroendocrine',
            'Montoro et al_Ionocyte',
            'Montoro et al_Tuft'
           )

data <- data %>% mutate(Cluster = factor(Cluster, levels=rev(cells.type.order)),
                       Geneset = factor(Geneset, levels=geneset.orders))

# Select colors
cols <- c(
 "Plasschaert et al" = "navy",
  "Montoro et al" = "red"
)


ref.color <- cols[term2ref[levels(data$Geneset)]]


# Plot Overlap Results

# Create a dummy dataset for the additional legend
dummy_data <- data.frame(
  Geneset = 1, # Dummy x-axis position
  Cluster = 1, # Dummy y-axis position
  label = factor(
    names(cols),
    levels = c("Plasschaert et al",
               "Montoro et al"
               )
  ) # Specify the desired order
)

# Original base plot
base_plot <- ggplot(
  data = data,
  mapping = aes(x = Geneset, y = Cluster, fill = -log10(p.adjust + 0.000000001))
) +
  geom_tile(col = "white", size = 1) +
  scale_fill_gradient2(low = "white", high = "darkred", 
#                           breaks = c(1, 3, 5),          # Specify custom breaks
                       name = "-log10(p.adjust)") +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
    axis.ticks.x = element_line(
      color = ref.color, size = 5, arrow = grid::arrow(
        length = grid::unit(8, "pt"), angle = 90
      )
    )
  ) +
  labs(x = NULL, y = "Current study")+
expand_limits(y = -0.2)


# Add the additional legend without affecting axis scales
plot_with_legend <- base_plot +
  # Add dummy data for legend
  geom_point(
    data = dummy_data,
    aes(x = Geneset, y = Cluster, colour = label), # Use fixed positions
    inherit.aes = FALSE) +
  scale_colour_manual(
    values = cols,
    name = "Studies") +
  guides(
    colour = guide_legend(override.aes = list(size = 5)) # Adjust legend icon size
  ) +
  coord_cartesian(clip = "off") # Prevent cropping of dummy points outside visible area


# Save plot
ggsave(plot_with_legend, file = "Mouse.NTCU.signatures.overlap.pdf", width = 16, height = 10)

