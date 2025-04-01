################################################################################################################################################
##                                                                                                                      
##  Human SINGLE-CELL GEX data - UMAP & Cell marker visualizations
##                                                                                                                      
##  Date: 05 Dec 2024                                                                                                                   
##  
##  Author: Ahmed Alhendi                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(data.table)
library(DESeq2)
library(reshape2)
library(ggsignif)
library(cowplot)
library(ggpubr)
library(tidyr)
library(RColorBrewer)
library(patchwork)
library(stringr)

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
