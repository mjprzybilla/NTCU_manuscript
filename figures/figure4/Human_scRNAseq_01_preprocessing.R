
################################################################################################################################################
##                                                                                                                      
##  Human SINGLE-CELL GEX data pre-processing
##                                                                                                                      
##  Date: 18 July 2024                                                                                                                   
##  
##  Author: Ahmed Alhendi                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

# Load all required packages

library(Seurat)
library(dplyr)
library(harmony)
library(data.table)
library(ggplot2)

# Function to load and label datasets
load_and_label_seurat <- function(file_path, dataset_name) {
  obj <- readRDS(file_path)
  obj$dataset <- dataset_name
  return(obj)
}

# Load datasets
seurat_list <- list(
  load_and_label_seurat("Deprez.never.smoker.tracheal.grch38.rds", "Deprez"),
  load_and_label_seurat("Goldfarbmuren.cluster.markers.harmony.cell.types.12July2024.rds", "Goldfarbmuren"),
  load_and_label_seurat("maral.treacheal.10PCs.Harmony.8July.v2.rds", "Maral")
)

# Merge all datasets
seurat.obj.filtered <- Reduce(merge, seurat_list)
seurat.obj.filtered <- JoinLayers(seurat.obj.filtered)

# Preprocess Seurat object
seurat.obj.filtered <- NormalizeData(seurat.obj.filtered) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Filter variable genes
genes <- VariableFeatures(seurat.obj.filtered)
genes <- genes[!grepl("^MT|^RP[SL]|^RPL", genes)]
VariableFeatures(seurat.obj.filtered) <- genes

# Scale and run PCA
seurat.obj.filtered <- seurat.obj.filtered %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(reduction.name = 'pca', verbose = FALSE)

# Run Harmony integration
seurat.obj.filtered <- RunHarmony(seurat.obj.filtered, group.by.vars = c("Smoking_status", "dataset"),
                                  assay.use = "RNA", reduction = "pca",
                                  reduction.save = "harmony_RNA", theta = c(1, 1))

# Run UMAP and clustering
seurat.obj.filtered <- RunUMAP(seurat.obj.filtered, reduction = "harmony_RNA", dims = 1:10, 
                               reduction.name = "umapAfterHarmony_RNA", reduction.key='umapAfterHarmony_RNA_')
seurat.obj.filtered <- FindNeighbors(seurat.obj.filtered, dims = 1:10, reduction = "harmony_RNA")
seurat.obj.filtered <- FindClusters(seurat.obj.filtered, algorithm = 1, resolution = c(0.6, 0.8))

# Find markers
Idents(seurat.obj.filtered) <- "RNA_snn_res.0.6"
markers <- FindAllMarkers(seurat.obj.filtered, only.pos = TRUE)

# Select top 5 markers per cluster
markers_filtered <- markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup()

# Save markers
fwrite(markers, "Human.finder.markers.by.clusterid.tsv", sep="\t")

# Annotate clusters based on identified markers
cluster_annotations <- c("0" = "Basal 1",
                         "1" = "Suprabasal 2",
                         "2" = "Basal 3",
                         "3" = "Secretory",
                         "4" = "Suprabasal 1",
                         "5" = "Basal 2",
                         "6" = "KRT13/KRT4",
                         "7" = "Ciliated",
                         "8" = "Basal cycling",
                         "9" = "Serous",
                         "10" = "Suprabasal 3",
                         "11" = "Basal SMG",
                         "12" = "Deuterosomal",
                         "13" = "Ionocyte",
                         "14" = "Neuroendocrine"
                        )

# Assign annotations
annotated_clusters <- cluster_annotations[as.character(Idents(seurat.obj.filtered))]
names(annotated_clusters) <- names(Idents(seurat.obj.filtered))
seurat.obj.filtered <- AddMetaData(seurat.obj.filtered, metadata = annotated_clusters, col.name = "Cluster_Annotations")


# Save final Seurat object
saveRDS(seurat.obj.filtered, "Human_tracheal_epithelial_harmony_seurat.rds")
