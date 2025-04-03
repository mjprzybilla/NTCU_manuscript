################################################################################################################################################
##                                                                                                                      
##  MOUSE NTCU MODEL SINGLE-CELL GEX DATA - EPITHELIAL SUBCLUSTERING
##                                                                                                                      
##  Date: 03 Jan 2024                                                                                                                   
##  
##  Author: Ahmed Alhendi                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

# Load all required packages

library(Seurat)
library(dplyr)
library(data.table)

options(Seurat.object.assay.version = "v3")


# Read in Seurat object
seurat_integrated_all <- readRDS("data/scRNAseq/NTCU_GEX_Trachea_processed_Celltype_seurat.rds")

# Subset to Epithelial compartment
seurat.epithelial <- subset(seurat_integrated_all, Cell_types %in% "Epithelial")

# Split the Seurat object 
split_seurat <- SplitObject(seurat.epithelial, split.by = "phenotype")


for (i in 1:length(split_seurat))
{
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- FindVariableFeatures(split_seurat[[i]], selection.method = "vst", nfeatures = 3000, verbose = FALSE)

    
}
features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = split_seurat,anchor.features = features, reduction = "cca", dims = 1:50)
seurat_integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

# run standard anlaysis workflow
seurat_integrated <- FindVariableFeatures(seurat_integrated)
seurat_integrated <- ScaleData(seurat_integrated)
seurat_integrated <- RunPCA(seurat_integrated)

seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:50)
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:50) %>% 
FindClusters(resolution = c(0.3,0.4,0.6))



# Find markers
Idents(seurat_integrated) <- "RNA_snn_res.0.4"

markers <- FindAllMarkers(object = seurat_integrated, 
                          only.pos = TRUE,
                          min.pct = 0.25,
		     			  logfc.threshold = 0.25)

# Save markers
fwrite(markers, "NTCU-treated.finder.markers.by.clusteris.tsv", sep="\t")


# Annotate clusters based on identified markers
cluster_annotations <- c("0" = "Basal Krt14+",
                         "1" = "Basal proliferative",
                         "2" = "Krt4/13+",
                         "3" = "Secretory",
                         "4" = "Secretory Mecom+",
                         "5" = "Basal Tgm2+",
                         "6" = "Basal",
                         "7" = "Basal Mecom+",
                         "8" = "Deuterosomal cells",
                         "9" = "Ciliated Cells",
                         "10" = "Ionocytes",
                         "11" = "Tuft",
                         "12" = "Neuroendocrine"
                        )


# Assign annotations
annotated_clusters <- cluster_annotations[as.character(Idents(seurat_integrated))]
names(annotated_clusters) <- names(Idents(seurat_integrated))
seurat_integrated <- AddMetaData(seurat_integrated, metadata = annotated_clusters, col.name = "celltypes")


# Filter by age (keep only 15 weeks)
seurat.obj.filtered <- subset(seurat_integrated, age %in% "15 weeks")

# Save final Seurat object
saveRDS(seurat.obj.filtered, "NTCU_Trachea_Epithelial_integrated_seurat.rds")
