################################################################################################################################################
##                                                                                                                      
##  HUMAN WES data DATA - Generate SAMPLE-SAMPLE correlation
##                                                                                                                      
##  Date: 18 july 2024                                                                                                                   
##  
##  Author: Ahmed Alhendi                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

# List of required packages
library(dplyr)
library(data.table)
library(tibble)
library(stringr)
library(maftools)
library(circlize)
library(ggplot2)
library(cowplot)
library(ggsci)
library(ComplexHeatmap)

# Define file paths
sample_info_file <- "data/WES/summar.plot.annotations.tsv"
path_bams <- "/path/to/bams/"

# Load sample information
sample_info <- fread(sample_info_file)
sample_info$bams <- paste0(path_bams, sample_info$Tumor_Sample_Barcode, ".dedup.BQSR.sorted")

# Mapping for sample names
mapping <- sample_info %>% select(bams, SAMPLE)

# Perform sample swaps correlation analysis
res <- maftools::sampleSwaps(bams = mapping$bams, build = "hg19", prefix="chr", add=FALSE, ncores = 10)

# Save results
write.table(res$pairwise_comparison, "Multi_region.sampleSwaps.corr.txt", sep="\t", quote=FALSE)
write.csv(res$AF_table, "Multi_region.sampleSwaps.AFtable.csv")

# Compute correlation matrix
cor_table <- cor(res$AF_table, use = "pairwise.complete.obs")

# Define colour scheme
colors <- colorRamp2(seq(0, 1, by = 0.01), 
                     colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(101))

# Create heatmap 
ht <- Heatmap(cor_table,
              col = colors,
              name = "Correlation",
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              column_title = "Correlation Heatmap",
              column_names_rot = 90,
              heatmap_legend_param = list(
                title = "Pearson\ncorrelation R",
                at = seq(0, 1, 0.2),
                labels = seq(0, 1, 0.2),
                legend_height = unit(6, "cm"),
                title_gp = gpar(fontsize = 12, fontface = "bold"),
                labels_gp = gpar(fontsize = 10)
              ),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", cor_table[i, j]), x, y, 
                          gp = gpar(fontsize = 10, col = "white"))
              }
)

# Save the heatmap
pdf("Multi_region.sampleSwaps.Pearson.corr.pdf", width=8, height=8)
draw(ht)
dev.off()