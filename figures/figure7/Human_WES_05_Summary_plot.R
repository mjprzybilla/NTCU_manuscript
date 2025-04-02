################################################################################################################################################
##                                                                                                                      
##  HUMAN WES data DATA - Generate SUMMARY PLOT 
##                                                                                                                      
##  Date: 14 August 2024                                                                                                                   
##  
##  Author: Ahmed Alhendi                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

# Load required libraries
library(dplyr)
library(data.table)
library(tibble)
library(stringr)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer) 

# Load annotation file
anno <- fread("data/WES/summar.plot.annotations.tsv")

# Prepare mutation data for heatmap
tmb.df <- anno %>% select(SAMPLE, Truncal.SNV, Subclonal.SNV) %>% 
  column_to_rownames(var="SAMPLE")
mutations <- list(tmb.df)

# Define color palettes
muts_anno_col <- c("#377EB8", "#9E2536")
hist_col <- c("MoD" = "#FFF5F0", "CIS" = "#FC9272")
smoking_col <- c("Ex-Smoker" = "grey", "Current Smoker" = "black")

# Define site colors
sites <- c("Trachea", "RIB", "RMB", "RIBRUL", "RUL", "RLL", "LMB", "LUL", "LULLLL")
site_col <- c(brewer.pal(8, "Accent"), "purple")
names(site_col) <- sites

# Define column annotations
column_ha <- HeatmapAnnotation(
  Mutations = anno_barplot(mutations, gp = gpar(fill = muts_anno_col), 
                           border = TRUE, bar_width = 1, height = unit(4, "cm")),
  annotation_name_side = "left",
  annotation_name_align = FALSE,
  show_annotation_name = TRUE,
  border = TRUE,
  show_legend = FALSE
)

column_ha2 <- HeatmapAnnotation(
  Smoking_status = as.vector(anno$Smoking_status),
  Histology = as.vector(anno$Histology),
  Site = as.vector(anno$Site),
  col = list(Site = site_col, Smoking_status = smoking_col, Histology = hist_col),
  annotation_name_side = "left",
  annotation_name_align = FALSE,
  show_annotation_name = TRUE,
  border = TRUE,
  show_legend = FALSE
)

# Define legends
lgd_list <- list(
  Legend(labels = c("Truncal", "Subclonal"), title = "Clonality", type = "grid",
         legend_gp = gpar(fill = muts_anno_col)),
  Legend(labels = c("Ex-Smoker", "Current Smoker"), title = "Smoking status", type = "grid",
         legend_gp = gpar(fill = c("grey", "black"))),
  Legend(labels = c("MoD", "CIS"), title = "Histology", type = "grid",
         legend_gp = gpar(fill = c("#FFF5F0", "#FC9272"))),
  Legend(labels = sites, title = "Sites", type = "grid",
         legend_gp = gpar(fill = site_col), nrow = 2)
)

# Prepare matrix for heatmap
mat <- anno %>% select(SAMPLE) %>% 
  column_to_rownames(var="SAMPLE") %>% 
  as.matrix() %>% t()


# Create heatmap
ht1 <- Heatmap(
  mat, 
  top_annotation = column_ha,
  bottom_annotation = column_ha2,
  column_order = seq_along(colnames(mat)),
  column_split = factor(anno$PID, levels = unique(anno$PID)),
  column_title = NULL,
  show_column_names = TRUE,
  column_names_side = "bottom"
)

# Draw and save heatmap
pdf("Sample_summary.human.CIS.plot.pdf", width = 6, height = 7)
draw(ht1, merge_legend = TRUE, heatmap_legend_side = "top", annotation_legend_list = lgd_list)
dev.off()