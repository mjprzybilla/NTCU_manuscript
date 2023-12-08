################################################################################################################################################
##                                                                                                                      
##  COMPARE DIFFERENTIALLY EXPRESSED GENES ACROSS NTCU VS CONTROL IN CELL TYPES
##                                                                                                                      
##  Date: 18 JULY 2023                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

# clear workspace beforehand
rm(list = ls())

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("Matrix", "tidyverse", "stringr", "biomaRt", "data.table", "corrplot", "patchwork", "ComplexHeatmap", "gridExtra", "tidyr", 
                      "BoutrosLab.plotting.general")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

############################################################################
##                    READ IN AND VISUALISE DEG DATA
############################################################################

# list all the cell type results
deg.files <- list.files("/Users/mp34/sanger/team154pc/mp34/lung/NTCU/scRNA", pattern = "NTCU_DESeq2", all.files = T, full.names = T)
deg.files <- deg.files[-grep("/._|Condition|Control_vs_NTCU|Proliferative", deg.files)]
deg.data <- lapply(deg.files, read.table, header = T, sep = ",")

# get the cell types
celltype.ids <- str_split_fixed(basename(deg.files), "_", 5)[,3]

i <- 1
deg.up.list <- list()
deg.down.list <- list()

# iterate over dataframes and subset to significant genes
for (i in 1:length(celltype.ids)){
  
  data.tmp <- deg.data[[i]]
  data.tmp <- data.tmp[data.tmp$padj < 0.1 & abs(data.tmp$log2FoldChange) > 0.5, ]
  data.tmp <- data.tmp[!is.na(data.tmp$X),]
  data.tmp$celltype <- celltype.ids[i]
  
  deg.down.list[[celltype.ids[i]]] <- data.tmp[data.tmp$log2FoldChange < 0, "X"]
  deg.up.list[[celltype.ids[i]]] <- data.tmp[data.tmp$log2FoldChange > 0, "X"]
  
}

# create comparison matrix
comp.up.mtx <- make_comb_mat(deg.up.list, mode = "intersect")
comp.down.mtx <- make_comb_mat(deg.down.list, mode = "intersect")

# filter matrices
comp.up.mtx <- comp.up.mtx[comb_degree(comp.up.mtx) > 1]
comp.down.mtx <- comp.down.mtx[comb_degree(comp.down.mtx) > 1]

# visualise up-regulated genes
pdf("/Users/mp34/team154_campbell/plots/NTCU/shared_DEG_up_UpsetPlot_wo_proliferative.pdf", width = 5, height = 3)
UpSet(comp.up.mtx, comb_order = order(comb_size(comp.up.mtx)), comb_col = "darkmagenta", bg_col = "lightgrey", bg_pt_col = "white", 
      right_annotation = upset_right_annotation(comp.up.mtx, gp = gpar(fill = "darkmagenta", color = "darkmagenta")))
dev.off()

# create core, semi-unique and unique gene lists
core.genes <- extract_comb(comp.up.mtx, "111")
semi.unique.genes <- unique(c(extract_comb(comp.up.mtx, "110"), extract_comb(comp.up.mtx, "101"), extract_comb(comp.up.mtx, "011")))
semi.unique.genes <- semi.unique.genes[-grep(paste0(core.genes, collapse = "|"), semi.unique.genes)]

core.up.df <- data.frame("gene" = core.genes, "direction" = "up", "status" = "core")
semi.unique.up.df <- data.frame("gene" = semi.unique.genes, "direction" = "up", "status" = "semi-unique")

# visualise down-regulated genes
pdf("/Users/mp34/team154_campbell/plots/NTCU/shared_DEG_down_UpsetPlot_wo_proliferative.pdf", width = 5, height = 3)
UpSet(comp.down.mtx, comb_order = order(comb_size(comp.down.mtx)), comb_col = "darkblue", bg_col = "lightgrey", bg_pt_col = "white", 
      right_annotation = upset_right_annotation(comp.down.mtx, gp = gpar(fill = "darkblue", color = "darkblue")))
dev.off()

# create core, semi-unique and unique gene lists
core.genes <- extract_comb(comp.down.mtx, "111")
semi.unique.genes <- unique(c(extract_comb(comp.down.mtx, "110"), extract_comb(comp.down.mtx, "101"), extract_comb(comp.down.mtx, "011")))
semi.unique.genes <- semi.unique.genes[-grep(paste0(core.genes, collapse = "|"), semi.unique.genes)]

core.down.df <- data.frame("gene" = core.genes, "direction" = "down", "status" = "core")
semi.unique.down.df <- data.frame("gene" = semi.unique.genes, "direction" = "down", "status" = "semi-unique")

# combine to overall dataframe
all.deg.df <- rbind(core.up.df, semi.unique.up.df, core.down.df, semi.unique.down.df)
write.table(all.deg.df, "/Users/mp34/team154_campbell/data/NTCU_core_DEG_wo_proliferative.txt", col.names = T, row.names = F, sep = "\t", quote = F)

############################################################################
##            VISUALISE VOLCANO PLOT FOR CLUB/SECRETORY CELLS
############################################################################

secr.deg.data <- deg.data[[3]]

# make dataframe and order according to log2FC and pvalue
colnames(secr.deg.data)[1] <- "hgnc_symbol" 
secr.deg.data <- secr.deg.data[!is.na(secr.deg.data$hgnc_symbol),]
secr.deg.data <- secr.deg.data[complete.cases(secr.deg.data),]
secr.deg.data <- secr.deg.data[order(secr.deg.data$padj, decreasing =  F),]

# 
secr.deg.data$color <- "NS"
secr.deg.data[abs(secr.deg.data$log2FoldChange) > 0.5, "color"] <- "Log2FC"
secr.deg.data[secr.deg.data$padj < 0.01, "color"] <- "p-value"
secr.deg.data[abs(secr.deg.data$log2FoldChange) > 0.5 & secr.deg.data$padj < 0.01, "color"] <- "Log2FC & p-value"
secr.deg.data$color <- factor(secr.deg.data$color, levels = c("NS", "Log2FC", "p-value", "Log2FC & p-value"))

col.values <- c("lightgrey", pal_lancet("lanonc")(3))
col.values <- c("lightgrey", "darkblue", "darkmagenta", "darkgreen")
options(ggrepel.max.overlaps = Inf)
p1 <- ggplot(secr.deg.data, aes(x=log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = color), size = 1) +
  geom_vline(xintercept=c(-0.5, 0.5), linetype="dotted", color = "black") +
  geom_hline(yintercept = -log10(0.01), linetype="dotted", color = "black") +
  theme_classic() +
  # geom_smooth(method=lm, color="darkblue") +
  labs(x="Log2FC",
       y="-Log10(p-value)") + 
  xlim(-8, 8) +
  theme(legend.position="bottom") +
  scale_color_manual(values = col.values) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.75),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold")) +
  geom_text_repel(data = secr.deg.data %>% filter(padj < 10e-25 , log2FoldChange <= -0.5),
                  aes(label = paste0(hgnc_symbol)) ,
                  hjust = -.35,
                  nudge_x = -0.75,
                  direction = "y",
                  fontface = "bold",
                  size = 2) + 
  geom_text_repel(data = secr.deg.data %>% filter(padj < 10e-25 , log2FoldChange >= 0.5),
                  aes(label = paste0(hgnc_symbol)) ,
                  hjust = -.35,
                  nudge_x = 0.5,
                  direction = "y",
                  fontface = "bold",
                  size = 2)

pdf("/Users/mp34/team154_campbell/plots/NTCU/SecretoryCells_volcano_DESeq2.pdf", width = 5, height = 4)
print(p1)
dev.off()


############################################################################
##            VISUALISE VOLCANO PLOT FOR BASAL CELLS
############################################################################

basal.deg.data <- deg.data[[1]]

# make dataframe and order according to log2FC and pvalue
colnames(basal.deg.data)[1] <- "hgnc_symbol" 
basal.deg.data <- basal.deg.data[!is.na(basal.deg.data$hgnc_symbol),]
basal.deg.data <- basal.deg.data[complete.cases(basal.deg.data),]
basal.deg.data <- basal.deg.data[order(basal.deg.data$padj, decreasing =  F),]

# 
basal.deg.data$color <- "NS"
basal.deg.data[abs(basal.deg.data$log2FoldChange) > 0.5, "color"] <- "Log2FC"
basal.deg.data[basal.deg.data$padj < 0.01, "color"] <- "p-value"
basal.deg.data[abs(basal.deg.data$log2FoldChange) > 0.5 & basal.deg.data$padj < 0.01, "color"] <- "Log2FC & p-value"
basal.deg.data$color <- factor(basal.deg.data$color, levels = c("NS", "Log2FC", "p-value", "Log2FC & p-value"))

col.values <- c("lightgrey", pal_lancet("lanonc")(3))
col.values <- c("lightgrey", "darkblue", "darkmagenta", "darkgreen")
options(ggrepel.max.overlaps = Inf)
p1 <- ggplot(basal.deg.data, aes(x=log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = color), size = 1) +
  geom_vline(xintercept=c(-0.5, 0.5), linetype="dotted", color = "black") +
  geom_hline(yintercept = -log10(0.01), linetype="dotted", color = "black") +
  theme_classic() +
  # geom_smooth(method=lm, color="darkblue") +
  labs(x="Log2FC",
       y="-Log10(p-value)") + 
  xlim(-7, 7) +
  theme(legend.position="bottom") +
  scale_color_manual(values = col.values) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.75),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold")) +
  geom_text_repel(data = basal.deg.data %>% filter(padj < 10e-35 , log2FoldChange <= -0.5),
                  aes(label = paste0(hgnc_symbol)) ,
                  hjust = -.35,
                  nudge_x = -0.75,
                  direction = "y",
                  fontface = "bold",
                  size = 2) + 
  geom_text_repel(data = basal.deg.data %>% filter(padj < 10e-35 , log2FoldChange >= 0.5),
                  aes(label = paste0(hgnc_symbol)) ,
                  hjust = -.35,
                  nudge_x = 0.5,
                  direction = "y",
                  fontface = "bold",
                  size = 2)

pdf("/Users/mp34/team154_campbell/plots/NTCU/BasalCells_volcano_DESeq2.pdf", width = 5, height = 4)
print(p1)
dev.off()


############################################################################
##            VISUALISE VOLCANO PLOT FOR CILIATED CELLS
############################################################################

ciliated.deg.data <- deg.data[[2]]

# make dataframe and order according to log2FC and pvalue
colnames(ciliated.deg.data)[1] <- "hgnc_symbol" 
ciliated.deg.data <- ciliated.deg.data[!is.na(ciliated.deg.data$hgnc_symbol),]
ciliated.deg.data <- ciliated.deg.data[complete.cases(ciliated.deg.data),]
ciliated.deg.data <- ciliated.deg.data[order(ciliated.deg.data$padj, decreasing =  F),]

# 
ciliated.deg.data$color <- "NS"
ciliated.deg.data[abs(ciliated.deg.data$log2FoldChange) > 0.5, "color"] <- "Log2FC"
ciliated.deg.data[ciliated.deg.data$padj < 0.01, "color"] <- "p-value"
ciliated.deg.data[abs(ciliated.deg.data$log2FoldChange) > 0.5 & ciliated.deg.data$padj < 0.01, "color"] <- "Log2FC & p-value"
ciliated.deg.data$color <- factor(ciliated.deg.data$color, levels = c("NS", "Log2FC", "p-value", "Log2FC & p-value"))

col.values <- c("lightgrey", pal_lancet("lanonc")(3))
col.values <- c("lightgrey", "darkblue", "darkmagenta", "darkgreen")
options(ggrepel.max.overlaps = Inf)
p1 <- ggplot(ciliated.deg.data, aes(x=log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = color), size = 1) +
  geom_vline(xintercept=c(-0.5, 0.5), linetype="dotted", color = "black") +
  geom_hline(yintercept = -log10(0.01), linetype="dotted", color = "black") +
  theme_classic() +
  # geom_smooth(method=lm, color="darkblue") +
  labs(x="Log2FC",
       y="-Log10(p-value)") + 
  xlim(-5, 5) +
  theme(legend.position="bottom") +
  scale_color_manual(values = col.values) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.75),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold")) +
  geom_text_repel(data = ciliated.deg.data %>% filter(padj < 10e-10, log2FoldChange <= -0.5),
                  aes(label = paste0(hgnc_symbol)) ,
                  hjust = -.35,
                  nudge_x = -0.75,
                  direction = "y",
                  fontface = "bold",
                  size = 2) + 
  geom_text_repel(data = ciliated.deg.data %>% filter(padj < 10e-10, log2FoldChange >= 0.5),
                  aes(label = paste0(hgnc_symbol)) ,
                  hjust = -.35,
                  nudge_x = 0.5,
                  direction = "y",
                  fontface = "bold",
                  size = 2)

pdf("/Users/mp34/team154_campbell/plots/NTCU/CiliatedCells_volcano_DESeq2.pdf", width = 5, height = 4)
print(p1)
dev.off()


############################################################################
##                  READ IN AND VISUALISE COLLECT TRI DATA
############################################################################

# list all the cell type results
collecttri.files <- list.files("/Users/mp34/sanger/team154pc/mp34/lung/NTCU/scRNA", pattern = "collectTRI", all.files = T, full.names = T)
collecttri.files <- collecttri.files[-c(1, 4)]
collecttri.data <- lapply(collecttri.files, read.table, header = T, sep = ",")

# get the cell types
celltype.ids <- str_split_fixed(basename(collecttri.files), "_", 5)[,1]

i <- 1
collecttri.up.list <- list()
collecttri.down.list <- list()

# iterate over dataframes and subset to significant genes
for (i in 1:length(celltype.ids)){
  
  data.tmp <- data.frame(t(collecttri.data[[i]]))
  data.tmp <- data.tmp[-c(1,2), ]
  colnames(data.tmp) <- c("activities", "padj")
  data.tmp$padj <- as.numeric(as.character(data.tmp$padj)) 
  data.tmp$activities <- as.numeric(as.character(data.tmp$activities)) 
  data.tmp <- data.tmp[data.tmp$padj < 0.1, ]
  data.tmp$celltype <- celltype.ids[i]
  
  collecttri.up.list[[celltype.ids[i]]] <- rownames(data.tmp[data.tmp$activities > 0, ])
  collecttri.down.list[[celltype.ids[i]]] <- rownames(data.tmp[data.tmp$activities < 0, ])
  
}

# create comparison matrix
comp.collecttri.up.mtx <- make_comb_mat(collecttri.up.list, mode = "intersect")
comp.collecttri.down.mtx <- make_comb_mat(collecttri.down.list, mode = "intersect")

# filter matrices
comp.collecttri.up.mtx <- comp.collecttri.up.mtx[comb_degree(comp.collecttri.up.mtx) > 1]
comp.collecttri.down.mtx <- comp.collecttri.down.mtx[comb_degree(comp.collecttri.down.mtx) > 1]

# visualise up-regulated genes
pdf("/Users/mp34/team154_campbell/plots/NTCU/shared_TFs_collectTRI_up_UpsetPlot_wo_proliferative.pdf", width = 5, height = 3)
UpSet(comp.collecttri.up.mtx, comb_order = order(comb_size(comp.collecttri.up.mtx)), comb_col = "darkmagenta", bg_col = "lightgrey", bg_pt_col = "white", 
      right_annotation = upset_right_annotation(comp.collecttri.up.mtx, gp = gpar(fill = "darkmagenta", color = "darkmagenta")))
dev.off()

# create core, semi-unique and unique gene lists
core.tfs <- extract_comb(comp.collecttri.up.mtx, "111")
semi.unique.tfs <- unique(c(extract_comb(comp.collecttri.up.mtx, "110"), extract_comb(comp.collecttri.up.mtx, "101"), extract_comb(comp.collecttri.up.mtx, "011")))
semi.unique.tfs <- semi.unique.tfs[-grep(paste0(core.tfs, collapse = "|"), semi.unique.tfs)]

core.up.df <- data.frame("gene" = core.tfs, "direction" = "up", "status" = "core")
semi.unique.up.df <- data.frame("gene" = semi.unique.tfs, "direction" = "up", "status" = "semi-unique")

# visualise down-regulated genes
pdf("/Users/mp34/team154_campbell/plots/NTCU/shared_TFs_collectTRI_down_UpsetPlot_wo_proliferative.pdf", width = 5, height = 3)
UpSet(comp.collecttri.down.mtx, comb_order = order(comb_size(comp.collecttri.down.mtx)), comb_col = "darkblue", bg_col = "lightgrey", bg_pt_col = "white", 
      right_annotation = upset_right_annotation(comp.collecttri.down.mtx, gp = gpar(fill = "darkblue", color = "darkblue")))
dev.off()

# create core, semi-unique and unique gene lists
core.tfs <- extract_comb(comp.collecttri.down.mtx, "111")
semi.unique.tfs <- unique(c(extract_comb(comp.collecttri.down.mtx, "110"), extract_comb(comp.collecttri.down.mtx, "101"), extract_comb(comp.collecttri.down.mtx, "011")))
semi.unique.tfs <- semi.unique.tfs[-grep(paste0(core.tfs, collapse = "|"), semi.unique.tfs)]

core.down.df <- data.frame("gene" = core.tfs, "direction" = "down", "status" = "core")
semi.unique.down.df <- data.frame("gene" = semi.unique.tfs, "direction" = "down", "status" = "semi-unique")

# combine to overall dataframe
all.tfs.df <- rbind(core.up.df, semi.unique.up.df, core.down.df, semi.unique.down.df)
write.table(all.tfs.df, "/Users/mp34/team154_campbell/data/NTCU_core_TFs_collectTRI_wo_proliferative.txt", col.names = T, row.names = F, sep = "\t", quote = F)

############################################################################
##              READ IN AND VISUALISE PROGENY RESULTS
############################################################################

# list all the cell type results
progeny.files <- list.files("/Users/mp34/sanger/team154pc/mp34/lung/NTCU/scRNA", pattern = "progeny", all.files = T, full.names = T)
progeny.files <- progeny.files[-c(1, 4)]
progeny.data <- lapply(progeny.files, read.table, header = T, sep = ",")

# get the cell types
celltype.ids <- str_split_fixed(basename(progeny.files), "_", 5)[,1]

i <- 2
progeny.act.list <- list()
progeny.fdr.list <- list()

# iterate over dataframes and subset to significant genes
for (i in 1:length(celltype.ids)){
  
  # wrangle into shape
  data.tmp <- data.frame(t(progeny.data[[i]]))
  colnames(data.tmp) <- c(paste0(c("activities_", "padj_"), celltype.ids[i]))
  data.tmp$pathways <- rownames(data.tmp)
  data.tmp <- data.tmp[-c(1,2), ]
  data.tmp[,1] <- as.numeric(as.character(data.tmp[,1]))
  data.tmp[,2] <- as.numeric(as.character(data.tmp[,2]))
  
  data.tmp[which(data.tmp[,2] == 0),2] <- min(data.tmp[data.tmp[,2] != 0, 2])
  
  progeny.act.list[[celltype.ids[i]]] <- data.tmp[, c("pathways", paste0("activities_", celltype.ids[i]))]
  progeny.fdr.list[[celltype.ids[i]]] <- data.tmp[, c("pathways", paste0("padj_", celltype.ids[i]))]
    
}

plot_data_acts <- purrr::reduce(progeny.act.list, full_join, by = "pathways")
plot_data_fdr <- purrr::reduce(progeny.fdr.list, full_join, by = "pathways")

rownames(plot_data_acts) <- rownames(plot_data_fdr) <- plot_data_fdr$pathways
plot_data_fdr$pathways <- NULL
plot_data_acts$pathways <- NULL
colnames(plot_data_acts) <- colnames(plot_data_fdr) <- str_split_fixed(colnames(plot_data_acts), "_", 2)[,2]

# set yaxis lab 
yaxis.lab <- rownames(plot_data_acts)

# set spot size and colour functions
spot.size.function <- function(x) {abs(x)/2}
spot.colour.function <- function(x) {
  colours <- rep('white', length(x));
  colours[x > 0] <- 'firebrick3';
  colours[x < 0] <- 'dodgerblue3';
  colours[x == 0] <- 'transparent';
  return(colours);
}

key_sizes <- seq(8, -8, -4);
key_cov <- draw.key(
  list(
    space = 'right',
    points = list(
      cex = spot.size.function(key_sizes),
      col = 'black',
      fill = spot.colour.function(key_sizes),
      pch = 21
    ),
    text = list(
      lab = c('4','2','0','-2','-4'),
      cex = 1.5,
      adj = 1.5,
      fontface = 'bold'
    ),
    title = expression(bold("Observed\nScore")),
    cex.title = 1.8,
    background = 'white',
    padding = 8
  )
)


# create dotmap
create.dotmap(
  x = t(plot_data_acts),
  bg.data = t(-log10(plot_data_fdr)),
  filename = '/Users/mp34/team154_campbell/plots/NTCU/NTCU_Progeny_heatmap_wo_proliferative.pdf',
  xaxis.cex = 1.2,
  xaxis.rot = 45,
  yaxis.tck = 0,
  xaxis.tck = 0,
  spot.colour.function = spot.colour.function,
  spot.size.function = spot.size.function,
  xaxis.lab = paste('   ', yaxis.lab),
  yaxis.lab = paste('', colnames(plot_data_acts)),
  colour.scheme = c('white','black'),
  na.spot.size = 0,
  bg.alpha = 1,
  colourkey = TRUE, 
  colourkey.labels = c(
    expression(1),
    expression(10^-2),
    expression(10^-3), 
    expression(10^-4),
    expression(10^-5),
    expression(10^-6),
    expression(10^-7)
  ),
  at = c(0, seq(1, 6.1, 0.2)),
  colourkey.labels.at = seq(0, 6, 1),
  colourkey.cex = 1.5,
  axis.top = 1.5,
  key = NULL,
  legend = list(
    inside = list(
      fun = key_cov,
      x = 1.05,
      y = 0.95
    ),
    inside = list(
      fun = draw.key(list(text = list(lab = ''), title = expression(bold('FDR')), cex.title = 1.8)),
      x = -0.1,
      y = -0.29
    )
  ),
  width = 13,
  height = 4,
  bottom.padding = 3,
  right.padding = 20,
  left.padding = 25,
  resolution = 300
)


############################################################################
##              READ IN AND VISUALISE GSEA RESULTS
############################################################################

# list all the cell type results
gsea.files <- list.files("/Users/mp34/sanger/team154pc/mp34/lung/NTCU/scRNA", pattern = "GSEA", all.files = T, full.names = T)
gsea.files <- gsea.files[-1]
gsea.data <- lapply(gsea.files, read.table, header = T, sep = ",")

# get the cell types
celltype.ids <- str_split_fixed(basename(gsea.files), "_", 5)[,1]

i <- 1
gsea.enr.list <- list()

# iterate over dataframes and subset to significant genes
for (i in 1:length(celltype.ids)){
  
  # wrangle into shape
  data.tmp <- data.frame(t(gsea.data[[i]]))
  data.tmp$pathways <- rownames(data.tmp)
  data.tmp <- data.tmp[-c(1,2), ]
  colnames(data.tmp) <- c(paste0(celltype.ids[i], "_pvals"), "pathways")
  data.tmp[,paste0(celltype.ids[i], "_pvals")] <- as.numeric(as.character(data.tmp[,paste0(celltype.ids[i], "_pvals")] )) 
  
  # transform p values 
  data.tmp[data.tmp[,paste0(celltype.ids[i], "_pvals")] == 0, paste0(celltype.ids[i], "_pvals")] <- min(data.tmp[data.tmp[,paste0(celltype.ids[i], "_pvals")] != 0, paste0(celltype.ids[i], "_pvals")])

  # save to combine
  gsea.enr.list[[celltype.ids[i]]] <- data.tmp
  
}

plot_data_log <- purrr::reduce(gsea.enr.list, full_join, by = "pathways")
plot_data_log[which(is.na(plot_data_log), arr.ind = T)] <- 1 

colnames(plot_data_log) <- str_split_fixed(colnames(plot_data_log), "_", 2)[,1]
rownames(plot_data_log) <- plot_data_log$pathways
plot_data_log$pathways <- NULL

# set yaxis lab 
rownames(plot_data_log) <- gsub("HALLMARK_", "", rownames(plot_data_log))

log.matrix <- as.matrix(t(-log10(plot_data_log)))
col_fun = colorRamp2(c(min(log.matrix), 25, 50), c("darkgreen","white", "#B63679FF"))
ht1 <- Heatmap(log.matrix, name = "Missense dN/dS", col = col_fun,
               # cell_fun = function(j, i, x, y, width, height, fill) {
               #   if(qmis.matrix[i, j] < 0.05)
               #     grid.text(sprintf("%.2f", sel_cv.matrix[i, j]), x, y, gp = gpar(fontsize = 10, fontface = "bold", col = "black"))
               # },
               cluster_rows = T, 
               cluster_columns = T,
               show_row_names = TRUE, 
               row_names_side = "left",
               show_row_dend = FALSE,
               show_column_names = TRUE,
               column_names_rot = 45,
               row_names_gp = gpar(fontsize = 6, fontface = "bold"), 
               column_names_gp = gpar(fontsize = 6, fontface = "bold"), 
               rect_gp = gpar(col = "white", lwd = 1), border_gp = gpar(col = "black", lwd = 2),
               column_title_gp = gpar(fontsize = 15, fontface = "bold"),
               heatmap_legend_param = list(color_bar = "continous",
                                           at = c(min(log.matrix), 25, 50),
                                           title = "-log10(p)"), 
               border = T)


pdf(paste0("/Users/mp34/team154_campbell/plots/NTCU/GSEA_heatmap_celltypes.pdf"), width=10, height = 3 , pointsize=0.1)
print(ht1)
dev.off()

