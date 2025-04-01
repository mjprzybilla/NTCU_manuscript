################################################################################################################################################
##                                                                                                                      
##  Signature Overlap with Human Epithelial from previous publications
##                                                                                                                      
##  Date: 18 July 2024                                                                                                                   
##  
##  Author: Ahmed Alhendi                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

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

