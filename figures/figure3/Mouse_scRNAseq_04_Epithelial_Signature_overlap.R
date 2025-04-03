
################################################################################################################################################
##                                                                                                                      
##  Signature Overlap with Human Epithelial from previous publications
##                                                                                                                      
##  Date: 08 Dec 2025                                                                                                                   
##  
##  Author: Ahmed Alhendi                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################


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
cells.type.order <- c('Basal proliferative','Basal','Basal Tgm2+','Basal Mecom+','Basal Krt14+','Krt4/13+','Secretory',
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

