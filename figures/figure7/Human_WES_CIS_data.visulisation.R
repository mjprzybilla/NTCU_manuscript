################################################################################################################################################
##                                                                                                                      
##  HUMAN WES data DATA VISULIZATION - FIGURE 7 & SUPP FIGURE 15 & SUPP FIGURE 13
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
  "dplyr", "data.table", "tibble", "stringr",
  
  # Genomic and Mutational Analysis
  "MutationalPatterns", "BSgenome", "BSgenome.Hsapiens.UCSC.hg19", "CONIPHER", "maftools"
  
  # Visualization
  "ggrepel","gridExtra", "circlize", "ggplot2", "cowplot", "ggsci", "ComplexHeatmap"
)

# Load all required packages
suppressMessages({
  lapply(required_packages, library, character.only = TRUE)
})


# Load custom plotting functions
source("scripts/WES_functions.R")
source("scripts/cloneMap.R")


############################################################################
##           Visualize of CONIPHER trees
############################################################################
# Define working directory
workDir <- paste0("data/WES/")
allseedingTable <- fread("data/WES/conipher_seedingTable.tsv")

# Function to process each patient
process_patient_tree <- function(patientid, data) {
    seedingTable <- allseedingTable %>% filter(patient_id == patientid)
    fullTreeOutput <- readRDS(file.path(paste0(workDir, patientid, ".tree.RDS")))
    fullTreeOutput$tumour_id <- patientid
    outDir <- file.path(workDir, patientid, "plots")
    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
    
    cols <- pal_cosmic()(4)
    clonalities <- fullTreeOutput$clonality_out$clonality_table_corrected
    cloneColorDF <- data.frame(
        tumour_id = patientid,
        clones = as.character(sort(unique(as.numeric(fullTreeOutput$graph_pyclone$Corrected_tree)))),
        colors = "gray",
        stroke = 1
    )
    
    timepoints <- colnames(clonalities)
    
    # Process each timepoint using a for loop
    for (timepoint in timepoints) {
        tmpdf <- data.frame(sample = timepoints, stringsAsFactors = FALSE) %>%
            mutate(regionOfOrigin = gsub("\\.P[0-9]\\.A[0-9]", "", sample)) %>%
            mutate(regionOfOrigin = sapply(strsplit(regionOfOrigin, "_"), function(x) x[3]))
        
        timepoint.to.use <- tmpdf %>% filter(sample == timepoint) %>% pull(sample)
        invasive.timepoint <- grep("INV", tmpdf$sample, value = TRUE)
        preinvasive.timepoint <- grep("INV", tmpdf$sample, value = TRUE, invert = TRUE)
        
        INVClonality <- get.tumLevel.clonality(clonalities, invasive.timepoint)
        PREClonality <- get.tumLevel.clonality(clonalities, preinvasive.timepoint)
        timepointClonality <- get.tumLevel.clonality(clonalities)
        
        INVPresence <- names(which(INVClonality != "absent"))
        PREPresence <- names(which(PREClonality != "absent"))
        timepointClonal <- names(which(timepointClonality != "absent" & timepointClonality == "clonal"))
        
        if (length(invasive.timepoint) >= 1) {
            sharedClusters <- intersect(INVPresence, PREPresence)
            invasiveUniq <- setdiff(INVPresence, PREPresence)
            preinvasiveUniq <- setdiff(PREPresence, INVPresence)
            
            cloneColorDF$colors[cloneColorDF$clones %in% sharedClusters] <- cols[2]
            cloneColorDF$colors[cloneColorDF$clones %in% invasiveUniq] <- cols[4]
            cloneColorDF$colors[cloneColorDF$clones %in% preinvasiveUniq] <- cols[3]
            cloneColorDF$colors[cloneColorDF$clones %in% timepointClonal] <- cols[1]
        } else {
            print("NO INV")
            cloneColorDF$colors[cloneColorDF$clones %in% PREPresence] <- cols[3]
            cloneColorDF$colors[cloneColorDF$clones %in% timepointClonal] <- cols[1]
        }
        
        timepointTreeOutput <- fullTreeOutput
        timepointTreeOutput$sample_id <- timepoint.to.use
        
        print(paste0("Processing Timepoint ", timepoint.to.use))
        pdf(file.path(outDir, paste0("tree.", timepoint.to.use, ".pdf")), useDingbats = FALSE)
        print(plottingTimepointTreeNodeHighlight(timepointTreeOutput, cloneColorDF, colorBy = "timepoints", timepoints.to.use = timepoint.to.use) + scale_y_reverse())
        dev.off()
    }
    
    # Generate clone color data
    cloneColorDF <- get.seeding.colors(fullTreeOutput, seedingTable)
    
    # Process each timepoint again for seeding colors using a for loop
    for (timepoint in timepoints) {
        timepointTreeOutput <- fullTreeOutput
        timepointTreeOutput$sample_id <- timepoint
        
        print(paste0("Processing Timepoint ", timepoint))
        pdf(file.path(outDir, paste0("tree.", timepoint, ".pdf")), useDingbats = FALSE)
        print(plottingTimepointTreeNodeHighlight(timepointTreeOutput, cloneColorDF, colorBy = "timepoints", timepoints.to.use = timepoint) + scale_y_reverse())
        dev.off()
    }
    
    # Generate and save overall patient tree plot
    pdf(file.path(outDir, paste0("tree.", patientid, ".pdf")), useDingbats = FALSE)
    print(plottingTimepointTreeNodeHighlight(fullTreeOutput, cloneColorDF, colorBy = "timepointAll") + scale_y_reverse())
    dev.off()
    
    # Generate clone maps for the patient
    plottingCloneMaps(fullTreeOutput, cloneColorDF, outDir = outDir)
}

# For loop to apply the function to all patients
for (patientid in unique(allseedingTable$patient_id)) {
    process_patient_tree(patientid, allseedingTable)
}


############################################################################
##           Mutational signature by clonality
############################################################################

# Define file paths
VCF.path <- "data/VCF_multiSite_split_by_clonality/"
SampleInfo.file <- "sample.names.tsv"

# Load sample info
SampleInfo <- fread(SampleInfo.file)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)


## library required

SampleInfo <- fread(SampleInfo.file)

# Process both Truncal and Subclonal mutations
#truncal_sigs <- process_vcf("clonal")
#subclonal_sigs <- process_vcf("subclonal")

truncal_sigs <- fread("data/WES/Top5.clonal.signatures.tsv")
subclonal_sigs <- process_vcf("data/WES/Top5.subclonal.signatures.tsv")


# Combine datasets for summary plot
truncal_sigs$clonality <- "Truncal"
subclonal_sigs$clonality <- "Subclonal"
all_sigs <- rbind(truncal_sigs, subclonal_sigs)

# Load sample info for summary plots
SampleInfo <- fread("sample.names.truncal_subclonal_counts.tsv")
SampleInfo$SampleID <- SampleInfo$SAMPLE
all_sigs <- left_join(all_sigs, SampleInfo, by="SampleID")

# Factor ordering
sigs.to.analyse <- c('SBS4','SBS92','SBS13','SBS2','SBS5','SBS1')
all_sigs$SampleID <- factor(all_sigs$SampleID, levels=unique(SampleInfo$SampleID))
all_sigs$PID <- factor(all_sigs$PID, levels=unique(SampleInfo$PID))
all_sigs$Site <- factor(all_sigs$Site, levels=unique(SampleInfo$Site))
all_sigs$signature <- factor(all_sigs$signature, levels=sigs.to.analyse)
all_sigs$clonality <- factor(all_sigs$clonality, levels=c("Truncal", "Subclonal"))

# Boxplot for clonality contribution
plot_clonality <- ggplot(all_sigs, aes(x=clonality, y=Sig_Contribution)) + 
  geom_boxplot(color="black") +
  geom_jitter(width=0.2, aes(color=signature)) +
  scale_colour_brewer(palette = "Paired") +
  theme_cowplot(font_size = 20) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none") +
  labs(x="", y="Signature Contribution") +
  facet_grid(~ signature, scales="free", space="free") +
  stat_compare_means(label.y = 310, aes(label = after_stat(p.format)))

ggsave(plot_clonality, file="clonality.absolute.contribution.pdf", width=7, height=5)

# Aggregate contributions
all_sigs_agg <- all_sigs %>%
  group_by(PID, Smoking_status, SampleID, signature) %>%
  summarise(Sig_Contribution = sum(Sig_Contribution)) %>%
  ungroup() %>%
  group_by(SampleID) %>%
  mutate(Sig_Contribution_sum = sum(Sig_Contribution),
         Relative_Contribution = Sig_Contribution / Sig_Contribution_sum) %>%
  ungroup() %>%
  mutate(clonality="All")

p1 <- all_sigs %>% ggplot(aes(x=SampleID, y=Sig_Contribution, fill=signature))+
geom_bar(stat="identity", color="black")+
theme_cowplot(font_size=22)+
facet_grid(clonality~PID, space = "free", scales = "free")+
theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
     legend.position="top")+
labs(x="", y="Absolute Contribution")+
scale_fill_brewer(palette = "Paired")+
 theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


p2 <- all_sigs_agg %>% ggplot(aes(x=SampleID, y=Relative_Contribution, fill=signature))+
geom_bar(stat="identity", color="black")+
theme_cowplot(font_size=22)+
facet_grid(clonality~PID, space = "free", scales = "free")+
theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
     legend.position="none")+
labs(x="", y="Relative Contribution")+
scale_fill_brewer(palette = "Paired")


p3 <- all_sigs_agg %>% ggplot(aes(x=SampleID, y=Sig_Contribution, fill=signature))+
geom_bar(stat="identity", color="black")+
theme_cowplot(font_size=22)+
facet_grid(Smoking_status~PID, space = "free", scales = "free")+
theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
     legend.position="none")+
labs(x="", y="Absolute Contribution")+
scale_fill_brewer(palette = "Paired")+
 theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

aggregated_plot <- plot_grid(p1, p3, p2, align='vh', vjust=1, hjust=-10, scale = 1, ncol=1)

ggsave(aggregated_plot, file="Mutational.signatures.plot.pdf", width=12, height=16)


############################################################################
##           Generate SUMMARY PLOT 
############################################################################

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


############################################################################
##           Generate dn/ds gene-level selection plot 
############################################################################

# Load input data
trunk.df <- fread("Truncal.DNDS.samples.txt")
sub.df <- fread("Subclonal.DNDS.samples.txt")
driver.genes <- fread("natMed.CIS.tracerx.LUSC.driverlist.csv")
gene.summary <- fread("all_summary_geneSummary.txt")

total_samples_n <- 12
gene.summary <- gene.summary %>% select(gene_name, MutatedSamples)

# Prepare truncal data
trunk.df <- trunk.df[order(trunk.df$qallsubs_cv), ]
trunk.df <- trunk.df %>% filter(gene_name %in% driver.genes$Gene.Symbol)
trunk.df <- left_join(trunk.df, gene.summary, by="gene_name")

trunk.df.sel_cv <- trunk.df %>% filter(wall_cv > 1) %>%
  mutate(Sites = MutatedSamples / total_samples_n * 100) %>%
  mutate(clonality = "Truncal") %>%
  mutate(Selection = "Truncal favoured")

# Prepare subclonal data
sub.df <- sub.df[order(sub.df$qallsubs_cv), ]
sub.df <- sub.df %>% filter(gene_name %in% driver.genes$Gene.Symbol)
sub.df <- left_join(sub.df, gene.summary, by="gene_name")

sub.df.sel_cv <- sub.df %>% filter(wall_cv > 1) %>%
  mutate(Sites = MutatedSamples / total_samples_n * 100) %>%
  mutate(clonality = "Subclonal") %>%
  mutate(Selection = "Subclonal favoured")

# Find shared truncal and subclonal genes
shared.tunck.sub.genes <- intersect(trunk.df.sel_cv$gene_name, sub.df.sel_cv$gene_name)

# Update Selection for shared genes
trunk.df.sel_cv$Selection[trunk.df.sel_cv$gene_name %in% shared.tunck.sub.genes] <- "Subclonal and truncal"
sub.df.sel_cv$Selection[sub.df.sel_cv$gene_name %in% shared.tunck.sub.genes] <- "Subclonal and truncal"

# Extract subclonal values for shared genes
shared.tunck.sub.genes.subclonal <- sub.df.sel_cv$wall_cv[sub.df.sel_cv$clonality == "Subclonal" &
                                                           sub.df.sel_cv$gene_name %in% shared.tunck.sub.genes]

# Remove shared genes from subclonal data
sub.df.sel_cv <- sub.df.sel_cv %>% filter(!gene_name %in% shared.tunck.sub.genes)

# Combine truncal and subclonal data
df.sel_cv1 <- rbind(trunk.df.sel_cv, sub.df.sel_cv, fill = TRUE)

# Spread data for plotting
df.sel_cv <- df.sel_cv1 %>% tidyr::spread(clonality, wall_cv, fill = 0)
df.sel_cv$Subclonal[df.sel_cv$gene_name %in% shared.tunck.sub.genes] <- shared.tunck.sub.genes.subclonal

# Update factor levels for Selection
df.sel_cv <- df.sel_cv %>% mutate(Selection = factor(Selection, 
                        levels = c("Truncal favoured", "Subclonal and truncal", "Subclonal favoured"))) %>%
  mutate(stroke = ifelse(qallsubs_cv < 0.05, 2, 0.5))

# Add qvalue category
df.sel_cv$qvalue_category <- ifelse(df.sel_cv$qallsubs_cv < 0.05, "qvalue < 0.05", "qvalue >= 0.05")

# Define plot colours
colors <- c("#377EB8", "#984EA3", "#9E2536")

# plotting function for the combined data
plot_dn_ds(df.sel_cv, colors)


############################################################################
##            Generate SAMPLE-SAMPLE CORRELATION PLOT  
############################################################################

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
pdf("Multi_region.sampleSwaps.Pearson.corr.v2.pdf", width=8, height=8)
draw(ht)
dev.off()
