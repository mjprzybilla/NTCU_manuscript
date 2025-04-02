
################################################################################################################################################
##                                                                                                                      
##  HUMAN WES data DATA - Generate dn/ds gene-level selection
##                                                                                                                      
##  Date: 14 August 2024                                                                                                                   
##  
##  Author: Ahmed Alhendi                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

# Load custom plotting functions
source("scripts/WES_functions.R")

# List of required packages

library(dndscv)
library(dplyr)
library(tibble)
library(data.table)
library(stringr)
library(seqinr)
library(Biostrings)
library(MASS)
library(ggrepel)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(ggsci)
library(GenomicRanges)


# Load input data
driver.genes <- fread("data/WES/natMed.CIS.tracerx.LUSC.driverlist.csv")
df.somatic_variants <- fread("data/WES/somatic_variants.hg19_multianno.clonality.drivers.tsv")

# Process for dnds cv analysis
#trunk.df <- run_dndscv_analysis(df.somatic_variants, driver.genes$Gene.Symbol, "Truncal")
#sub.df <- run_dndscv_analysis(df.somatic_variants, driver.genes$Gene.Symbol, "Subclonal")

trunk.df <- fread("data/WES/Truncal.DNDS.samples.txt")
sub.df <- fread("data/WES/Subclonal.DNDS.samples.txt")

gene.summary <- fread("data/WES/all_summary_geneSummary.txt")

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
