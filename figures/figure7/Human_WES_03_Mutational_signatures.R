
################################################################################################################################################
##                                                                                                                      
##  HUMAN WES data DATA - Mutational signature by clonality
##                                                                                                                      
##  Date: 14 August 2024                                                                                                                   
##  
##  Author: Ahmed Alhendi                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

source("scripts/WES_functions.R")

# List of required packages
library(dplyr)
library(data.table)
library(tibble)
library(stringr)
library(MutationalPatterns)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(cowplot)
library(ggsci)

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

