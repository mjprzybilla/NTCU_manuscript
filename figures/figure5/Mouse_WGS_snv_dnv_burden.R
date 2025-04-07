################################################################################################################################################
##                                                                                                                      
##  VISUALIZE SNVS, DBS AND INDELS PER SAMPLE FOR NTCU MICE
##                                                                                                                      
##  Date: 24 MAY 2023                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

# clear workspace beforehand
rm(list = ls())

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("BiocManager", "RColorBrewer", "data.table", "lsa", "lattice", "reshape2",
                      "ggrepel", "readr", "stringr", "tidyverse", "hdp", "sigfit", "BSgenome.Hsapiens.UCSC.hg38",
                      "nrmisc", "cowplot", "ggpubr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))


#####################################################################################
# SET PARAMETERS AND READ IN DATA
#####################################################################################

# define colours for signatures
sig.cols=c("grey80","peachpuff","forestgreen","firebrick","steelblue","pink2", "turquoise1",'orange2',"chartreuse","mediumorchid3","grey20", "yellow2","mediumaquamarine","tomato4")

# output directory
output.dir <- "/Users/mp34/team154_campbell/plots/NTCU"

# READ SNV INFORMATION
snv_counts <- read.table(paste0( "/Users/mp34/sanger/team154pc/mp34/lung/NTCU/WGS/HDP/clones", "/trinuc_mut_mat_full.txt"))
snv.signatures <- read.table(paste0("/Users/mp34/sanger/team154pc/mp34/lung/NTCU/WGS/HDP/clones/sample_signatures_Lung_2023_08_20.txt"))

# READ DBS INFORMATION
dbs_counts <- read.table("/Users/mp34/sanger/team154pc/mp34/lung/NTCU/WGS/data/processed/dbs_mutation_matrix_clones.tsv", sep = "\t")
rownames(dbs_counts) <- dbs_counts$MutationType
dbs_counts$MutationType <- NULL

dbs.signatures <- read.table("/Users/mp34/sanger/team154pc/mp34/lung/NTCU/WGS/data/processed/ALL_clone_FITTED_DBS_Signature_countribution_matrix.txt", sep = "\t")

# wrangle into shape
for (i in 1:ncol(dbs.signatures)){

  dbs.signatures[,i] <- dbs.signatures[,i] / colSums(dbs.signatures)[i]
}

# READ INDEL INFORMATION
indel_counts <- read.delim(paste0("/Users/mp34/sanger/team154pc/mp34/lung/NTCU/WGS/data/processed/ALL_clones_INDEL_count_matrix.txt"), header = T)
indel_counts <- as.data.frame(t(indel_counts))
indel_counts <- as.data.frame(rowSums(indel_counts))
indel_counts$sampleID <- rownames(indel_counts)
colnames(indel_counts) <- c("num_indels", "sampleID")

# # ID INFORMATION
# id.signatures <- fread(paste0("/Users/mp34/sanger/team154pc/mp34/lung/NTCU/WGS/MutationalPatterns/INDELS/SigFit_Indel_ExtractionClones_Exposures.tsv"))
# colnames(id.signatures)[1] <- c("sampleID")

id.signatures <- read.table(paste0("/Users/mp34/sanger/team154pc/mp34/lung/NTCU/WGS/data/processed/ALL_clones_FITTED_INDEL_Signature_countribution_matrix.txt"))

# wrangle into shape
for (i in 1:ncol(id.signatures)){

  id.signatures[,i] <- id.signatures[,i] / colSums(id.signatures)[i]
}

id.signatures <- t(id.signatures)

# READ IN NTCU INFORMATION
metadata <- fread("/Users/mp34/team154_campbell/data/NTCU_complete_metadata.txt", sep="\t", header=T, data.table=F)

# READ IN CLONALITY ESTIMATES
clonality.df <- fread("/Users/mp34/sanger/team154pc/mp34/lung/NTCU/WGS/clonality_estimation/NTCU_clonalityEstimate_table.txt")
clonality.df <- clonality.df[clonality.df$clonal == TRUE & clonality.df$Mutations > 500, ]

# save the sensitivity adjusted burden 
sensit_tip_tbl <- read.csv("/Users/mp34/team154_campbell/data/NTCU_per_clone_sensitivity_adjusted_lobe_tip_complete.csv")

#####################################################################################
# VISUALIZE SNVS FOR BOTH MICE
#####################################################################################

# check mutation numbers
# snv.melt <- data.frame("num_muts" = rowSums(snv_counts))
# snv.melt$sampleID <- rownames(snv.melt)

snv.melt <- sensit_tip_tbl[,c("patient_cluster", "tip_length")]
snv.melt$patient_cluster <- paste0(str_split_fixed(snv.melt$patient_cluster, "_", 4)[,1], "_",
                                   str_split_fixed(snv.melt$patient_cluster, "_", 4)[,2], "_", 
                                   str_split_fixed(snv.melt$patient_cluster, "_", 4)[,3], "_" , 
                                   "Cl.", str_split_fixed(snv.melt$patient_cluster, "_", 4)[,4])
colnames(snv.melt) <- c("sampleID", "num_muts")

# merge with sbs signatures
snv.signatures$sampleID <- rownames(snv.signatures)
snv.sig.melt <- reshape2::melt(snv.signatures)

combined.clone.snv.table <- merge(snv.melt, snv.sig.melt, by = "sampleID")
combined.clone.snv.table$patientID <- substr(combined.clone.snv.table$sampleID, 1, 6)
combined.clone.snv.table$sig_muts <- round(combined.clone.snv.table$num_muts * combined.clone.snv.table$value)

snv.data.of.interest <- combined.clone.snv.table[, c("sampleID", "patientID", "num_muts", "variable", "sig_muts")]
# snv.data.of.interest <- snv.data.of.interest[snv.data.of.interest$sampleID %in% clonality.df$Sample,]
snv.data.of.interest <- unique(snv.data.of.interest)

ordered.data <- snv.data.of.interest
ordered.data[, c("sig_muts", "variable")] <- NULL
ordered.data <- unique(ordered.data) %>% group_by(patientID) %>% arrange(patientID, num_muts) %>% mutate(index = row_number(patientID))

# create dataframe with index to order samples along x axis
snv.data.of.interest <- merge(snv.data.of.interest, ordered.data[, c("sampleID", "index")], by = c("sampleID"))

ordered.sample_snv_counts <- snv.data.of.interest %>% group_by(patientID) %>% mutate(mean_snvs = round(mean(num_muts)))
ordered.sample_snv_counts <- ordered.sample_snv_counts %>% group_by(patientID) %>% mutate(mean_samples = mean(index))
ordered.sample_snv_counts$variable <- factor(ordered.sample_snv_counts$variable, levels = c( "SBS18", "SBS11", "SBS32", "SBS36", "N2", "SBS5"))

# colours for the patients
cols=c(brewer.pal(length(unique(ordered.sample_snv_counts$variable)), "Paired"),"magenta","firebrick")
names(cols) <- c("SBS5", "SBS11", "SBS18", "SBS32", "SBS36", "N2")

ggplot(ordered.sample_snv_counts, aes(x=index, y=sig_muts, fill = variable)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_color_manual(values = "black") +
  scale_fill_manual(values = cols) +
  labs(y = "# SNV", x = "") +
  facet_grid(~ patientID, scales='free_x', space = "free_x") +
  theme_classic() + 
  geom_hline(aes(yintercept = mean_snvs, group = patientID), colour = 'darkgrey', size = 1) +
  geom_text(aes(mean_samples, mean_snvs, label = mean_snvs, hjust = 1.25, vjust = -1), size = 3) +
  scale_x_continuous(expand = c(0.005, 0.005)) +
  scale_y_continuous(breaks = c(1000, 5000, 10000, 25000, 40000), expand = c(0.005, 0.005)) +
  # scale_y_log10(expand = c(0.005, 0.005)) +
  theme(line = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.1),
        strip.text = element_text(face="bold", size=12, colour = "black",margin = margin(t = 5, r = 50, b = 5, l = 50)),
        # strip.text = element_blank(),
        # strip.background = element_rect(fill="white", colour="black", size=1),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 0, angle = 45, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
        axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, vjust = .5, face = "bold", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(color = "black", size = 6, hjust = .5, face = "bold"), 
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold"), 
        legend.position = "bottom")
ggsave(paste0(output.dir, "/NTCU_clones_based_snv_burden.pdf"), width = 6, height = 3, dpi = 600)

#####################################################################################
# VISUALIZE SBS SIGNATURES
#####################################################################################

# wrangle the signatures into shape
snv.signatures$sampleID <- rownames(snv.signatures)
snv.melt <- reshape2::melt(snv.signatures)
colnames(snv.melt) <- c("sampleID", "sbs_sigs", "sbs_contribution")

round(colMeans(snv.signatures[,1:(ncol(snv.signatures)-1)]),2)

# get the information of interest
snv.sigs.data.of.interest <- combined.clone.snv.table[, c("sampleID", "patientID", "variable", "value")]
colnames(snv.sigs.data.of.interest) <- c("sampleID", "patientID", "sbs_sigs", "sbs_contribution")
snv.sigs.data.of.interest <- unique(snv.sigs.data.of.interest)

# 
ordered.snv.melt <- merge(snv.sigs.data.of.interest, ordered.sample_snv_counts[,c("sampleID", "index", "mean_snvs")], by = c("sampleID"))
ordered.snv.melt <- unique(ordered.snv.melt)
ordered.snv.melt <- ordered.snv.melt[order(ordered.snv.melt$index),]

ordered.snv.melt$sbs_sigs <- factor(ordered.snv.melt$sbs_sigs, levels = c( "SBS18", "SBS11", "SBS32", "SBS36", "N2", "SBS5"))

ggplot(ordered.snv.melt, aes(x=index, y=sbs_contribution, fill = sbs_sigs)) + 
  geom_bar(stat="identity") + 
  scale_color_manual(values = "black") +
  scale_fill_manual(values = cols) +
  labs(y = "# SNV", x = "") +
  facet_grid(~ patientID, scales='free_x', space = "free_x") +
  theme_classic() + 
  scale_x_continuous(expand = c(0.005, 0.005)) +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  theme(line = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.1),
        # strip.text = element_text(face="bold", size=12, colour = "black",margin = margin(t = 5, r = 50, b = 5, l = 50)),
        strip.text = element_blank(),
        # strip.background = element_rect(fill="white", colour="black", size=1),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 0, angle = 45, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
        axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, vjust = .5, face = "bold", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(color = "black", size = 6, hjust = .5, face = "bold"), 
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold"), 
        legend.position = "bottom")
ggsave(paste0(output.dir, "/NTCU_sample_sbs_contributions.pdf"), width = 6, height = 3, dpi = 600)

#####################################################################################
# VISUALIZE DBS FOR BOTH MICE
#####################################################################################

# check mutation numbers
dbs.melt <- data.frame("num_muts" = rowSums(t(dbs_counts)))
dbs.melt$sampleID <- rownames(dbs.melt)

# merge with sbs signatures
dbs.signatures <- as.data.frame(t(dbs.signatures))
dbs.signatures$sampleID <- rownames(dbs.signatures)
dbs.sig.melt <- reshape2::melt(dbs.signatures)

combined.clone.dbs.table <- merge(dbs.melt, dbs.sig.melt, by = "sampleID")
combined.clone.dbs.table$patientID <- substr(combined.clone.dbs.table$sampleID, 1, 6)
combined.clone.dbs.table$sig_muts <- round(combined.clone.dbs.table$num_muts * combined.clone.dbs.table$value)

dbs.data.of.interest <- combined.clone.dbs.table[, c("sampleID", "patientID", "num_muts", "variable", "sig_muts")]
dbs.data.of.interest <- unique(dbs.data.of.interest)

#
dbs.data.of.interest <- merge(dbs.data.of.interest, ordered.sample_snv_counts[,c("sampleID", "index")], by = c("sampleID"))
dbs.data.of.interest <- unique(dbs.data.of.interest)
dbs.data.of.interest <- dbs.data.of.interest[order(dbs.data.of.interest$index),]

dbs.data.of.interest <- dbs.data.of.interest %>% group_by(patientID) %>% mutate(mean_snvs = round(mean(num_muts)))
dbs.data.of.interest <- dbs.data.of.interest %>% group_by(patientID) %>% mutate(mean_samples = mean(index))

#
dbs.data.of.interest$variable <- factor(dbs.data.of.interest$variable, levels = c(paste0("DBS", c(1:11))))


# colours for the patients
nb.cols <- length(unique(dbs.data.of.interest$variable))
mycolors <- colorRampPalette(ggsci::pal_lancet(palette = "lanonc")(9))(nb.cols)
names(mycolors) <- c(paste0("DBS", c(4, 5, 11)))

ggplot(dbs.data.of.interest, aes(x=index, y=sig_muts, fill = variable)) +
  geom_bar(position="stack", stat="identity") +
  scale_color_manual(values = "black") +
  scale_fill_manual(values = mycolors) +
  labs(y = "# DNVs", x = "") +
  facet_grid(~ patientID, scales='free_x', space = "free_x") +
  theme_classic() +
  geom_hline(aes(yintercept = mean_snvs, group = patientID), colour = 'darkgrey', size = 1) +
  geom_text(aes(mean_samples, mean_snvs, label = mean_snvs, hjust = 1.25, vjust = -1), size = 3) +
  scale_x_continuous(expand = c(0.005, 0.005)) +
  scale_y_continuous(breaks = c(0, 25, 50, 75), expand = c(0.005, 0.005)) +
  # scale_y_log10(expand = c(0.005, 0.005)) +
  theme(line = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.1),
        strip.text = element_text(face="bold", size=12, colour = "black",margin = margin(t = 5, r = 50, b = 5, l = 50)),
        # strip.text = element_blank(),
        # strip.background = element_rect(fill="white", colour="black", size=1),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 0, angle = 45, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 1, vjust = 0, face = "bold"),
        axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, vjust = .5, face = "bold",
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(color = "black", size = 6, hjust = .5, face = "bold"),
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold"),
        legend.position = "bottom")
ggsave(paste0(output.dir, "/NTCU_clones_based_dbs_burden.pdf"), width = 6, height = 3, dpi = 600)

#####################################################################################
# VISUALIZE DBS SIGNATURES
#####################################################################################

# wrangle the signatures into shape
colnames(dbs.sig.melt) <- c("sampleID", "sbs_sigs", "sbs_contribution")

# MEAN CONTRIBUTION PER SIGNATURE
round(colMeans(dbs.signatures[,1:(ncol(dbs.signatures)-1)]),2)
# SBS5 SBS40  SBS4 SBS12 SBS92  SBS2 SBS13
# 0.42  0.35  0.17  0.02  0.03  0.00  0.01

# get the information of interest
dbs.sigs.data.of.interest <- combined.clone.dbs.table[, c("sampleID", "patientID", "variable", "value")]
colnames(dbs.sigs.data.of.interest) <- c("sampleID", "patientID", "sbs_sigs", "sbs_contribution")
dbs.sigs.data.of.interest <- unique(dbs.sigs.data.of.interest)

#
ordered.dbs.melt <- merge(dbs.sigs.data.of.interest, ordered.sample_snv_counts[,c("sampleID", "index", "mean_snvs")], by = c("sampleID"))
ordered.dbs.melt <- unique(ordered.dbs.melt)
ordered.dbs.melt <- ordered.dbs.melt[order(ordered.dbs.melt$index),]

#
ordered.dbs.melt$sbs_sigs <- factor(ordered.dbs.melt$sbs_sigs, levels = c(paste0("DBS", c(1:11))))

ggplot(ordered.dbs.melt, aes(x=index, y=sbs_contribution, fill = sbs_sigs)) +
  geom_bar(stat="identity") +
  scale_color_manual(values = "black") +
  scale_fill_manual(values = mycolors) +
  labs(y = "# SNV", x = "") +
  facet_grid(~ patientID, scales='free_x', space = "free_x") +
  theme_classic() +
  scale_x_continuous(expand = c(0.005, 0.005)) +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  theme(line = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.1),
        # strip.text = element_text(face="bold", size=12, colour = "black",margin = margin(t = 5, r = 50, b = 5, l = 50)),
        strip.text = element_blank(),
        # strip.background = element_rect(fill="white", colour="black", size=1),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 0, angle = 45, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 1, vjust = 0, face = "bold"),
        axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, vjust = .5, face = "bold",
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(color = "black", size = 6, hjust = .5, face = "bold"),
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold"),
        legend.position = "bottom")
ggsave(paste0(output.dir, "/NTCU_clones_dbs_contributions.pdf"), width = 6, height = 3, dpi = 600)

#####################################################################################
# VISUALIZE INDELS FOR BOTH MICE
#####################################################################################

# merge with sbs signatures

id.signatures <- as.data.frame(id.signatures)
id.signatures$V1 <- rownames(id.signatures)
# rownames(id.signatures) <- id.signatures$V1
id.sig.melt <- reshape2::melt(id.signatures)
colnames(id.sig.melt)[1] <- "sampleID"

combined.clone.indel.table <- merge(indel_counts, id.sig.melt, by = "sampleID")
combined.clone.indel.table$patientID <- substr(combined.clone.indel.table$sampleID, 1, 6)
combined.clone.indel.table$indel_muts <- round(combined.clone.indel.table$num_indels * combined.clone.indel.table$value)

indel.data.of.interest <- combined.clone.indel.table[, c("sampleID", "patientID", "num_indels", "variable", "indel_muts")]
indel.data.of.interest <- unique(indel.data.of.interest)

# 
indel.data.of.interest <- merge(indel.data.of.interest, ordered.sample_snv_counts[,c("sampleID", "index")], by = c("sampleID"))
indel.data.of.interest <- unique(indel.data.of.interest)
indel.data.of.interest <- indel.data.of.interest[order(indel.data.of.interest$index),]

indel.data.of.interest <- indel.data.of.interest %>% group_by(patientID) %>% mutate(mean_indels = round(mean(num_indels)))
indel.data.of.interest <- indel.data.of.interest %>% group_by(patientID) %>% mutate(mean_samples = mean(index))

indel.data.of.interest$variable <- factor(indel.data.of.interest$variable, levels = c(paste0("ID", c(1:11))))

# colours for the patients
nb.cols <- length(unique(indel.data.of.interest$variable))
mycolors <- colorRampPalette(ggsci::pal_lancet(palette = "lanonc")(9))(nb.cols)
names(mycolors) <- c(paste0("ID", c(1, 2, 5, 8, 9)))

ggplot(indel.data.of.interest, aes(x=index, y=indel_muts, fill = variable)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_color_manual(values = "black") +
  scale_fill_manual(values = mycolors) +
  labs(y = "# DNVs", x = "") +
  facet_grid(~ patientID, scales='free_x', space = "free_x") +
  theme_classic() + 
  geom_hline(aes(yintercept = mean_indels, group = patientID), colour = 'darkgrey', size = 1) +
  geom_text(aes(mean_samples, mean_indels, label = mean_indels, hjust = 1.25, vjust = -1), size = 3) +
  scale_x_continuous(expand = c(0.005, 0.005)) +
  scale_y_continuous(breaks = c(0, 5, 10, 25, 50), expand = c(0.005, 0.005)) +
  # scale_y_log10(expand = c(0.005, 0.005)) +
  theme(line = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.1),
        strip.text = element_text(face="bold", size=12, colour = "black",margin = margin(t = 5, r = 50, b = 5, l = 50)),
        # strip.text = element_blank(),
        # strip.background = element_rect(fill="white", colour="black", size=1),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 0, angle = 45, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
        axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, vjust = .5, face = "bold", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(color = "black", size = 6, hjust = .5, face = "bold"), 
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold"), 
        legend.position = "bottom")
ggsave(paste0(output.dir, "/NTCU_sample_based_indel_burden.pdf"), width = 6, height = 3, dpi = 600)

#####################################################################################
# VISUALIZE DBS SIGNATURES
#####################################################################################

# wrangle the signatures into shape
colnames(id.sig.melt) <- c("sampleID", "indel_sigs", "indel_contribution")

# get the information of interest
indel.sigs.data.of.interest <- combined.clone.indel.table[, c("sampleID", "patientID", "variable", "value")]
colnames(indel.sigs.data.of.interest) <- c("sampleID", "patientID", "indel_sigs", "indel_contribution")
indel.sigs.data.of.interest <- unique(indel.sigs.data.of.interest)

# 
ordered.indel.melt <- merge(indel.sigs.data.of.interest, ordered.sample_snv_counts[,c("sampleID", "index")], by = c("sampleID"))
ordered.indel.melt <- unique(ordered.indel.melt)
ordered.indel.melt <- ordered.indel.melt[order(ordered.indel.melt$index),]

ggplot(ordered.indel.melt, aes(x=index, y=indel_contribution, fill = indel_sigs)) + 
  geom_bar(stat="identity") + 
  scale_color_manual(values = "black") +
  scale_fill_manual(values = mycolors) +
  labs(y = "# SNV", x = "") +
  facet_grid(~ patientID, scales='free_x', space = "free_x") +
  theme_classic() + 
  scale_x_continuous(expand = c(0.005, 0.005)) +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  theme(line = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.1),
        # strip.text = element_text(face="bold", size=12, colour = "black",margin = margin(t = 5, r = 50, b = 5, l = 50)),
        strip.text = element_blank(),
        # strip.background = element_rect(fill="white", colour="black", size=1),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 0, angle = 45, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
        axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, vjust = .5, face = "bold", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(color = "black", size = 6, hjust = .5, face = "bold"), 
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold"), 
        legend.position = "bottom")
ggsave(paste0(output.dir, "/NTCU_sample_indel_contributions.pdf"), width = 6, height = 3, dpi = 600)


#####################################################################################
# CREATE COMBINED PLOT
#####################################################################################

# colours for the patients
cols=c(brewer.pal(length(unique(ordered.sample_snv_counts$variable)), "Paired"),"magenta","firebrick")
names(cols) <- c("SBS5", "SBS11", "SBS18", "SBS32", "SBS36", "N2")

p1 <- ggplot(ordered.sample_snv_counts, aes(x=index, y=sig_muts, fill = variable)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_color_manual(values = "black") +
  scale_fill_manual(values = cols) +
  labs(y = "# SNV", x = "") +
  facet_grid(~ patientID, scales='free_x', space = "free_x") +
  theme_classic() + 
  geom_hline(aes(yintercept = mean_snvs, group = patientID), colour = 'darkgrey', size = 1) +
  geom_text(aes(mean_samples, mean_snvs, label = mean_snvs, hjust = 1.25, vjust = -1), size = 3) +
  scale_x_continuous(expand = c(0.005, 0.005)) +
  scale_y_continuous(breaks = c(1000, 5000, 10000, 25000, 40000), expand = c(0.005, 0.005)) +
  # scale_y_log10(expand = c(0.005, 0.005)) +
  theme(line = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.1),
        # strip.text = element_text(face="bold", size=12, colour = "black",margin = margin(t = 5, r = 50, b = 5, l = 50)),
        strip.text = element_blank(),
        # strip.background = element_rect(fill="white", colour="black", size=1),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(),
        axis.text.x = element_text(color = "black", size = 0, angle = 45, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
        axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, vjust = .5, face = "bold", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(color = "black", size = 6, hjust = .5, face = "bold"), 
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold"), 
        legend.position = "right",
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

p2 <- ggplot(ordered.snv.melt, aes(x=index, y=sbs_contribution, fill = sbs_sigs)) + 
  geom_bar(stat="identity") + 
  scale_color_manual(values = "black") +
  scale_fill_manual(values = cols) +
  labs(y = "SBS Contribution", x = "") +
  facet_grid(~ patientID, scales='free_x', space = "free_x") +
  theme_classic() + 
  scale_x_continuous(expand = c(0.005, 0.005)) +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  theme(line = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.1),
        # strip.text = element_text(face="bold", size=12, colour = "black",margin = margin(t = 5, r = 50, b = 5, l = 50)),
        strip.text = element_blank(),
        # strip.background = element_rect(fill="white", colour="black", size=1),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(),
        axis.text.x = element_text(color = "black", size = 0, angle = 45, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
        axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, vjust = .5, face = "bold", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(color = "black", size = 6, hjust = .5, face = "bold"), 
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold"), 
        legend.position = "right",
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

# colours for the patients
nb.cols <- length(unique(dbs.data.of.interest$variable))
mycolors <- colorRampPalette(ggsci::pal_lancet(palette = "lanonc")(9))(nb.cols)
names(mycolors) <- c(paste0("DBS", c(4, 5, 11)))

p3 <- ggplot(dbs.data.of.interest, aes(x=index, y=sig_muts, fill = variable)) +
  geom_bar(position="stack", stat="identity") +
  scale_color_manual(values = "black") +
  scale_fill_manual(values = mycolors) +
  labs(y = "# DNVs", x = "") +
  facet_grid(~ patientID, scales='free_x', space = "free_x") +
  theme_classic() +
  geom_hline(aes(yintercept = mean_snvs, group = patientID), colour = 'darkgrey', size = 1) +
  geom_text(aes(mean_samples, mean_snvs, label = mean_snvs, hjust = 1.25, vjust = -1), size = 3) +
  scale_x_continuous(expand = c(0.005, 0.005)) +
  scale_y_continuous(breaks = c(0, 25, 50, 75), expand = c(0.005, 0.005)) +
  # scale_y_log10(expand = c(0.005, 0.005)) +
  theme(line = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.1),
        # strip.text = element_text(face="bold", size=12, colour = "black",margin = margin(t = 5, r = 50, b = 5, l = 50)),
        strip.text = element_blank(),
        # strip.background = element_rect(fill="white", colour="black", size=1),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(),
        axis.text.x = element_text(color = "black", size = 0, angle = 45, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 1, vjust = 0, face = "bold"),
        axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, vjust = .5, face = "bold",
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(color = "black", size = 6, hjust = .5, face = "bold"),
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold"),
        legend.position = "right",
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

p4 <- ggplot(ordered.dbs.melt, aes(x=index, y=sbs_contribution, fill = sbs_sigs)) +
  geom_bar(stat="identity") +
  scale_color_manual(values = "black") +
  scale_fill_manual(values = mycolors) +
  labs(y = "DBS Contribution", x = "") +
  facet_grid(~ patientID, scales='free_x', space = "free_x") +
  theme_classic() +
  scale_x_continuous(expand = c(0.005, 0.005)) +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  theme(line = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.1),
        # strip.text = element_text(face="bold", size=12, colour = "black",margin = margin(t = 5, r = 50, b = 5, l = 50)),
        strip.text = element_blank(),
        # strip.background = element_rect(fill="white", colour="black", size=1),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(),
        axis.text.x = element_text(color = "black", size = 0, angle = 45, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 1, vjust = 0, face = "bold"),
        axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, vjust = .5, face = "bold",
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(color = "black", size = 6, hjust = .5, face = "bold"),
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold"),
        legend.position = "right",
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

# colours for the patients
nb.cols <- length(unique(indel.data.of.interest$variable))
mycolors <- colorRampPalette(ggsci::pal_lancet(palette = "lanonc")(9))(nb.cols)

p5 <- ggplot(indel.data.of.interest, aes(x=index, y=indel_muts, fill = variable)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_color_manual(values = "black") +
  scale_fill_manual(values = mycolors) +
  labs(y = "# Indels", x = "") +
  facet_grid(~ patientID, scales='free_x', space = "free_x") +
  theme_classic() + 
  geom_hline(aes(yintercept = mean_indels, group = patientID), colour = 'darkgrey', size = 1) +
  geom_text(aes(mean_samples, mean_indels, label = mean_indels, hjust = 1.25, vjust = -1), size = 3) +
  scale_x_continuous(expand = c(0.005, 0.005)) +
  scale_y_continuous(breaks = c(0, 5, 10, 25, 50), expand = c(0.005, 0.005)) +
  # scale_y_log10(expand = c(0.005, 0.005)) +
  theme(line = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.1),
        # strip.text = element_text(face="bold", size=12, colour = "black",margin = margin(t = 5, r = 50, b = 5, l = 50)),
        strip.text = element_blank(),
        # strip.background = element_rect(fill="white", colour="black", size=1),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(),
        axis.text.x = element_text(color = "black", size = 0, angle = 45, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
        axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, vjust = .5, face = "bold", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(color = "black", size = 6, hjust = .5, face = "bold"), 
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold"), 
        legend.position = "right", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

p6 <- ggplot(ordered.indel.melt, aes(x=index, y=indel_contribution, fill = indel_sigs)) + 
  geom_bar(stat="identity") + 
  scale_color_manual(values = "black") +
  scale_fill_manual(values = mycolors) +
  labs(y = "ID Contribution", x = "") +
  facet_grid(~ patientID, scales='free_x', space = "free_x") +
  theme_classic() + 
  scale_x_continuous(expand = c(0.005, 0.005)) +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  theme(line = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.1),
        # strip.text = element_text(face="bold", size=12, colour = "black",margin = margin(t = 5, r = 50, b = 5, l = 50)),
        strip.text = element_blank(),
        # strip.background = element_rect(fill="white", colour="black", size=1),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(),
        axis.text.x = element_text(color = "black", size = 0, angle = 45, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
        axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, vjust = .5, face = "bold", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(color = "black", size = 6, hjust = .5, face = "bold"), 
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold"), 
        legend.position = "right",
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

# manually setting the number of rows, auto-generate upper-case labels
pdf(paste0(output.dir, "/NTCU_snv_dbs_indel_visualisation_all.pdf"), width = 7, height = 7)
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow = 6, align = "v")
# cowplot::plot_grid(p1, p2, p5, p6, nrow = 4, align = "v")
dev.off()

