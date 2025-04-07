################################################################################################################################################
##                                                                                                                      
##  CREATE PHYLOGENETIC TREES WITH BRANCHES COLOURED ACCORDING TO SIGNATURES
##                                                                                                                      
##  Date: 17 AUGUST 2021                                                                                                                    
##  
##  Author: Tim Coorens adapted by Moritz Przybilla                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

# clear workspace beforehand
rm(list = ls())

options(stringsAsFactors = F)

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("ape", "reshape2", "ggrepel", "readr", "dplyr", "ggtree", "stringr", "RColorBrewer", "data.table",
                      "hdp", "tidytree", "tidyverse", "ggpubr", "sigfit", "ggplot2", "readxl",
                      "fitdistrplus", "ape", "phytools", "MCMCglmm", "adephylo", "phangorn", "ggridges", "viridis")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

#####################################################################################
# READ IN DATA
#####################################################################################
output.dir <- "/Users/mp34/sanger/team154pc/mp34/lung/NTCU/WGS/HDP/clones"

# read in the files from HDP
mut_example_multi= readRDS("/Users/mp34/sanger/team154pc/mp34/lung/NTCU/WGS/HDP/clones/HDP_multi_chain.Rdata")
mutations=read.table("/Users/mp34/sanger/team154pc/mp34/lung/NTCU/WGS/HDP/clones/trinuc_mut_mat_full.txt") # matrix with number of substitutions (all possible 96 as columns x clusters as rows)
key_table=read.table("/Users/mp34/sanger/team154pc/mp34/lung/NTCU/WGS/HDP/clones/key_table_full.txt", header = T) # table with two columns - Sample (the cluster) and Patient

# remove clones with less than 100 mutations
freq=table(key_table$Patient)

# determine distance 
dp_distn <- comp_dp_distn(mut_example_multi)

# number of samples and signatures
ndp <- nrow(dp_distn$mean)
ncomp <- ncol(dp_distn$mean)

# matrix with signatures as rows and clusters for each sample as columns
# exposures <- t(dp_distn$mean[length(freq)+1+1:nrow(mutations),,drop=FALSE])
# colnames(exposures)=rownames(mutations)

# read in signatures
sample_signatures <- read.table("/Users/mp34/sanger/team154pc/mp34/lung/NTCU/WGS/HDP/clones/sample_signatures_Lung_2023_08_20.txt")
sample_signatures <- t(sample_signatures) 
sample_signatures <- sample_signatures[,rownames(mutations)]
# sample_signatures <- sample_signatures[c("SBS1", "SBS2", "SBS4", "SBS5", "SBS9", "SBS13", "SBS16", "SBS18", "SBS40", "SBS92"),]

# get the patients of interest
patients <- unique(str_split_fixed(key_table$Sample, "_Cl", 2)[,1])

#  signatures of interest
# sigs=rownames(exposures)
sigs=rownames(sample_signatures)
sig_profiles=mut_example_multi@comp_categ_distn$mean
# title.sigs <- c("N0", "SBS5", "N1", "SBS4", "SBS40", "N2", "SBS39", "SBS2/SBS13", "SBS9")
title.sigs <- sigs
title.sigs <- factor(title.sigs, levels = c(paste0("SBS", c(1:96)), paste0("N", c(1:6))))
title.sigs <- title.sigs[order(title.sigs)]
sigs <- as.character(title.sigs)

# get the phylo or tree files
tree.files <- list.files("/Users/mp34/sanger/team154pc/mp34/lung/NTCU/WGS/NDP", pattern = "cluster_tree.phylo", all.files = T, full.names = T, recursive = T)
tree.files <- tree.files[grep("lobe", tree.files)]
tree.files <- tree.files[-grep("snvs_indels", tree.files)]
patient <- patients[1]

for (patient in patients){
  
  print(patient)
  tree = read.tree(tree.files[grep(patient, tree.files)][1])
  tree_df=fortify(tree)
  cols=c(brewer.pal(length(sigs), "Paired"),"magenta","firebrick")
  
  # cols=c("grey80","peachpuff","forestgreen","firebrick","steelblue","pink2",
  #       "turquoise1",'orange2',"chartreuse") # "mediumorchid3","grey20",
  #       # "yellow2","mediumaquamarine","tomato4")
  names(cols)=sigs
  samples=colnames(sample_signatures)[grepl(patient,colnames(sample_signatures))]
  samples=samples[grep("Cl.", samples)]
  branches= str_split_fixed(samples, "_Cl.", 2)[,2]
  branches <- branches[as.character(branches) %in% as.character(tree_df$label)]
  
  pdf(paste0("/Users/mp34/team154_campbell/plots/NTCU/",patient,"_tree_with_hdp_signatures_branch_length_2023_08_20.pdf"))
  plot(tree,label.offset=0.01*max(tree_df$x))
  k <- 1
  for (k in 1:length(branches)){
    n=as.numeric(branches[k])
    row = rownames(tree_df)[tree_df$label == n]
    x_end=as.numeric(tree_df[row, "x"])
    parent = as.numeric(tree_df[row, "parent"])
    x_start=tree_df$x[parent]
    x_intv=x_end-x_start
    y=node.height(tree)[as.numeric(row)]
    tipnum=sum(tree_df$isTip)
    s <- "SBS1"
    for (s in sigs){
      # x_end=x_start+exposures[s,samples[k]]*x_intv
      x_end=x_start+sample_signatures[s,paste0(patient, "_Cl.", n)]*x_intv
      rect(ybottom=y-min(0.02*tipnum,0.2),ytop=y+min(0.02*tipnum,0.2),xleft=x_start,xright=x_end,col=cols[s])
      x_start=x_end
    }
  }
  nodelabels(text=tree_df[tree_df$isTip ==FALSE,]$label,node=tree_df[tree_df$isTip ==FALSE,]$node, col = "black", bg = "white", cex = 1, font = 7)
  axisPhylo(side = 1,backward=F)
  legend("topright",title="Signatures", legend=paste0(title.sigs), 
         fill=cols, bty="n",cex=0.8, ncol=1, xjust=0.5)
  dev.off()
  
}

#####################################################################################################
# COUNT THE SNVS AGAIN FROM ROOT TO TIP AFTER CORRECTION
#####################################################################################################

tree.dir <- "/Users/mp34/team154_campbell/plots/NTCU/Sensitivity_adjusted_trees/no_root"
tree_out_dir = "/Users/mp34/team154_campbell/plots/NTCU/Sensitivity_adjusted_trees/no_root/"
edge.files <- list.files(tree.dir, full.names = T, pattern = ".edge_tbl.csv", recursive = T, all.files = T)

# set the input dir to the NDP results
ndp.dir <- "/Users/mp34/sanger/team154pc/mp34/lung/NTCU/WGS/NDP"

# get clonal prevalence files
clonal.prev.files <- list.files(ndp.dir, pattern = "cluster_and_samples.csv", full.names = T, recursive = T)
clonal.prev.files <- clonal.prev.files[grep("lobe", clonal.prev.files)]
clonal.prev.files <- clonal.prev.files[-grep("snvs_indels", clonal.prev.files)]

script.dir = "/Users/mp34/sanger/team154pc/mp34/src/breast-ndp-tree-building/"
source(paste0(script.dir,"/ndp_tree_functions.R"))
source(paste0(script.dir,"/phylogeny_functions.R"))

i <- 1
for (i in 1:length(edge.files)) {
  
  patientID = str_split_fixed(basename(edge.files[i]), "\\.edge", 2)[,1]
  this_edgefile = fread(edge.files[i])
  
  if (ncol(this_edgefile) > 3) {
    this_edgefile <- this_edgefile %>%
      dplyr::select("From" = old_parent, "To" = old_child, Evidence)
  }
  
  this_branch_length <- fread(paste0(tree_out_dir, "/", patientID, "_branch_lengths.csv"))
  tr = convert_edges_to_phylo(tree_tbl =  this_edgefile, muts_per_cluster =  this_branch_length)
  
  ### save treefile
  write.tree(tr, file = paste(tree_out_dir, paste(patientID, 'cluster_tree.phylo',sep='.'),sep='/'))
  
  # this_tree <- convert_edges_to_phylo(this_edgefile, this_edgelength)
  # this_edgelength <- cbind(patientID, this_edgelength)
  # this_edgelength$cluster_id <- gsub("Cl.", "", this_edgelength$cluster_id)
  
  # this_edgelength <- fread(clonal.prev.files[grep(patientID, clonal.prev.files)])[,c("cluster_id", "no.of.mutations.assigned")]
  # colnames(this_edgelength) <- c("cluster_id", "num.muts")
  # 
  # this_edgelength <- fread(gsub(".edge_tbl.csv", "_branch_lengths.csv", edge.files[i]))
  # this_tree <- convert_edges_to_phylo(this_edgefile, this_edgelength)
  # this_edgelength <- cbind(patientID, this_edgelength)
  # this_edgelength$cluster_id <- gsub("Cl.", "", this_edgelength$cluster_id)
  
  tree_df=fortify(tr)
  cols=c(brewer.pal(length(sigs), "Paired"),"magenta","firebrick")
  
  # cols=c("grey80","peachpuff","forestgreen","firebrick","steelblue","pink2",
  #       "turquoise1",'orange2',"chartreuse") # "mediumorchid3","grey20",
  #       # "yellow2","mediumaquamarine","tomato4")
  names(cols)=sigs
  samples=colnames(sample_signatures)[grepl(patientID,colnames(sample_signatures))]
  samples=samples[grep("Cl.", samples)]
  branches= str_split_fixed(samples, "_Cl.", 2)[,2]
  branches <- branches[as.character(branches) %in% as.character(tree_df$label)]
  
  pdf(paste0(tree_out_dir, "/", patientID, "_tree_sensitivity_adjusted_filtered.pdf"))
  par(xpd=TRUE, mar=c(12.1, 4.1, 4.1, 2.1), family = "sans")
  # First draw a simple tree
  draw_nice_tree(tr, patientID=patientID) 
  dev.off()
  
  pdf(paste0("/Users/mp34/team154_campbell/plots/NTCU/",patientID,"_tree_with_hdp_signatures_branch_length_sensitivity_corrected_2023_08_22.pdf"))
  plot(tr,label.offset=0.01*max(tree_df$x))
  axisPhylo(side = 1,backward=F)
  legend("topright",title="Signatures", legend=paste0(title.sigs), 
         fill=cols, bty="n",cex=0.8, ncol=1, xjust=0.5)
  
  k <- 1
  for (k in 1:length(branches)){
    n=as.numeric(branches[k])
    row = rownames(tree_df)[tree_df$label == n]
    x_end=as.numeric(tree_df[row, "x"])
    parent = as.numeric(tree_df[row, "parent"])
    x_start=tree_df$x[parent]
    x_intv=x_end-x_start
    y=node.height(tr)[as.numeric(row)]
    tipnum=sum(tree_df$isTip)
    s <- "SBS1"
    for (s in sigs){
      # x_end=x_start+exposures[s,samples[k]]*x_intv
      x_end=x_start+sample_signatures[s,paste0(patientID, "_Cl.", n)]*x_intv
      rect(ybottom=y-min(0.02*tipnum,0.2),ytop=y+min(0.02*tipnum,0.2),xleft=x_start,xright=x_end,col=cols[s])
      x_start=x_end
    }
  }
  nodelabels(text=tree_df[tree_df$isTip ==FALSE,]$label,node=tree_df[tree_df$isTip ==FALSE,]$node, col = "black", bg = "white", cex = 1, font = 7)
  dev.off()
  

}




