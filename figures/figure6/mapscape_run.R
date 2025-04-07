args <- commandArgs(trailingOnly = T)

library(rmarkdown)
library(readxl)
library(stringr)
library(data.table)
library(dplyr)
library(funr)
library(phylobase)
library(RColorBrewer)
library(ape)
library(gplots)
library(wesanderson)

gitdir = funr::get_script_path()
gitdir = "/Users/mp34/team154_campbell/mapscape-generator"
setwd(gitdir)

## -- Receive command line args --
args <- commandArgs(trailingOnly = TRUE)
args


#### get command line arguments
# required input files:
#   1) sample name
# 2) static image file (.png)
# 3) table of pixel locations of microbiopsies on image
# 4) clonal prevalence in each microbiopsy
# 5) edge table of ape tree of clonal structure
# 6) whether low prevelence clones should be displayed, if not, tree filtered to only display clones present
# 7) output directory
# 8) (optional) table to mutations present in each clone
####

# Check for command line args, else set defaults
if (length(args)==7) {
  sample_name = args[1] # sample_name = "MD7047_Right_lobe"
  image.file = args[2] # image.file = "NTCU/MD7047_Right_lobe/MD7047_img.png"
  location.table.file = args[3] # location.table.file = "NTCU/MD7047_Right_lobe/MD7047_loc_tbl.xlsx"
  clonal.prevelance.file = args[4] # clonal.prevelance.file = "NTCU/MD7047_Right_lobe/cluster_and_samples.csv"
  tree.edge.file = args[5] # tree.edge.file = "NTCU/MD7047_Right_lobe/MD7047_Right_lobe.edge_tbl.csv"
  low_prev = args[6] # low_prev = "FALSE"
  output.dir = args[7] # output.dir = "NTCU/MD7047_Right_lobe/"
  message("Only 8 arguments, mutation file not provided")
} else if (length(args)==8) {
  sample_name = args[1]
  image.file = args[2]
  location.table.file = args[3]
  clonal.prevelance.file = args[4]
  tree.edge.file = args[5]
  low_prev = args[6]
  output.dir = args[7]
  mut.file = args[8] # mut.file = "NTCU/MD7047_Right_lobe/MD7047_ndp_assigned_muts_final.csv"
  message("Mutation File Provided")
} else {
  stop("Script takes either exactly 7 or 8 arguments.", call.=FALSE)
}

false_list <- c("f", "F", "False", "FALSE", "false")
true_list <- c("t", "T", "True", "TRUE", "true")

if (low_prev %in% false_list) {
  low_prev = FALSE
} else if (low_prev %in% true_list) {
  low_prev = TRUE
}

# Create output file
current_date = Sys.Date()
analysis.date = format(current_date, format = "%y%m%d")
out.file = paste0(output.dir, "/", sample_name, "_mapscape_", analysis.date, ".html")

#### Read and format inputs for running mapscape
if (exists("mut.file")) {
  print("mutation table found")
  
  location_tbl <- read_excel(location.table.file)
  included_samples <- location_tbl$sample_id
  
  pat_age <- unique(location_tbl$patient_id)
  
  edge_tbl <- fread(tree.edge.file) %>%
    dplyr::select("source" = old_parent, "target"= old_child)
  
  included_clones <- unique(edge_tbl$target)
  
  mut_tbl <- fread(mut.file)
  mut_tbl_short <- mut_tbl %>%
    dplyr::select(Chrom, Pos, Ref, Alt, Gene, Protein, Effect, vaf, sampleID, cluster_id, Driver)
  mut_tbl_short$cluster_id <- gsub("Cl.", "", mut_tbl_short$cluster_id)
  mut_tbl_short$sampleID <- str_sub(mut_tbl_short$sampleID, -8, -1)
  mut_tbl_short$sampleID <- gsub("lo", "", mut_tbl_short$sampleID)
  
  mut_tbl_final <- mut_tbl_short %>%
    dplyr::select("chrom" = Chrom, "coord" = Pos, "clone_id" = cluster_id, "sample_id"= sampleID, "VAF" = vaf, Gene, Protein, Effect, Driver) %>%
    dplyr::filter(sample_id %in% included_samples) %>%
    dplyr::filter(clone_id %in% included_clones) %>%
    dplyr::filter(Protein != "-")
  
  prev_tbl <- fread(clonal.prevelance.file)
  prev_tbl <- prev_tbl %>% filter(apply(prev_tbl[, c(1:(ncol(..prev_tbl)-2))], 1, function(x) any (x >= 0.1)))
  
  mut_counts <- prev_tbl[,c((ncol(..prev_tbl)-1):ncol(..prev_tbl))]
  colnames(mut_counts) <- c("Cluster ID", "Mutations Assigned")
  
  prev_tbl$cluster_id <- gsub("Cl.", "", prev_tbl$cluster_id)
  
  sample_names_long <- colnames(prev_tbl)

  # change sample name to 8 digit numbers for readability
  sample_names_short <- str_sub(sample_names_long[1:(length(sample_names_long)-2)], -8, -1)
  
  if(str_split_fixed(sample_name, "_", 2)[,1] == "MD6812"){
    
    sample_names_short[2:length(sample_names_short)] <- paste0("a", sample_names_short[2:length(sample_names_short)])
    
  } else if (str_split_fixed(sample_name, "_", 2)[,1] == "MD7047"){
    
    sample_names_short[2:length(sample_names_short)] <- paste0("d", sample_names_short[2:length(sample_names_short)])
  }

  sample_names_short <- gsub("lo", "", sample_names_short)
  colnames(prev_tbl) <- c(sample_names_short, "cluster_id", "mut_count")
  # drop last 2 columns, convert to long format
  prev_tbl_melt <- melt(prev_tbl[,1:(length(sample_names_long)-1)], id.vars = "cluster_id")
  colnames(prev_tbl_melt) <- c("clone_id", "sample_id", "clonal_prev")
  # remove samples not present in image (required by software)
  prev_tbl_filt <- prev_tbl_melt %>%
    dplyr::filter(sample_id %in% included_samples)
  prev_tbl_filt$clonal_prev <- format(prev_tbl_filt$clonal_prev, scientific = FALSE)
  prev_tbl_filt <- prev_tbl_filt[prev_tbl_filt$clonal_prev > 0.025, ]
  
} else {
  print("no mutation table found")
  location_tbl <- read_excel(location.table.file)
  included_samples <- location_tbl$sample_id
  edge_tbl <- fread(tree.edge.file) %>%
    dplyr::select("source" = old_parent, "target"= old_child)
  
  prev_tbl <- fread(clonal.prevelance.file)
  prev_tbl$cluster_id <- gsub("Cl.", "", prev_tbl$cluster_id)
  sample_names_long <- colnames(prev_tbl)
  # change sample name to 8 digit numbers for readability
  sample_names_short <- str_sub(sample_names_long[1:(length(sample_names_long)-2)], -8, -1)
  colnames(prev_tbl) <- c(sample_names_short, "cluster_id", "mut_count")
  # drop last 2 columns, convert to long format
  prev_tbl_melt <- melt(prev_tbl[,1:(length(sample_names_long)-1)], id.vars = "cluster_id")
  colnames(prev_tbl_melt) <- c("clone_id", "sample_id", "clonal_prev")
  prev_tbl_melt$sample_id <- as.character(prev_tbl_melt$sample_id)
  # remove samples not present in image (required by software)
  prev_tbl_filt <- prev_tbl_melt %>%
    dplyr::filter(sample_id %in% included_samples)
  prev_tbl_filt$clonal_prev <- format(prev_tbl_filt$clonal_prev, scientific = FALSE)
}

### get static tree image, assumed to be in same directory as H&E image, named first 7 characters of image_file + _static_tree.png
static_tree_file <- paste(dirname(image.file), paste0(substr(basename(image.file), 1, 7), "_SigFit_annotated.png"), sep = "/")

### get ndp_output_dir
ndp_dir <- dirname(clonal.prevelance.file)

### assign colors, each major clone (first division from root) and its descendant all get the same color scheme
# get tree
phylo_files <- list.files(dirname(tree.edge.file), pattern = "phylo", full.names = T, recursive = F)
ape.tree <- phylo_files[str_detect(phylo_files, substr(sample_name, 1, 7))]

tree <- read.tree(ape.tree)
tree4 <- as(tree, "phylo4")

# Get major clones
top_clones <- edge_tbl %>%
  dplyr::filter(source == "PD")

# ID which clones actually present in these samples (will get most distinct colours)
max_prev <- prev_tbl_filt %>%
  dplyr::group_by(clone_id) %>%
  dplyr::summarise(max_prev = max(clonal_prev)) %>%
  dplyr::filter(max_prev >= 0.1)
clones_present <- max_prev$clone_id

top_clones$present <- top_clones$target %in% clones_present

# Assign color for top clones, if 8 or fewer, distinct colours chosen, 
# if 8 or fewer are present, distince dark2 pallet used, rest randomly assigned
# if 12 of fewer are present, disticne Paried pallet used, rest randomly assign
# otherwise randomly assigned
if ( nrow(top_clones) < 3) {
  top_clone_col <- RColorBrewer::brewer.pal(8, "Dark2")
  top_clone_col <- sample(top_clone_col, nrow(top_clones))
} else if ( nrow(top_clones) >= 3 &&  nrow(top_clones) <= 8) {
  top_clone_col <- RColorBrewer::brewer.pal(nrow(top_clones), "Dark2")
  top_clone_col <- sample(top_clone_col, length(top_clone_col))
}  else if (sum(top_clones$present) < 3) {
  pres_col <- RColorBrewer::brewer.pal(8, "Dark2")
  pres_col <- sample(pres_col,sum(top_clones$present))
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(distinct = T), invert = T)]
  abs_col <- sample(color, nrow(top_clones)-length(pres_col) )
} else if (sum(top_clones$present) >= 3 && sum(top_clones$present) <= 8) {
  pres_col <- RColorBrewer::brewer.pal(sum(top_clones$present), "Dark2")
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(distinct = T), invert = T)]
  abs_col <- sample(color, nrow(top_clones)-length(pres_col) )
} else if (sum(top_clones$present) <= 12) {
  pres_col <- RColorBrewer::brewer.pal(sum(top_clones$present), "Paired")
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(distinct = T), invert = T)]
  abs_col <- sample(color, nrow(top_clones)-length(pres_col) )
} else {
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(distinct = T), invert = T)]
  top_clone_col <- sample(color, nrow(top_clones))
}

top_clones$top_col <- ""

if (exists("pres_col")) {
  top_clones[top_clones$present == TRUE,]$top_col <- sample(col2hex(pres_col), length(pres_col), replace = F)
  top_clones[top_clones$present == FALSE,]$top_col <- col2hex(abs_col)
} else {
  top_clones$top_col <- col2hex(top_clone_col)
}

clones <- top_clones$target
# Get descendants for top clones, ramp colour towards white for descendante
# For clones with 9 or more descedants second colour is added otherwise too difficult to distinguish
# For clones with 15 or more descendants, third colour is added
color_list <- list()
for (i in 1:length(clones)) {
  this_clone =  clones[i]
  this_color = top_clones[top_clones$target == this_clone, top_col]
  nodes <- phylobase::nodeLabels(tree4)
  node_tbl <- tibble(node = names(nodes), label = nodes)
  new_node <- node_tbl[node_tbl$label == this_clone,]$node 
  new_descendants <- descendants(tree4, as.numeric(new_node), type = "ALL")
  if (length(new_descendants) >= 1 && length(new_descendants) <= 9) {
    branch_clones <- names(new_descendants)
    ramp_colors <- colorRampPalette(c(this_color, "white"))
    clone_colors <- ramp_colors(length(branch_clones)+1)
    clone_colors <- clone_colors[-(length(branch_clones)+1)]
    color_tbl <- tibble(clone_id = branch_clones, colour = clone_colors )
    color_list[[i]] <- color_tbl
  } else if (length(new_descendants) >= 10 && length(new_descendants) < 15) {
    branch_clones <- names(new_descendants)
    # ramp towards white for half of clones, towards black for the other half
    # black ramp is then reveserved and combined with white ramp
    ramp_white <- colorRampPalette(c(this_color, "white"))
    ramp_grey <- colorRampPalette(c(this_color, "black"))
    clone_white <- ramp_white(round(length(branch_clones)/2)+1)
    clone_white <- clone_white[-(length(clone_white))]
    
    clone_grey <- ramp_grey(round(length(branch_clones)/2)+6)
    clone_grey <- clone_grey[-c((length(clone_grey)-3):(length(clone_grey)))]
    
    remain_cols <- length(new_descendants) - length(clone_white)
    
    clone_colors <- c(clone_white, clone_grey[1:remain_cols])
    
    color_tbl <- tibble(clone_id = branch_clones, colour = clone_colors )
    color_list[[i]] <- color_tbl
  } else if (length(new_descendants) >= 15) {
    childs <- descendants(tree4, as.numeric(new_node), type = "children") # get first level of descendants for new colour scheme
    child_clones <- names(childs)
    if (length(childs) <= 10) {
      pal1 <- wes_palette("Darjeeling1")
      pal2 <- wes_palette("BottleRocket2")
      pals <- c(pal1, pal2)
      child_col <- sample(pals, length(childs))
      child_col <- col2hex(child_col)
    } else {
      color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(distinct = T), invert = T)]
      child_col <- sample(color, length(childs))
      child_col <- col2hex(child_col)
    }
    parent_tbl <- tibble(clone_id = this_clone, colour = this_color)
    child_tbl <- tibble(clone_id = child_clones, colour = child_col)
    child_tbl <- rbind(parent_tbl, child_tbl)
    child_list <- list()
    for (j in 1:length(child_clones)) {
      this_child = child_clones[j]
      this_color = child_tbl[child_tbl$clone_id == this_child,]$colour
      this_child_node = node_tbl[node_tbl$label == this_child ,]$node 
      if (length(this_child_node)>0) {
        child_descendants <- descendants(tree4, as.numeric(this_child_node), type = "all")
        branch_children <-  names(child_descendants)
        ramp_colors <- colorRampPalette(c(this_color, "white"))
        clone_colors <- ramp_colors(length(branch_children)+2)
        clone_colors <- clone_colors[-c(1,(length(branch_children)+2))]
        child_tibble <- tibble(clone_id = branch_children, colour = clone_colors )
        child_list[[j]] <- child_tibble
      }
      
    }
    child_color_table <- bind_rows(child_list)
    child_color_table <- rbind(child_color_table, child_tbl)
    color_list[[i]] <- child_color_table
    } else {
    color_tbl <- tibble(clone_id = as.character(this_clone), colour = this_color)
    color_list[[i]] <- color_tbl
  }
}

clone_colour_tbl <- bind_rows(color_list)
PD_tbl <- tibble(clone_id = "PD", colour = "000000")
clone_colour_tbl <- bind_rows(clone_colour_tbl, PD_tbl)
clone_colour_tbl$clone_id <- as.character(clone_colour_tbl$clone_id)

color_tbl = clone_colour_tbl

# SigFit.dir = "~/nfs_tb14/scratch_117/Breast_Project/SigFit/220117/"

# Run Mapscape_Generation.Rmd
if (exists("mut_tbl_final")) {
  rmarkdown::render(input = "Mapscape_Generation.Rmd", output_format = "html_document",  
                    output_file = out.file, params = list(prev_tbl_filt = prev_tbl_filt,
                                                          edge_tbl = edge_tbl,
                                                          location_tbl = location_tbl,
                                                          image.file = image.file,
                                                          static_tree_file = static_tree_file,
                                                          sample_name = sample_name,
                                                          analysis.date = analysis.date,
                                                          low_prev = low_prev,
                                                          color_tbl = color_tbl,
                                                          mut_tbl_final = mut_tbl_final,
                                                          gitdir = gitdir,
                                                          pat_age = pat_age,
                                                          ndp_dir = ndp_dir,
                                                          mut_counts = mut_counts,
                                                          SigFit.dir = SigFit.dir))
  
} else {
  rmarkdown::render(input = "Mapscape_Generation.Rmd", output_format = "html_document",  
                    output_file = out.file, params = list(prev_tbl_filt = prev_tbl_filt,
                                                          edge_tbl = edge_tbl,
                                                          location_tbl = location_tbl,
                                                          image.file = image.file,
                                                          static_tree_file = static_tree_file,
                                                          sample_name = sample_name,
                                                          analysis.date = analysis.date,
                                                          low_prev = low_prev,
                                                          color_tbl = color_tbl,
                                                          gitdir = gitdir,
                                                          pat_age = pat_age,
                                                          ndp_dir = ndp_dir,
                                                          mut_counts = mut_counts))

}
