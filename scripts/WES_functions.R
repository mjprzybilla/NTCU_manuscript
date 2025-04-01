
# Function to process VCF files for mutation analysis
process_vcf <- function(clonality) {
  SampleInfo$files <- paste0(VCF.path, "/", SampleInfo$Tumor_Sample_Barcode, clonality, ".vcf")
  
  vcf_files <- SampleInfo$files
  sample_names <- SampleInfo$SAMPLE
  sample_Status <- SampleInfo$Smoking_status
  
  snv_grl <- read_vcfs_as_granges(vcf_files, sample_names, "BSgenome.Hsapiens.UCSC.hg19", type = "snv", predefined_dbs_mbs = TRUE)
  
  # Type occurrences
  type_occurrences <- mut_type_occurrences(snv_grl, "BSgenome.Hsapiens.UCSC.hg19")
  fwrite(type_occurrences, paste0(clonality, ".SNV.type_occurrences.txt"), sep="\t")
  
  # Spectrum plot
  plot_spectrum(type_occurrences, by = sample_Status, CT = TRUE, legend = TRUE) %>%
    ggsave(paste0(clonality, ".SNV.spectrum.splitbystatus.pdf"), width=7, height=6)
  
  # Mutation matrix
  mut_mat <- mut_matrix(vcf_list = snv_grl, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
  
  # 96-profile plot
  plot_96_profile(mut_mat, condensed = FALSE, ymax = 0.05) %>%
    ggsave(paste0(clonality, ".SNV.profiles.plot.pdf"), width=7, height=14)
  
  # Fit known signatures
  signatures <- get_known_signatures()
  fit_res <- fit_to_signatures(mut_mat, signatures)
  
  # Selected top signatures
  sigs.to.analyse <- c('SBS4','SBS92','SBS13','SBS2','SBS5','SBS1')
  cancer_signatures.used <- signatures[, sigs.to.analyse]
  fit_res_selected <- fit_to_signatures(mut_mat, cancer_signatures.used)
  
  # Contribution plots
  plot_contribution(fit_res_selected$contribution, coord_flip = FALSE, mode = "relative") +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) %>%
    ggsave(paste0(clonality, ".SBS.signatures.relative.pdf"), width=14, height=7)
  
  plot_contribution(fit_res_selected$contribution, coord_flip = FALSE, mode = "absolute") +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) %>%
    ggsave(paste0(clonality, ".SBS.signatures.absolute.pdf"), width=14, height=7)
  
  # Save top signatures
  top_sigs <- as.data.frame(fit_res_selected$contribution) %>%
    rownames_to_column("signature") %>%
    pivot_longer(-signature, names_to = "SampleID", values_to = "Sig_Contribution")
  
  fwrite(top_sigs, paste0(clonality, ".top.signatures.txt"), sep="\t")
  
  # Relative signature calculation
  sample_sum <- top_sigs %>%
    group_by(SampleID) %>%
    summarise(Sig_Contribution_sum = sum(Sig_Contribution))
  
  top_sigs <- top_sigs %>%
    left_join(sample_sum, by = "SampleID") %>%
    mutate(Relative_Contribution = Sig_Contribution / Sig_Contribution_sum)
  
  fwrite(top_sigs, paste0("Top5.", clonality, ".signatures.tsv"), sep="\t")
  
  return(top_sigs)
}

# Function to get the seeding colors
get.seeding.colors <- function (fullTreeOutput, seedingTable) {
    print(fullTreeOutput$tumour_id)
    
    tree <- matrix(as.numeric(fullTreeOutput$graph_pyclone$Corrected_tree), ncol = 2)
    trunk <- fullTreeOutput$graph_pyclone$trunk
    clonalities <- fullTreeOutput$clonality_out$clonality_table_corrected
    ccfTable <- fullTreeOutput$nested_pyclone$ccf_cluster_table
    
    clusters <- sort(unique(as.numeric(tree)))
    regions <- colnames(clonalities)
    invasive.timepoint <- grep("INV", regions, value = TRUE)
    preinvasive.timepoints <- grep("INV", regions, value = TRUE, invert = TRUE)
    
    clonalities <- clonalities[rownames(clonalities) %in% clusters, , drop = FALSE]
    ccfTable <- ccfTable[rownames(ccfTable) %in% clusters, , drop = FALSE]
    ccfTable[which(clonalities == "clonal")] <- 100
    
    clonalities <- data.frame(clonalities, stringsAsFactors = FALSE) %>% mutate(clones = rownames(.))
    ccfTable <- data.frame(ccfTable, stringsAsFactors = FALSE) %>% mutate(clones = rownames(.))
    
    overviewClonalities <- clonalities %>%
        mutate(tumour_id = fullTreeOutput$tumour_id) %>%
        select(tumour_id, clones) %>%
        mutate(INVClonality = get.tumLevel.clonality(clonalities, timepoints.to.use = invasive.timepoint)) %>%
        mutate(PREClonality = get.tumLevel.clonality(clonalities, timepoints.to.use = preinvasive.timepoints)) %>%
        mutate(colors = ifelse(INVClonality == "clonal" | PREClonality == "clonal", "black", "transparent"))
    
    seedingClusters <- seedingTable %>% 
        filter(tumour_id == fullTreeOutput$tumour_id) %>% 
        pull(seedingCluster) %>%
        unique()
    
    treeStructure <- get.treeStructure(tree, trunk)
    treeStructureOverview <- lapply(strsplit(treeStructure$structure, split = ":"), function(y) {
        return(setNames(y, paste0("pos", 1:length(y))))
    })
    fullClusterOrder <- as.character(unlist(apply(bind_rows(treeStructureOverview), 2, unique)))
    
    seedingClusters <- seedingClusters[order(match(seedingClusters, fullClusterOrder))]
    
    if (length(invasive.timepoint) >= 1) {                              
        INV.seeding.cluster <- seedingTable %>%
            filter(grepl('INV', sample_id)) %>%
            pull(seedingCluster) %>%  
            unique() %>%
            sort()
        
        print(INV.seeding.cluster)
        seedingClusters <- seedingClusters[!seedingClusters %in% INV.seeding.cluster]
        seedingClusters <- c(INV.seeding.cluster, seedingClusters)
        
        overviewClonalities$stroke <- 0.5
        tmp <- overviewClonalities$clones[which(overviewClonalities$INVClonality != "absent")]
        overviewClonalities$stroke[overviewClonalities$clones %in% tmp] <- 3
    }
    
    seeding.clusters.colors <- brewer.pal(7, "Set1")
    
    if (length(invasive.timepoint) < 1) {                              
        seeding.clusters.colors <- seeding.clusters.colors[-1]
        overviewClonalities$stroke <- 0.5
    }
    
    print(seeding.clusters.colors)
    col.palette <- setNames(seeding.clusters.colors, seedingClusters)
    
    for (x in seedingClusters) {
        overviewClonalities$colors[which(overviewClonalities$clones == x)] <- col.palette[as.character(x)]
        
        seedingCloneBranches <- lapply(strsplit(treeStructure$structure, split = ":"), function(y) {
            if (any(grepl(paste0("^", x, "$"), y))) {
                foo <- y[seq(grep(paste0("^", x, "$"), y), length(y))]
                return(setNames(foo, paste0("pos", 1:length(foo))))
            }
            return(NULL)
        })
        seedingCloneBranches <- seedingCloneBranches[!sapply(seedingCloneBranches, is.null)]
        
        col.palette.tmp <- colorRampPalette(c(col.palette[as.character(x)], "white"))(length(unique(unlist(seedingCloneBranches))))
        clusterOrder <- as.character(unlist(apply(bind_rows(seedingCloneBranches), 2, unique)))
        clusterOrder <- clusterOrder[!is.na(clusterOrder)]
        
        overviewClonalities$colors[match(clusterOrder, overviewClonalities$clones)] <- col.palette.tmp
    }
    
    return(overviewClonalities)
}

# Function to work out the clonality
get.tumLevel.clonality <- function(clonalities, timepoints.to.use = NULL) {
        ### only use timepoints of interest
        ### if timepoints.to.use == NULL then assume using all timepoints in the clonality table
        if (!is.null(timepoints.to.use)) {
            if (length(timepoints.to.use) == 0) {
                tmp <- rep("NA", nrow(clonalities))
                names(tmp) <- rownames(clonalities)
                return(tmp)
            }
            clonalities <- clonalities[timepoints.to.use]
        }
        ### sanity checks
        if (class(clonalities) != "data.frame") error("Clonalities have wrong format.")
        if (nrow(clonalities) == 0) return(NULL)

        tmp <- sapply(seq(1, nrow(clonalities)), function(i) {
            if (all(clonalities[i,] == "clonal")) return("clonal")
            if (any(clonalities[i,] == "subclonal")) return("subclonal")
            if (all(clonalities[i,] == "absent")) return("absent")
            if (any(clonalities[i,] == "absent")) return("subclonal")
            return(NA)
        })
        names(tmp) <- rownames(clonalities)
        return(tmp)
} 

# Function to get the tree structure
get.treeStructure <- function(tree, trunk) {
    tree.structure.final <- c()
    tree.trunk.final     <- as.numeric(trunk)
    treeInfo.tmp         <- tree
    
    # If there is only one node and it matches the trunk, return the trunk as both the trunk and tips
    if (nrow(treeInfo.tmp) == 1 & tree.trunk.final == treeInfo.tmp[1, 1] & tree.trunk.final == treeInfo.tmp[1, 2]) {
        return(list(
            trunk = tree.trunk.final,
            tips = tree.trunk.final,
            structure = tree.trunk.final,
            allTreeClones = tree.trunk.final
        ))
    }
    
    allNodes             <- unique(as.numeric(tree))
    tree.tips.final      <- allNodes[which(!allNodes %in% treeInfo.tmp[, 1])]
    
    while(length(which(treeInfo.tmp[, 1] == tree.trunk.final)) != 0) {
        idx <- idx.rm <- min(which(treeInfo.tmp[, 1] == tree.trunk.final))
        node <- treeInfo.tmp[idx, 2]
        
        # Traverse the tree until a tip is found
        while (!node %in% tree.tips.final) {
            idx <- c(idx, min(which(treeInfo.tmp[, 1] == node)))
            if (length(which(treeInfo.tmp[, 1] == node)) > 1) {
                idx.rm <- idx[length(idx)]
            } else {
                idx.rm <- c(idx.rm, min(which(treeInfo.tmp[, 1] == node)))
            }
            node <- treeInfo.tmp[min(which(treeInfo.tmp[, 1] == node)), 2]
        }
        
        # Collapse by row to simplify removal of duplicate nodes
        tree.structure.final <- c(tree.structure.final, paste0(
            apply(treeInfo.tmp[idx, , drop = FALSE], 1, paste0, collapse = ":"), collapse = ":"
        ))
        
        # Remove branches already covered
        treeInfo.tmp <- treeInfo.tmp[-idx.rm, , drop = FALSE]
    }
    
    # Remove duplicate nodes
    tree.structure.final <- unlist(lapply(lapply(strsplit(tree.structure.final, split = ":"), unique), paste0, collapse = ":"))
    
    return(list(
        trunk = tree.trunk.final,
        tips = tree.tips.final,
        structure = tree.structure.final,
        allTreeClones = allNodes
    ))
}

# Function to get data from tree
getDataForTree <- function(tree, trunk, node.fill = "gray", node.shape = 21, node.size = 13, stroke.size = 1) {
    # Add the trunk node to the tree
    tree <- rbind(c("0", trunk), tree)

    # Create an undirected graph from the tree data
    g.tree <- igraph::as.undirected(igraph::graph.data.frame(tree))

    # Get the names of the nodes
    node.indx <- igraph::V(g.tree)$name

    # Layout the tree
    l.tree <- igraph::layout_as_tree(g.tree, root = '0', flip.y = FALSE)

    # Prepare the nodes for plotting
    tree.plot.nodes <- data.frame(
        Node = node.indx[-1], 
        x = l.tree[-1, 1], 
        y = l.tree[-1, 2],
        size = node.size,
        shape = node.shape,
        fill = node.fill,
        color = "black",
        stroke = stroke.size,
        node.label = "", 
        stringsAsFactors = FALSE
    )

    # Prepare the edges for plotting
    tree.plot.edges <- data.frame(
        Edge = apply(tree, 1, function(x) paste0(x, collapse = "-")),
        xstart = l.tree[match(tree[, 1], node.indx), 1],
        xend = l.tree[match(tree[, 2], node.indx), 1],
        ystart = l.tree[match(tree[, 1], node.indx), 2],
        yend = l.tree[match(tree[, 2], node.indx), 2],
        color = "black",
        size = stroke.size, 
        stringsAsFactors = FALSE
    )

    # Set up additional edge attributes
    tree.plot.edges$edge.label <- ""
    tree.plot.edges$label.x <- (tree.plot.edges$xstart + tree.plot.edges$xend) / 2
    tree.plot.edges$label.y <- (tree.plot.edges$ystart + tree.plot.edges$yend) / 2
    tree.plot.edges$highlight.color <- "transparent"

    # Return the nodes and edges data for plotting
    return(list(
        tree.plot.nodes = tree.plot.nodes, 
        tree.plot.edges = tree.plot.edges
    ))
}

# Function to plot tree
plottingTxTree <- function(tree.plot.nodes, tree.plot.edges, tree.title = NULL, tree.subtitle = NULL, node.label.type = "none", add.edge.labels = FALSE, edge.highlights = NULL) {
    # Start with an empty ggplot object
    p <- ggplot()

    # If edge highlights are provided, add them to the plot
    if (!is.null(edge.highlights)) {
        p <- p + geom_segment(
            data = tree.plot.edges, 
            aes(x = xstart, y = ystart, xend = xend, yend = yend),
            color = tree.plot.edges$highlight.color, 
            linewidth = tree.plot.edges$size * 5
        )
    }

    # Plot the edges and nodes, including node styling
    p <- p + 
        geom_segment(
            data = tree.plot.edges, 
            aes(x = xstart, y = ystart, xend = xend, yend = yend), 
            color = tree.plot.edges$color, 
            linewidth = tree.plot.edges$size
        ) +
        geom_point(
            data = tree.plot.nodes, 
            aes(x = x, y = y), 
            fill = "white", color = "black", 
            size = tree.plot.nodes$size, 
            shape = tree.plot.nodes$shape, 
            stroke = tree.plot.nodes$stroke
        ) +
        geom_point(
            data = tree.plot.nodes, 
            aes(x = x, y = y), 
            fill = tree.plot.nodes$fill, 
            color = tree.plot.nodes$color, 
            size = tree.plot.nodes$size, 
            shape = tree.plot.nodes$shape, 
            stroke = tree.plot.nodes$stroke
        ) +
        theme_void()  # Remove background and axes

    # Add title and subtitle if provided
    if (!is.null(tree.title)) {
        p <- p + ggtitle(tree.title) + theme(plot.title = element_text(hjust = 0.5))
    }

    if (!is.null(tree.subtitle)) {
        p <- p + labs(subtitle = tree.subtitle) + theme(plot.subtitle = element_text(hjust = 0.5))
    }

    # If edge labels should be added, plot them
    if (add.edge.labels) {
        p <- p + geom_text(
            data = tree.plot.edges, 
            aes(x = label.x, y = label.y, label = stringr::str_wrap(label, 20)), 
            hjust = 0, nudge_x = 0.05, vjust = 1, size = 3
        )
    }

    # Handle node labels based on the specified type
    if (node.label.type == "none") {
        return(p)  # Return the plot without node labels
    } else if (node.label.type == "nodeNumber") {
        p <- p + geom_text(
            data = tree.plot.nodes, 
            aes(x = x, y = y, label = Node), 
            size = 7
        )
        return(p)
    } else if (node.label.type == "custom") {
        p <- p + geom_point(
            data = tree.plot.nodes %>% filter(node.label != ""), 
            aes(x = x, y = y), 
            fill = "white", 
            color = "white", 
            shape = 24, 
            size = tree.plot.nodes %>% filter(node.label != "") %>% pull(size) * 0.5
        )
        return(p)
    } else {
        warning(paste0("Node label type ", node.label.type, " unknown. Returning basic plot."))
        return(p)
    }
}

# Function to highlight nodes while plotting tree
plottingTimepointTreeNodeHighlight <- function(fullTreeOutput, colorTable, colorBy = "timepointAll", timepoints.to.use = NULL) {
        sampleID <- fullTreeOutput$sample_id
        tumorID     <- fullTreeOutput$tumour_id
        tree        <- fullTreeOutput$graph_pyclone$Corrected_tree
        trunk       <- fullTreeOutput$graph_pyclone$trunk
        edgelength  <- fullTreeOutput$graph_pyclone$edgelength

        tmp <- getDataForTree(tree, trunk)
        tree.plot.nodes <- tmp$tree.plot.nodes
        tree.plot.edges <- tmp$tree.plot.edges
        rm(tmp)

        colorTable <- colorTable %>% filter(tumour_id == tumorID) %>% rename(color = colors)

        if (colorBy == "timepoints") {
            if (is.null(timepoints.to.use)) {
                print("No timepoints specified. Using tissueAll")
                colorBy <- "tissueAll"
            }
        }

        if (colorBy == "timepointAll") {
            plot.title.label <- tumorID
            colorTable <- colorTable %>% mutate(line.color = "black") %>% mutate(node.size = 13)
            tmp <- get.tumLevel.clonality(fullTreeOutput$clonality_out$clonality_table_corrected)
            colorTable <- colorTable %>% left_join(data.frame(clones = names(tmp), timepointsOfInterest = as.character(tmp), stringsAsFactors = FALSE), by = "clones")
            colorTable <- colorTable %>% mutate(color = ifelse(timepointsOfInterest %in% c("clonal", "subclonal"), as.character(color), "#FFFFFF")) %>%
                                           mutate(line.color = ifelse(timepointsOfInterest %in% c("clonal", "subclonal"), as.character(line.color), "gray")) %>%
                                           mutate(node.size = ifelse(timepointsOfInterest %in% c("clonal", "subclonal"), as.numeric(as.character(node.size)), 6))

            
        } else if (colorBy == "timepoints") {
            plot.title.label <- sampleID
            colorTable <- colorTable %>% mutate(line.color = "black") %>% mutate(node.size = 13)
            tmp         <- get.tumLevel.clonality(fullTreeOutput$clonality_out$clonality_table_corrected, timepoints.to.use)
            colorTable <- colorTable %>% left_join(data.frame(clones = names(tmp), timepointsOfInterest = as.character(tmp), stringsAsFactors = FALSE), by = "clones")
            colorTable <- colorTable %>% mutate(color = ifelse(timepointsOfInterest %in% c("clonal", "subclonal"), as.character(color), "#FFFFFF")) %>%
                                           mutate(line.color = ifelse(timepointsOfInterest %in% c("clonal", "subclonal"), as.character(line.color), "gray")) %>%
                                           mutate(node.size = ifelse(timepointsOfInterest %in% c("clonal", "subclonal"), as.numeric(as.character(node.size)), 6))
        }
    
        tree.plot.nodes <- tree.plot.nodes %>% 
            select(-fill, -color, -size, -stroke) %>% 
            left_join(colorTable %>% 
                select(Node = clones, fill = color, color = line.color, size = node.size, stroke) %>% 
                unique(), by = "Node") %>% 
            select(Node, x, y, size, shape, fill, color, stroke, node.label)

        tree.plot.edges <- tree.plot.edges %>% 
            mutate(daughter = sapply(strsplit(Edge, split = "-"), function(x) x[2])) %>%
            select(-color) %>%
            left_join(colorTable %>% 
                select(daughter = clones, color = line.color) %>%
                unique(), by = "daughter") %>%
            select(Edge, xstart, xend, ystart, yend, color, size, edge.label, label.x, label.y, highlight.color)

        # pat.df <- subset(mutTableAll, mutTableAll$tumour_id == tumorID)
        # driverTable <- subset(pat.df, pat.df$DriverMut == TRUE)
        # driverTable <- driverTable[, c("PyCloneCluster_SC", "Hugo_Symbol")]
        # ### remove drivers that aren't clustered
        # driverTable <- subset(driverTable, !is.na(driverTable$PyCloneCluster_SC))
        # ### remove drivers that aren't in tree clusters
        # driverTable <- subset(driverTable, driverTable$PyCloneCluster_SC %in% unique(c(tree)))
        # driverTable$Edge <- sapply(driverTable$PyCloneCluster_SC, function(x) {
        #     if (x == trunk) return(paste0("0-", trunk))
        #     return(paste0(tree[grep(paste0("^", x, "$"), tree[, 2]), 1], "-", x))
        # })
        # driverTable <- driverTable %>% group_by(PyCloneCluster_SC) %>% mutate(label = paste0(Hugo_Symbol, collapse = ", ")) %>% select(PyCloneCluster_SC, label) %>% unique()
        # node.annotation <- setNames(rep("*", nrow(driverTable)), driverTable$PyCloneCluster_SC)

        p <- plottingTxTree(tree.plot.nodes, 
                            tree.plot.edges, 
                            tree.title = plot.title.label, 
                            tree.subtitle = NULL, 
                            node.label.type = "none", 
                            add.edge.labels = FALSE, 
                            edge.highlights = NULL)
        return(p)
}

# Function to plot the CloneMaps for each sample
plottingCloneMaps <- function(fullTreeOutput, cloneColorDF, outDir = "") {
        tumour_id       <- fullTreeOutput$tumour_id
        tree            <- fullTreeOutput$graph_pyclone$Corrected_tree
        clonalities     <- fullTreeOutput$clonality_out$clonality_table_corrected
        ccfTable        <- fullTreeOutput$nested_pyclone$ccf_cluster_table

        ### change tree to numeric
        tree <- matrix(as.numeric(tree), ncol = 2)

        clusters        <- sort(unique(as.numeric(tree)))
        timepoints         <- colnames(clonalities) 
        invasive.timepoint <- grep("INV", timepoints, value = TRUE)
        preinvasiveTP     <- grep("INV", timepoints, value = TRUE, invert = TRUE)


        clonalities     <- clonalities[rownames(clonalities) %in% clusters, , drop = FALSE]
        ccfTable        <- ccfTable[rownames(ccfTable) %in% clusters, , drop = FALSE]
       
        ### adjust clonal CCF to 100
        ccfTable[which(clonalities == "clonal")] <- 100

        clonalities     <- data.frame(clonalities, stringsAsFactors = FALSE)
        ccfTable        <- data.frame(ccfTable, stringsAsFactors = FALSE)

        clonalities     <- clonalities %>% mutate(clones = rownames(.))
        ccfTable        <- ccfTable %>% mutate(clones = rownames(.))

        overviewClonalities <- cloneColorDF %>% filter(tumour_id == !!tumour_id)

        col.palette     <- setNames(overviewClonalities$colors, overviewClonalities$clones)

        for (region in colnames(ccfTable)[-which(colnames(ccfTable) == "clones")]) {
            print(region)
            regionCCF <- ccfTable[, c("clones", region)]
            colnames(regionCCF)[2] <- "CCF"
            regionCCFfix <- make.CCFs.tree.consistant(tree, regionCCF, increase.parents = TRUE)
            pdf(paste0(outDir, region, ".pdf"), useDingbats = FALSE)
            cloneMap(tree, regionCCFfix, border.thickness = 3, clone.cols = col.palette, output.Clone.map.obj = FALSE)
            dev.off()
        }
}


# Plotting dn_ds
plot_dn_ds <- function(df.sel_cv, colors) {
  dn.ds.plot <- ggplot(df.sel_cv, aes(x = Truncal, y = Subclonal, fill = Selection)) +
    geom_point(aes(size = Sites, stroke = stroke), col = "black", shape = 21) +
    scale_size(range = c(3, 10)) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    guides(fill = guide_legend(override.aes = list(size = 5))) +
    scale_x_continuous(limits = c(-4000, 20000), breaks = c(1, 10000, 20000)) +
    scale_y_continuous(limits = c(-500, 2000), breaks = c(1, 1000, 2000)) +
    geom_hline(yintercept = 1, color = "grey") +
    geom_vline(xintercept = 1, color = "grey") +
    
    # Label truncal genes at the bottom
    geom_text_repel(data = subset(df.sel_cv, Selection == "Truncal favoured"), 
                    aes(label = gene_name), 
                    force_pull = 0.5,
                    nudge_y = -400 ,
                    direction = "x",
                    angle = 90,
                    hjust = 2,
                    segment.size = 0.2,
                    size = 6,
                    segment.color = "black") +
    
    # Label other genes at the top
    geom_text_repel(data = subset(df.sel_cv, Selection == "Subclonal favoured"), 
                    aes(label = gene_name), 
                    force_pull = 0.5,
                    nudge_y = 400,
                    direction = "y",
                    angle = 0,
                    hjust = 2,
                    segment.size = 0.2,
                    size = 6,
                    segment.color = "black") +

    # Label both genes at the top
    geom_text_repel(data = subset(df.sel_cv, Selection == "Subclonal and truncal"), 
                    aes(label = gene_name),
                    force_pull = 0.5,
                    nudge_x = -400,
                    direction = "y",
                    hjust = 0,
                    vjust = -2,
                    segment.size = 0.5,
                    size = 6,
                    segment.color = "black") +
    cowplot::theme_cowplot(font_size = 22) +
    labs(size = "Regions (%)", fill = "", x = "dN/dS truncal", y = "dN/dS subclonal")
  
  # Save plot
  ggsave(file = "dn.ds.truncal.clonal.pdf", plot = dn.ds.plot, width = 10, height = 6)
}
