

################################################################################################################################################
##                                                                                                                      
##  HUMAN WES data DATA - Visualize of CONIPHER trees
##                                                                                                                      
##  Date: 18 july 2024                                                                                                                   
##  
##  Author: Ahmed Alhendi                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################


# Load custom plotting functions
source("scripts/WES_functions.R")
source("scripts/cloneMap.R")

# List of required packages
library(dplyr)
library(data.table)
library(tibble)
library(stringr)
library(CONIPHER)
library(ggrepel)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(ggsci)
library(RColorBrewer) 

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

