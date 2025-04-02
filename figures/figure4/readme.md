# Figure 4 -  Epithelial cell shift during early human lung squamous cell carcinogenesis.

**Human_scRNAseq_01_preprocessing.R** - Script to perform the preprocessing of human airway scRNA-Seq datasets. It includes data integration steps for merging multiple datasets, normalisation, and QC filtering to ensure high-quality data for downstream analysis.

**Human_scRNAseq_02_Signature_overlap.R** - Script tp perform a signature overlap analysis, comparing gene expression signatures between the integrated human airway dataset and previous studies of human epithelial cells.

**Human_scRNAseq_03_data_visulisation.R** - Script to generate various visualisations, including UMAP plots, heatmaps, and marker dot plots, to visualise the expression patterns and cell clusters in the human airway datasets.

**Human_scRNAseq_04_Milo_cell_abundance.R** - Script to use the Milo package to assess differential cell abundance across experimental conditions. It performs a statistical analysis to identify cell types or subclusters that are significantly enriched or depleted.

**Human_scRNAseq_05_Slingshot_trajectory.R** - Script to apply the Slingshot package to infer developmental trajectories in human airway epithelial cells. It models the pseudotemporal progression of the cells based on their gene expression profiles, highlights changes between current-smokers vs non-somkers.

