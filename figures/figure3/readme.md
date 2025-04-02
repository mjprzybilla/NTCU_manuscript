# Figure 3 -  Single cell profiling reveals an epithelial cell fate shift following NTCU treatment.

**scRNA_01_scanpy_NTCU_GEX_analysis.ipynb** - Script to perform standard analysis workflow of single-cell RNA-sequencing data from the trachea.

**Mouse_scRNAseq_02_Epithelial_subclustering.R** - Script to perform subclustering of epithelial cells identified from the trachea scRNA-Seq dataset.

**Mouse_scRNAseq_03_data_visualization.R** - Script to generate various visualisations, including UMAP plots, marker dotplots, and cell proportion bar plots. These visualisations provide insights into the cellular heterogeneity within the epithelial population.

**Mouse_scRNAseq_04_Epithelial_Signature_overlap.R** - Script to compare the epithelial cell signatures derived from the current dataset with previously published signatures from mouse epithelial studies. It uses gene set overrepresentation analysis (ORA) to explore similarity between signatures.

**Mouse_scRNAseq_05_Epithelial_scCODA.ipynb** - Script to focuses on the implementation of Bayesian model for compositional single-cell data analysis to identify cell type changes between NTCU model and Control subjects.

**Mouse_scRNAseq_06_Epithelial_Monocle_trajectory.R** - Script to uses the Monocle package to perform comparative trajectory analysis on the epithelial cells from trachea between NTCU and Control.
