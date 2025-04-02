# Figure 7 -  Clonal relatedness between anatomically distinct human preinvasive lung lesions.

**Human_WES_01_Germline_Correlation.R** - Script to generate correlation plot based on germline mutation data. It calculates pairwise correlations between different samples to explore similarity or dissimilarity in mutation profiles (Figure S13).

**Human_WES_02_CONIPHER.R** - Script to use CONIPHER to prefom hierarchical clustering of samples of multiple annotomical sites and visualises the relationships between them using tree-based structures.

**Human_WES_03_Mutational_signatures.R** - Script to perform mutational signature analysis, identifying distinct mutational patterns within the WES data from CIS samples (Figure S15). Results were sorted by clonality drived from CONIPHER.

**Human_WES_04_dndscv.R** - Script uses the dndscv package to perform dN/dS analysis on gene-level focusing on lung cancer driver for clonal and subclonal mutations, a test for selective pressure in the WES data (Figure S14).

**Human_WES_05_Summary_plot.R** - Script to generate summary plot of preinvasive samples and patient characteristics used for assessment of clonality
