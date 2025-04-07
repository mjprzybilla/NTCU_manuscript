# FIGURE 5 - LCM-WGS

## Data

**sample_signatures_Lung_2023_08_20.txt** - Fitted SBS signature file per clone in NTCU mouse data.

**trinuc_mut_mat_full.txt** - SBS counts per clone in NTCU mouse data.

**ALL_clone_FITTED_DBS_Signature_countribution_matrix.txt** - Fitted DBS signature file per clone in NTCU mouse data.

**ALL_clone_FITTED_DBS_count_matrix.txt** - DBS counts per clone in NTCU mouse data.

**ALL_clones_FITTED_INDEL_Signature_countribution_matrix.txt** - Fitted ID signature file per clone in NTCU mouse data.

**ALL_clones_FITTED_INDEL_count_matrix.txt** - Indel counts per clone in NTCU mouse data.

**NTCU_clonalityEstimate_table.txt** - Clonality estimations per LCM microbiopsy across each NTCU-treated mouse. 

## Scripts

**Mouse_WGS_create_multiChain_hdp_object_NTCU_clones.R** - Script to extract SBS mutational signatures using the NTCU mouse data.

**Mouse_WGS_snv_dnv_burden.R** - Script to visualise extracted SBS and ID signatures per clone from NTCU mouse data. 

**Mouse_WGS_ndp_tree_generation_NTCU.R** - Script to generate and visualise phylogenetic trees from NDP-evaluated clones in NTCU mouse data. The respective phylogenetic tree objects are given in the folder for Figure 6.

**Mouse_WGS_create_trees_with_signatures.R** - Script to visualise phylogenetic trees from NDP-evaluated clones with mutational signatures in NTCU mouse data.

