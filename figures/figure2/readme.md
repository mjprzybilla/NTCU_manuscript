# FIGURE 2

This GitHub repository contains the files to reproduce Figures 2E,F and G.

Content:
```
└─ Plot_Figure_2EGF.m: load the raw data in directory data, computes and plot the chefs shown in Fig 2E, F and G. 
└─ data/: Directory containing experimental and numerical data, required by Plot_Figure_2EFG.m
|   └─ area_hole_cells.mat: area of voids (units of cells) obtained for the 3 conditions: two controls, and long term NTCU treatment, measured at 24 weeks (used for reproducing Figs. 2E,F and G).
|   └─ neutral_sim_params.mat and neutral_data_holes.mat: simulation parameters and area of voids (units of cells) obtained from 200 repeats of the neutral model (used for reproducing Fig. 2F).
|   └─ non_neutral_sim_params.mat and non_neutral_data_holes.mat: simulation parameters and area of voids (units of cells) obtained from 200 repeats of the non neutral model (used for reproducing Fig. 2G).
```
