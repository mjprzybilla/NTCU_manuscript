# FIGURE 2

This GitHub repository contains the files to reproduce Figures 2E,F and G.

Content:
```
└─ Plot_Figure_2EGF: load the raw data in data_voidarea_cells.mat, computes and plot the chefs shown in Fig 2E, F and G. Curves for the theoretical model can be obtained from running the scripts in Models/neutral_model (for 2F) and Models/non_neutral_model (for 2G).
└─ data
| └─ area_hole_cells.mat: area of voids (units of cells) obtained for the 3 conditions: two controls, and long term NTCU treatment, measured at 24 weeks (used for reproducing Figs. 2E,F and G).
| └─ neutral_sim_params.mat and neutral_data_holes.mat: simulation parameters and area of voids (units of cells) obtained from 200 repeats of the neutral model (used for reproducing Fig. 2F).
| └─ non_neutral_sim_params.mat and non_neutral_data_holes.mat: simulation parameters and area of voids (units of cells) obtained from 200 repeats of the non neutral model (used for reproducing Fig. 2G).
```