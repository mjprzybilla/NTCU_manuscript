# FIGURE 2

This GitHub repository contains the following folders:
```
└─ Panels e and f: directory containing scripts to reproduce panels e and f (6 files, 3 data files and 3 scripts)
|   └─ data_voidarea_cells (*.mat and *.xlsx): both contain the same data, corresponding to the area of voids (units of cells) obtained for the 3 conditions: control, short term and long term NTCU treament, measured at 24 weeks.
|   └─ data_votervoid_size: corresponding to the area of voids (units of cells) obtained from 50 simulations of the neutral voter model.
|   └─ plot_Figure_2e_f: loads the void size data and plots the cumulative distribution of void sizes.
|   └─ script_voter_model: corresponds to the stochastic simulation of the voter model
|   └─ run_script_voter_model: compiles script_voter_model into a Matlab *.mex file and runs a single realization of the model, it allows to plot a timelapse of the dynamics and the distribution of void sizes.
