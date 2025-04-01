# FIGURE 2

This GitHub repository contains the scripts to run the numerical simulations of the neutral model, the non-neutral model, and the non-neutral model with driver mutations.

Content:
```
└─ data_sim_initial_condition: Images obtained from the experimental Tam Ctrl (4 days), which are used as initial conditions for the numerical simulations. 
└─ neutral_model: scripts to run simulations of the neutral model.
| └─ sim_neutral_model.m: script with numerical implementation of the model
| └─ main_run_sim_neutral_model.m: script to run the simulations and analyse the outputs
└─ non_neutral_model: scripts to run simulations of the non neutral model.
| └─ sim_non_neutral_model.m: script with numerical implementation of the non neutral model
| └─ main_run_sim_non_neutral_model.m: script to run the simulations and analyse the outputs
└─ non_neutral_model_w_driver_mutation: scripts to run simulations of non neutral model in which one clone has an advantage over neighbouring clones.
| └─ sim_non_neutral_w_driver: script with numerical implementation of the model
| └─ main_run_simulations_non_neutral_w_drive: script to run the simulations and analyse the outputs
| └─ sim_output_350x100_n16_mutant_n1_w_driver_mutation.mat: output of a single simulation of the model. Used to reproduce Supp Fig 12E
| └─ Plot_Supp_Fig_12E.m: script loads sim_output_350x100_n16_mutant_n1_w_driver_mutation.mat and produces the sequence of panels shown in Supp Fig 12E