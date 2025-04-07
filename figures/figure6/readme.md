# Figure 6

## Scripts

**mapscape_run.R** - R script to create interactive mapscape visualisation based on input files provided in each folder (e.g. MD7047_Right_Lung). 

**Mapscape_Generation.Rmd** - R markdown script to create interactive mapscape visualisation based on input files provided in each folder (e.g. MD7047_Right_Lung). 

**mapscape_generation.sh** - Bash script executed as given below in order to create the interactive mapscape visualisation based on input files provided in each folder (e.g. MD7047_Right_Lung). 

```
./mapscape_generation.sh \
MD7047_Right_Lung \
figure6/MD7047_Right_Lung/MD7047_img.png \
figure6/MD7047_Right_Lung/MD7047_loc_tbl.xlsx \
figure6/MD7047_Right_Lung/cluster_and_samples.csv \
figure6/MD7047_Right_Lung/MD7047_Right_Lung.edge_tbl.csv \
FALSE \
figure6/MD7047_Right_Lung/ \
figure6/MD7047_Right_Lung/MD7047_Right_Lung_ndp_assigned_muts_final.csv
```
