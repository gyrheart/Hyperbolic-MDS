# Hyperbolic-MDS
Performing non-metric multi-dimensional scaling in hyperbolic space

Saved_data: this folder contains all the saved .mat files required. The Lukk’s dataset was not included here. Before you run the codes for real data analysis, you need to go to /Saved_data/SaveLukkData.m, download the data from:
https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-62/E-MTAB-62.processed.2.zip
Then transform and save data using the commands in SaveLukkData.m

Hyperbolic_functions: this folder contains all the functions required to run the codes.

Hyperbolic_mds:  this folder contains the code that tests hyperbolic MDS with synthetic results (Fig_mds_syn.m) and the code that reproduces Figure 3 in the manuscript (Fig_mds_Lukk.m). Note that the codes are commented so that you can directly generate the figures, if you want to run the hyperbolic MDS functions, just uncomment the codes.
