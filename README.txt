2024/02

Miles Miller, Thomas Ng, 2024
This work is licensed under CC BY-NC-SA 4.0 

README
Requirements:
MATLAB R2020b or higher.

The scripts in this directory performs 3 functions:

1) Simulate BRAFi/MEKi drug penetration in tumor.

Files:
- Generate_BRAFi_MEKi_TUMOR_profiles.m

This file contains the parameters for instantiating the tumor and the pertinent parameters for each BRAFi/MEKi drug, and runs the simulation. Edit the relevant parameter for simulation. 

The output of the file is 'OUTPUT_DRUG_CONC.mat', which contains the drug concentration and occupancy information for the BRAFi and MEKi being considered in this study.

Dependent files:

- simulate_TUMOR.m, ODEs2solve_MM1.m, main_solveODE.m

Files for simulation. 

- redblue.m
Redblue color map

Expected runtime: 5-10 minutes depending on size of tumor/length of simulation.

2) Generate time vs. concentration/fraction occupancy profiles for BRAFi/MEKi pairs

- timeline_conc.m, timeline_frac.m
Files for generating figures.
- arrowPlot_CC.m, arrowPlot.m, ccplot.m
Files generating functions for making the line graphs

Run after loading 'OUTPUT_DRUG_CONC.mat'.

3) Provides output for the integrated model
- PROLIF_DANA.mat
In vitro BRAFi/MEKi data from Schulz et al (PMID: 36230853) for Malme3M cell line. Load this file prior to running optimizecost_prolif.m

- optimizecost_prolif.m
Loads the output from section 1 to generate cell proliferation estimations across simulated tumors. 
