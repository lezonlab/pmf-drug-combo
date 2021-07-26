# pmf-drug-combo
Predicts effects of two-drug combinations using probabilistic matrix factorization

----------
File descriptions:

all60graph.m -- Generates plots for all 60 cell lines using NCI ALMANAC data
 
ComboDrugReader.java -- Reads raw NCI ALMANAC data and generates 60 files containing ComboScores on a per-cell line basis.

DrugEfficacyReader.java -- Reads raw NCI ALMANAC data and generates 60 files containing efficacies on a per-cell line basis.

HubBuilder.m -- Matlab function for generating random graphs using the hub method

pmfnest.m -- Matlab code for PMF

pmfSimulatedExperiment.m -- Matlab code for simulating PMF-guided experiments

README.md -- This file

ScaleFreeBuilder.m -- Matlab function for generating random scale-free graphs

SurfacePlotPMF.m -- Matlab code for generating Fig.2

WattsStrogatz.m -- Matlab code for generating random Watts-Strogatz graphs

----------
Usage: 
Step 1: Download NCI ALMANAC data as ComboDrugGrowth_Nov2017.csv
Step 2: In the command line, compile all java code with command "javac *.java"
Step 3: To run ComboDrugReader.java and create ComboScore matrices, use the command "java ComboDrugReader". Likewise, to run DrugEfficacyReader.java and create combination efficacy matrices use the command "java DrugEfficacyReader".
Step 4: For input into the PMF algorithm pmfnest.m, data must be normalized to mean-zero and unit variance.
Step 5: As in SurfacePlotPMF, reorganize data into a n^2 by 3 matrix, such that columns 1 and 2 contain the original row and column of each value respectively, and column 3 contains the old value.
Step 6: Run pmfnest with specified hyperparameters and batch size. Output matrices D and T are the matrix factor representations of PMF, and B = D*T' gives the PMF approximation of the original matrix.

