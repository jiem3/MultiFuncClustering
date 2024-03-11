The version of R packages used in the codes are:

mvtnorm: 1.1.3
xts: 0.13.1
stringr: 1.5.0
fpca: 0.2.1
MFPCA: 1.3.10
Rmpfr: 0.9.2
fda: 6.1.4
funData: 1.3.8
clues: 0.6.2.2
ggplot2: 3.4.1
Hmisc: 5.1.0
reshape2: 1.4.4
snowfall: 1.84-6.2


There are eight files and two folders under this directory, which are listed below.

Files:

README.md

funcs.r: contains functions used in this paper. 

simures_func.r: contains functions related to simulation results evaluation. 

script_simu.r: major R codes used for simulation in the paper.

script_sysB.r: major R codes used for System B data analysis in the paper.

sysBpca.rds: functional principal component scores of the System B data used in the paper.

fpca_0.2-1.tar.gz: source for installing the fpca package.

clues_0.6.2.2.tar.gz: source for installing the clues pcakge.

Folder:

pca3: contains EM results with 10 different starting points in System B data analysis.

simulation_code: contains codes used for simulation analysis. There are 2 more folders inside it, which are listed below.

nopenalty: contains simulation codes related to no variable selection clustering

penalty: contains simulation codes related to clustering with variable selection, using the three penalties mentioned in the paper.

There are 10 files in each of the two folders.

Simulation_N50.r,Simulation_N350.r, and Simulation_N500.r: codes for simulation with different sample sizes.
(The codes for setting N200 is the same as Simulation_noise16.r.)

Simulation_nois8.r, Simulation_noise16.r, Simulation_noise32.r, and Simulation_nois64.r are for simulation with different number of noisy sensors.

Simulation_var1.r, Simulation_var2.r, and Simulation_var2.5.r are for simulation with different signal strength.
(The codes for var1.5 is the same as Simulation_noise16.r)






 
