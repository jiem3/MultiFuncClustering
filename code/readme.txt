R program for simulation part of 
"Multivariate Functional Clustering with Variable Selection and Application to Sensor Data from Engineering Systems"
by Zhongnan Jin, Jie Min, Yili Hong, Pang Du, and Qingyu Yang

There are four files and one folder under this directory, which are listed below.

Files:

readme.txt 

funs.r: contains functions used in this paper. 

simures_func.r: contains functions related to simulation results evaluation. 

script_simu.r: major R codes used for simulation in the paper.

Folder:

simulation_code: contains codes used for simulation analysis. There are 2 more folders inside it, which are listed below.

nopenalty: contains simulation codes related to no variable selection clustering

penalty: contains simulation codes related to clustering with variable selection, using the three penalties mentioned in the paper.

There are 10 files in each of the two folders, with the same file names.

Simulation_N50.r,Simulation_N350.r, and Simulation_N500.r: codes for simulation with different sample sizes.
(The codes for setting N200 is the same as Simulation_noise16.r.)

Simulation_nois8.r, Simulation_noise16.r, Simulation_noise32.r, and Simulation_nois64.r are for simulation with different number of noisy sensors.

Simulation_var1.r, Simulation_var2.r, and Simulation_var2.5.r are for simulation with different signal strength.
(The codes for var1.5 is the same as Simulation_noise16.r)


 
