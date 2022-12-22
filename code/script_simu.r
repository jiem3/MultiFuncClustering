########################################################################################
#This file contains the major R code for the simulation in
#paper "Multivariate Functional Clustering with Variable Selection and Application to Sensor Data from Engineering Systems"
#by Zhongnan Jin, Jie Min, Yili Hong, Pang Du, and Qingyu Yang
#########################################################################################

setwd("your directory path")


source("funcs.r")
source("simures_func.R")



##Figure2(a)
figure2.plot(x='a')
##Figure2(b)
figure2.plot(x='b')
##Figure2(c)
figure2.plot(x='c')
##The data for drawing figure 2 is randomly generated.
##If run the code multiple times, the generated figures are slightly different


##Please use the code in folder 'simulation_code' for the simulation analysis###
##The files in 'simulation_code/penalty' are codes for clustering with penalties
##The files in 'simulation_code/nopenalty' are for clustering without penalty


##Simulation_N50.r -- Simulation_N500.r are for simulation with different sample size settings.
##(the codes for setting N200 is the same as Simulation_noise16.r)
##Simulation_nois8.r -- Simulation_nois64.r are for simulation with different number of noisy sensors
##Simulation_var1.r -- Simulation_var2.5.r are for simulation with different signal strength
##(the codes for var1.5 is the same as Simulation_noise16.r)

##It takes a long time to finish the analysis, 
##so we suggest to run the codes on a server using parallel computing.
##An example of running the simulation codes in command line is:
#for((i=1; i<= 200; i++)) 
#do
#  R CMD BATCH "--args SEED=$i" Simulation_nois8.R nois8_${i}.Rout &
#  done
#wait

#################################Simulation Results#############################
##Simulation for sample size n.
##Inputs of the function are the directory path of the stored simulation results.
##path1, path2, path3, path4 are directories for n=500, n=350, n=200, n=50 with penalties
##nppath1, nppath2, nppath3, nppath4 are directories for n=500, n=350, n=200, n=50 without penalty

##Table 1
simu.collect.N(path1,path2,path3,path4,nppath1,nppath2,nppath3,nppath4)
##Figure 3
simu.collect.N(path1,path2,path3,path4,nppath1,nppath2,nppath3,nppath4,
               TABLE=FALSE,PLOT=TRUE)

##simulation for number of niosy sensors pn
##Inputs of the function are the directory path of the stored simulation results.
##path1, path2, path3, path4 are directories for pn=64, pn=32, pn=16, pn=8 with penalties
##nppath1, nppath2, nppath3, nppath4 are directories for pn=64, pn=32, pn=16, pn=8 without penalty

##Table 2
simu.collect.Noisy(path1,path2,path3,path4,nppath1,nppath2,nppath3,nppath4)
##Figure 4
simu.collect.Noisy(path1,path2,path3,path4,nppath1,nppath2,nppath3,nppath4,
                   TABLE=FALSE,PLOT=TRUE)

##simulation for signal strength
##Inputs of the function are the directory path of the stored simulation results.
##path1, path2, path3, path4 are directories for delta=1, delta=1.5, delta=2, delta=2.5 with penalties
##nppath1, nppath2, nppath3, nppath4 are directories for delta=1, delta=1.5, delta=2, delta=2.5  without penalty

##Table 3
simu.collect.SigS(path1,path2,path3,path4,nppath1,nppath2,nppath3,nppath4)
##Figure 5
simu.collect.SigS(path1,path2,path3,path4,nppath1,nppath2,nppath3,nppath4,
                  TABLE=FALSE,PLOT=TRUE)

