########################################################################################
#This file contains the major R code for the System B data analysis in
#paper "Multivariate Functional Clustering with Variable Selection and Application to Sensor Data from Engineering Systems"
#by Zhongnan Jin, Jie Min, Yili Hong, Pang Du, and Qingyu Yang
#########################################################################################
setwd("your directory path")

source("funcs.r")

##Read the FPCA scores, run EM algorithm with 10 different starting points, and pick the best one using BIC
##It takes a long time to finish the analysis, 
##so we suggest to run the codes on a server using parallel computing.
library(snowfall)

no.nodes=1
no.cores=10

sfInit(parallel=TRUE, cpus=no.cores,type="SOCK")
sfExportAll()
sfClusterSetupRNG(seed=round(2^32*runif(1)))
sfSource("funcs.R")

tmp=sfClusterApplyLB(1:10,real.onerun)
sfStop()

##After running the codes above, there are saved rds files for the three different
##penalties with 10 different starting points in the pca3 folder
EM_res_l_ind<-EM_res_l_var<-EM_res_l_group<-vector(mode="list",length=10)
BIC_indi<-BIC_var<-BIC_group<-rep(NA,10)

##Read the saved files
for(i in 1:10){
  EM_res_l_ind[[i]]<-readRDS(file=paste("pca3/",i,"EM_result_individual_PCA.rds",sep=""))
  BIC_indi[i] <- as.numeric(EM_res_l_ind[[i]]$BIC[1])
  EM_res_l_var[[i]]<-readRDS(file=paste("pca3/",i,"EM_result_variable_PCA.rds",sep=""))
  BIC_var[i] <- as.numeric(EM_res_l_var[[i]]$BIC)
  EM_res_l_group[[i]]<-readRDS(file=paste("pca3/",i,"EM_result_grouped_PCA.rds",sep=""))
  BIC_group[i] <- as.numeric(EM_res_l_group[[i]]$BIC)
  
}

idx1<-which.min(unlist(BIC_indi))
idx2<-which.min(unlist(BIC_var))
idx3<-which.min(unlist(BIC_group))

EM_result_individual <- EM_res_l_ind[[idx1]]
EM_result_variable <- EM_res_l_var[[idx2]]
EM_result_grouped <- EM_res_l_group[[idx3]]

EM_result_individual$BIC 
EM_result_variable$BIC 
EM_result_grouped$BIC  #The group penalty with 4 clusters has the best BIC

##Calculate number of clusters using the three different penalties
table(EM_result_individual$Cluster)
table(EM_result_variable$Cluster)
table(EM_result_grouped$Cluster)

##Calculate number of censors removed using group penalty
rm.var.group<-EM_result_grouped$VAR_removed##The PC components that needs to be removed

##Check the number of removed sensors
num_PCA<-3
remain.sensor.length<-42
sensor.id.start<-seq(1,remain.sensor.length*num_PCA,by=num_PCA)
sensor.id.end<-sensor.id.start+(num_PCA-1)

rm.sensor.group<-NULL
check.start<-sensor.id.start
check.end<- sensor.id.end

for(i in 1:length(check.start)){
  tmpseq<-check.start[i]:check.end[i]
  if(all(tmpseq %in% rm.var.group)){

    sensor.id<-check.end[i]/num_PCA
    rm.sensor.group<-c(rm.sensor.group, sensor.id)
  }

}

length(rm.sensor.group)

