########################################################################################
#This file contains the major R code for the simulation and data analysis in
#paper "Multivariate Functional Clustering with Variable Selection and Application to Sensor Data from Engineering Systems"
#by Zhongnan Jin, Jie Min, Yili Hong, Pang Du, and Qingyu Yang
#########################################################################################

#setwd("E:/Dr.Hong/Zhongnanclustering/final_code/Rcode")
setwd("your directory path")


source("funcs.r")
source("simures_func.R")

################################################################################
#start data analysis
################################################################################
##Read the rescaled system B data
data.30.before.xts<-readRDS('sysBrescale.rds')
colnames(data.30.before.xts)[67]='Stop_ID'

##Figure1(a)
par(mai=c(0.85, 0.85, .05, .05))
VAR <- 64 
L <- 419
minY <- min(data.30.before.xts[,VAR])
maxY <- max(data.30.before.xts[,VAR])
plot(1:30,data.30.before.xts[data.30.before.xts[,67]==1,VAR],type='l', ylim=c(minY, maxY),
     xlab = "Measure Points", ylab = "Measurements", main="",cex.axis=1.5,cex.lab = 1.3)
for(i in 1:L){
  cur<-data.30.before.xts[data.30.before.xts[,67]==i,]
  lines(1:30,as.numeric(coredata(cur[,VAR],type="l")),
        col = i,lwd = 0.5
  )

}

##Figure1(b)
par(mai=c(0.85, 0.85, .05, .05))
VAR <- 32 
L <- 419 
minY <- min(data.30.before.xts[,VAR])
maxY <- max(data.30.before.xts[,VAR])
plot(1:30,data.30.before.xts[data.30.before.xts[,67]==1,VAR],type='l', ylim=c(minY, maxY),
     xlab = "Measure Points", ylab = "Measurements", main="",cex.axis=1.5, cex.lab = 1.3)
for(i in 1:L){
  cur<-data.30.before.xts[data.30.before.xts[,67]==i,]
  lines(1:30,as.numeric(coredata(cur[,VAR],type="l")),
        col = i,lwd = 0.5
  )
}

##Figure1(c)
par(mai=c(0.85, 0.85, .05, .05))
VAR <- 2 
L <- 419 
minY <- min(data.30.before.xts[,VAR])
maxY <- max(data.30.before.xts[,VAR])
plot(1:30,data.30.before.xts[data.30.before.xts[,67]==1,VAR],type='l', ylim=c(minY, maxY),
     xlab = "Measure Points", ylab = "Measurements", main="",cex.axis=1.5, cex.lab=1.3)
for(i in 1:L){
  cur<-data.30.before.xts[data.30.before.xts[,67]==i,]
  lines(1:30,as.numeric(coredata(cur[,VAR],type="l")),
        col = i,lwd = 0.5
  )
}

##Figure1(d)
par(mai=c(0.85, 0.85, .05, .05))
VAR <- 38 
L <- 419 
minY <- min(data.30.before.xts[,VAR])
maxY <- max(data.30.before.xts[,VAR])
plot(1:30,data.30.before.xts[data.30.before.xts[,67]==1,VAR],type='l', ylim=c(minY, maxY),
     xlab = "Measure Points", ylab = "Measurements", main="",cex.axis=1.5, cex.lab =1.3)
for(i in 1:L){
  cur<-data.30.before.xts[data.30.before.xts[,67]==i,]
  lines(1:30,as.numeric(coredata(cur[,VAR],type="l")),
        col = i,lwd = 0.5
  )
}


################################################################################
#################################SIMULATION#####################################
################################################################################

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

##It is better to run simulations on a server using parallel computing
##An example of running the simulation code in command line is:
#for((i=1; i<= 200; i++)) 
#do
#  R CMD BATCH "--args SEED=$i" Simulation_nois8.R nois8_${i}.Rout &
#  done
#wait

#################################Simulation Results#############################
##simulation for sample size n
##input the directory path for the stored simulation results
##path1, path2, path3, path4 are directories for n=500, n=350, n=200, n=50 with penalties
##nppath1, nppath2, nppath3, nppath4 are directories for n=500, n=350, n=200, n=50 without penalty

##Table 1
simu.collect.N(path1,path2,path3,path4,nppath1,nppath2,nppath3,nppath4)
##Figure 3
simu.collect.N(path1,path2,path3,path4,nppath1,nppath2,nppath3,nppath4,
               TABLE=FALSE,PLOT=TRUE)

##simulation for number of niosy sensors pn
##input the directory path for the stored simulation results
##path1, path2, path3, path4 are directories for pn=64, pn=32, pn=16, pn=8 with penalties
##nppath1, nppath2, nppath3, nppath4 are directories for pn=64, pn=32, pn=16, pn=8 without penalty

##Table 2
simu.collect.Noisy(path1,path2,path3,path4,nppath1,nppath2,nppath3,nppath4)
##Figure 4
simu.collect.Noisy(path1,path2,path3,path4,nppath1,nppath2,nppath3,nppath4,
                   TABLE=FALSE,PLOT=TRUE)

##simulation for signal strength
##input the directory path for the stored simulation results
##path1, path2, path3, path4 are directories for delta=1, delta=1.5, delta=2, delta=2.5 with penalties
##nppath1, nppath2, nppath3, nppath4 are directories for delta=1, delta=1.5, delta=2, delta=2.5  without penalty

##Table 3
simu.collect.SigS(path1,path2,path3,path4,nppath1,nppath2,nppath3,nppath4)
##Figure 5
simu.collect.SigS(path1,path2,path3,path4,nppath1,nppath2,nppath3,nppath4,
                  TABLE=FALSE,PLOT=TRUE)


################################################################################
#####################Application to System B Data###############################
################################################################################
##run EM algorithm with 10 different starting point and pick the best one use BIC
library(snowfall)
no.nodes=1
no.cores=10

sfInit(parallel=TRUE, cpus=no.cores,type="SOCK")
sfExportAll()
sfClusterSetupRNG(seed=round(2^32*runif(1)))
tmp=sfClusterApplyLB(1:10,real.onerun)
sfStop()
##After running the codes above, there are saved rds files for the three different
##penalties with 10 different starting point

EM_res_l_ind<-EM_res_l_var<-EM_res_l_group<-vector(mode="list",length=10)
BIC_indi<-BIC_var<-BIC_group<-rep(NA,10)

##read the saved files
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
EM_result_grouped$BIC  ##The group penalty with 4 clusters has the best BIC

##Calculate number of clusters using the three different penalties
table(EM_result_individual$Cluster)
table(EM_result_variable$Cluster)
table(EM_result_grouped$Cluster)

##Calculate number of censors removed using group penalty
rm.var.group<-EM_result_grouped$VAR_removed##The PC components that needs to be removed


rm.var.index<-c(29,30,31,32,49,50,97,98,99,100)#index from system B data
rm.var.index<-sort(c(rm.var.index,c(7,8,9,10,17,18,21,22,23,24,25,26,
                                    29,30,33,34,37,38,39,40,41,42,
                                    43,44,47,48,55,56,61,62,101,102,
                                    107,108,109,110,111,112)))#index from system B data
rm.sensor.index<-round(rm.var.index[seq(1,length(rm.var.index),by=2)+1]/2)

num_PCA<-3
remain.sensor.length<-66-length(rm.sensor.index)
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



##Figure 6:
##Read the rescaled system B data
data.30.before.xts<-readRDS('sysBrescale.rds')

par(oma = c(5,5,2,1.5), mfrow = c(4,4),mar=c(1,0, 0.5,1.5))
clusterlab <-c(2,4,3,1)
Varlist<-c(64,32,41,43)
sensorlist<-c(1,2,5,6)
for(z in 1:4){
  
  for(j in c(4,1,3,2)){

    
    VAR <- Varlist[z] 
    L <- which(EM_result_grouped$Cluster==j) 
    
    minY <- min(data.30.before.xts[,VAR])
    maxY <- max(data.30.before.xts[,VAR])+1
    plot(1:30,data.30.before.xts[data.30.before.xts[,67]==L[1],VAR],type='l', ylim=c(minY, maxY),
         xlab="", ylab="",axes=FALSE,col = j+1)
    
    
    
    axis(1,cex.axis=1.25, at= c(0,10,20,30),labels = c(0,10,20,30))
    
    
    if(j == 4 & z == 1){
      axis(2, at=c(-1,1,3),labels = c(-1,1,3),cex.axis=1.2)
    }
    if(j == 4 & z == 4){
      axis(2, at=c(-3,-1,1),labels = c(-3,-1,1),cex.axis=1.2)
    }
    if(j == 4 & z == 2){
      axis(2, at=c(-3,0,3),labels = c(-3,0,3),cex.axis=1.2)
    }
    if(j == 4 & z == 3){
      axis(2, at=c(-6,-1,4),labels = c(-6,-1,4),cex.axis=1.2)
    }
    
    
    corners = par("usr") 
    par(xpd = TRUE) #Draw outside plot area

    
    if(z == 1){
      text(y = corners[4]-0.5, x = mean(corners[1:2]), paste("Cluster", clusterlab[j]),  cex = 1.35)
    }
    if(j == 2){
      corners = par("usr") 
      par(xpd = TRUE) #Draw outside plot area
      text(x = corners[2]+2.5, y = mean(corners[3:4]), paste("Sensor", sensorlist[z]), srt = -90, cex = 1.35)
      
    }
    
    if(j == 4 & z == 2){
      mtext(text = "Measurements", side=2, line = 3, adj = -5, cex=1.1)
    }
    if(z == 4 & j == 1){
      mtext(text = "Measure Points",side=1, adj = 5,line = 3, cex=1.1)
    }
    
    tmp<-data.30.before.xts[data.30.before.xts[,67]==L[1],VAR]
    for(i in 2:length(L)){
      
      cur<-data.30.before.xts[data.30.before.xts[,67]==L[i],VAR]
      tmp = tmp + cur
      lines(1:30,as.numeric(coredata(cur,type="l"))
            ,col =j+1,lwd=0.5
      )
      
    }
    
    lines(1:30, tmp/length(L),type='l', ylim=c(minY, maxY),
          col = 1,lwd=2)
   
  }

}

##Figure 7:
par(oma = c(5,5,2,1.5), mfrow = c(2,4),mar=c(1,0, 0.5,1.5))
Varlist<-c(38,66)
clusterlab <-c(2,4,3,1)
sensorlist<-c(4,7)
for(z in 1:2){
  for(j in c(4,1,3,2)){
    VAR <- Varlist[z] # which sensor/varible we want to plot
    L <- which(EM_result_grouped$Cluster==j) # how many observations 
    
    minY <- min(data.30.before.xts[,VAR])
    maxY <- max(data.30.before.xts[,VAR])
    plot(1:30,data.30.before.xts[data.30.before.xts[,67]==L[1],VAR],type='l', ylim=c(minY, maxY),
         main="",xlab="",ylab="",axes=FALSE,col=j+1)
    tmp<-data.30.before.xts[data.30.before.xts[,67]==L[1],VAR]
    for(i in 1:length(L)){
      cur<-data.30.before.xts[data.30.before.xts[,67]==L[i],VAR]
      tmp <- tmp + cur
      lines(1:30,as.numeric(coredata(cur,type="l"))
            ,col =j+1,lwd=0.5)
    }
    
    lines(1:30, tmp/length(L),type='l', 
          col = 1,lwd=2)
    
    axis(1,cex.axis=1.25, at= c(0,10,20,30),labels = c(0,10,20,30))
    if(j == 4 & z == 2){
      axis(2, at=c(-6,0,5),labels = c(-6,0,5),cex.axis=1.2)
    }
    
    if(j == 4 & z == 1){
      axis(2, at=c(-6,-1,4),labels = c(-6,-1,4),cex.axis=1.2)
    }
   
    corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
    par(xpd = TRUE) #Draw outside plot area
    
    if(z == 1){
      text(y = corners[4]+0.1, x = mean(corners[1:2]), paste("Cluster", clusterlab[j]),  cex = 1.35)
    }
    if(j == 2){
      corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
      par(xpd = TRUE) #Draw outside plot area
      text(x = corners[2]+2.5, y = mean(corners[3:4]), paste("Sensor", sensorlist[z]), srt = -90, cex = 1.35)
      
    }
    if(j == 4 & z == 2){
      mtext(text = "Measurements", side=2, line = 3, adj = -5, cex=1.1)
    }
    if(z == 2 & j == 1){
      mtext(text = "Measure Points",side=1, adj = 5,line = 3, cex=1.1)
    }
    
  }
}





