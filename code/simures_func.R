library(ggplot2)
library(Hmisc)
library(reshape2)





simu.collect.N<-function(path1="../results/16_500",
                         path2="../results/16_350",
                         path3="../results/16_200",
                         path4="../results/16_50",
                         nppath1="../results/16_500_np",
                         nppath2="../results/16_350_np",
                         nppath3="../results/16_200_np",
                         nppath4="../results/16_50_np",
                         TABLE=TRUE,
                         PLOT=FALSE
                         ){
  
  
  
if(TABLE==TRUE){
  n<-200
  
  
  
  
  File_N500 <- list.files(path = path1, pattern = NULL,
                          all.files = FALSE,
                          full.names = TRUE, recursive = FALSE,
                          ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_N500  <- combine_simu_data_cluster(File_names = File_N500 )
  
  
  File_N350 <- list.files(path = path2, pattern = NULL,
                          all.files = FALSE,
                          full.names = TRUE, recursive = FALSE,
                          ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_N350  <- combine_simu_data_cluster(File_names = File_N350 )
  
  File_N200 <- list.files(path = path3, pattern = NULL,
                          all.files = FALSE,
                          full.names = TRUE, recursive = FALSE,
                          ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_N200  <- combine_simu_data_cluster(File_names = File_N200 )
  
  
  File_N50 <- list.files(path = path4, pattern = NULL,
                         all.files = FALSE,
                         full.names = TRUE, recursive = FALSE,
                         ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_N50  <- combine_simu_data_cluster(File_names = File_N50 )
  
  
  
  
  
  
  ##500######
  Table_ss <- Table_N500[1:n,]
  
  
  
  a1<-colMeans(Table_ss[,c(17,18,19)])
  a2<- sapply(Table_ss[,c(17,18,19)], sd)
  
  var.rm.mean.500<-colMeans(Table_ss[,c(17,18,19)])
  
  
  
  a1<-Table_ss[,23] # individual penalty
  a2<-Table_ss[,24] # variable penalty
  a3<-Table_ss[,25] # Number of sensors removed
  
  mae.500<-cbind( sum(abs(a1-3))/length(a1), sum(abs(a2-3))/length(a2) ,sum(abs(a3-3))/length(a3))

  
  
  ##Here: only when 3 PCA components from one sensor is removed, the sensor is removed
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  s1<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sensor.rm.mean.500 <- c(m1,m2,m3)

  
  
  #fasly remove
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<- mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  
  sensor.rm.f.mean.500<-c(m1,m2,m3)
  
  
  
  ##350#######
  
  Table_ss <- Table_N350[1:n,]
  
  
  
  a1<-colMeans(Table_ss[,c(17,18,19)])
  a2<- sapply(Table_ss[,c(17,18,19)], sd)
  
  var.rm.mean.350<-colMeans(Table_ss[,c(17,18,19)])
  
  
  
  
  a1<-Table_ss[,23] # individual penalty
  a2<-Table_ss[,24] # variable penalty
  a3<-Table_ss[,25] # Number of sensors removed
  
  mae.350<-cbind( sum(abs(a1-3))/length(a1), sum(abs(a2-3))/length(a2) ,sum(abs(a3-3))/length(a3))

  
  
  ##Here: only when 3 PCA components from one sensor is removed, the sensor is removed
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  s1<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sensor.rm.mean.350 <- c(m1,m2,m3)

  
  #fasly remove
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<- mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  
  sensor.rm.f.mean.350<-c(m1,m2,m3)
  
  
  ##200####
  
  Table_ss <- Table_N200[1:n,]
  
  
  
  a1<-colMeans(Table_ss[,c(17,18,19)])
  a2<- sapply(Table_ss[,c(17,18,19)], sd)
  
  
  var.rm.mean.200<-colMeans(Table_ss[,c(17,18,19)])

  
  
  
  a1<-Table_ss[,23] # individual penalty
  a2<-Table_ss[,24] # variable penalty
  a3<-Table_ss[,25] # Number of sensors removed
  
  mae.200<-cbind( sum(abs(a1-3))/length(a1), sum(abs(a2-3))/length(a2) ,sum(abs(a3-3))/length(a3))

  
  
  ##Here: only when 3 PCA components from one sensor is removed, the sensor is removed
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  
  
  sensor.rm.mean.200 <- c(m1,m2,m3)

  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<- mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  
  sensor.rm.f.mean.200<-c(m1,m2,m3)

  #50#######
  
  Table_ss <- Table_N50[1:n,]
  var.rm.mean.50<-colMeans(Table_ss[,c(17,18,19)])
  var.rm.sd.50<- sapply(Table_ss[,c(17,18,19)], sd)
  
  
  # To calculate MAE for K
  a1<-Table_ss[,23] # individual penalty
  a2<-Table_ss[,24] # variable penalty
  a3<-Table_ss[,25] # Number of sensors removed
  
  mae.50<-cbind( sum(abs(a1-3))/length(a1), sum(abs(a2-3))/length(a2) ,sum(abs(a3-3))/length(a3))
  colnames(mae.50)<-c("individual","variable","group")
  
  
  # Number of sensors removed
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  s1<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  s3<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sensor.rm.mean.50 <- c(m1,m2,m3)

  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  s1<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  s3<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  sensor.rm.f.mean.50<-c(m1,m2,m3)
  
  mat.500<-cbind(c(mae.500),c(var.rm.mean.500),c(sensor.rm.mean.500),
                 c(sensor.rm.f.mean.500))
  mat.350<-cbind(c(mae.350),c(var.rm.mean.350),c(sensor.rm.mean.350),
                 c(sensor.rm.f.mean.350))
  mat.200<-cbind(c(mae.200),c(var.rm.mean.200),c(sensor.rm.mean.200),
                 c(sensor.rm.f.mean.200))
  mat.50<-cbind(c(mae.50),c(var.rm.mean.50),c(sensor.rm.mean.50),
                c(sensor.rm.f.mean.50))
  
  mat<-rbind(mat.50,mat.200,mat.350,mat.500) #mat.500,mat.1000
  
  ##############################NO penalty#####################################

  
  
  File_N500 <- list.files(path = nppath1, pattern = NULL,
                          all.files = FALSE,
                          full.names = TRUE, recursive = FALSE,
                          ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_N500  <- combine_simu_data_cluster(File_names = File_N500 )
  
  File_N350 <- list.files(path = nppath2, pattern = NULL,
                          all.files = FALSE,
                          full.names = TRUE, recursive = FALSE,
                          ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_N350  <- combine_simu_data_cluster(File_names = File_N350 )
  
  
  File_N200 <- list.files(path = nppath3, pattern = NULL,
                          all.files = FALSE,
                          full.names = TRUE, recursive = FALSE,
                          ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_N200  <- combine_simu_data_cluster(File_names = File_N200 )
  
  
  File_N50 <- list.files(path = nppath4 , pattern = NULL,
                         all.files = FALSE,
                         full.names = TRUE, recursive = FALSE,
                         ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_N50  <- combine_simu_data_cluster(File_names = File_N50 )
  
  
  
  n<-200
  
  

  
  ##8
  Table_ss <- Table_N500[1:n,]
  
  
  
  a1<-colMeans(Table_ss[,c(17,18,19)])
  a2<- sapply(Table_ss[,c(17,18,19)], sd)
  
  var.rm.mean.500<-colMeans(Table_ss[,c(17,18,19)])

  
  
  a1<-Table_ss[,23] # individual penalty
  a2<-Table_ss[,24] # variable penalty
  a3<-Table_ss[,25] # Number of sensors removed
  
  mae.500<-cbind( sum(abs(a1-3))/length(a1), sum(abs(a2-3))/length(a2) ,sum(abs(a3-3))/length(a3))

  
  ##Here: only when 3 PCA components from one sensor is removed, the sensor is removed
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:12, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:12, sensor_val = 1:6, col_index = 22)$result_vec)
  
  s1<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 20)$result_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 21)$result_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sensor.rm.mean.500 <- c(m1,m2,m3)

  
  #fasly remove
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<- mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 22)$False_vec)
  
  
  sensor.rm.f.mean.500<-c(m1,m2,m3)
  
  ##350
  
  Table_ss <- Table_N350[1:n,]
  
  
  
  a1<-colMeans(Table_ss[,c(17,18,19)])
  a2<- sapply(Table_ss[,c(17,18,19)], sd)
  
  var.rm.mean.350<-colMeans(Table_ss[,c(17,18,19)])

  
  
  
  
  a1<-Table_ss[,23] # individual penalty
  a2<-Table_ss[,24] # variable penalty
  a3<-Table_ss[,25] # Number of sensors removed
  
  mae.350<-cbind( sum(abs(a1-3))/length(a1), sum(abs(a2-3))/length(a2) ,sum(abs(a3-3))/length(a3))

  
  
  
  ##Here: only when 3 PCA components from one sensor is removed, the sensor is removed
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  s1<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sensor.rm.mean.350 <- c(m1,m2,m3)

  
  #fasly remove
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<- mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  
  sensor.rm.f.mean.350<-c(m1,m2,m3)
  
  
  
  ##200
  
  Table_ss <- Table_N200[1:n,]
  
  
  
  a1<-colMeans(Table_ss[,c(17,18,19)])
  a2<- sapply(Table_ss[,c(17,18,19)], sd)
  
  
  var.rm.mean.200<-colMeans(Table_ss[,c(17,18,19)])
  #var.rm.sd.4<- sapply(Table_ss[,c(17,18,19)], sd)
  
  
  
  a1<-Table_ss[,23] # individual penalty
  a2<-Table_ss[,24] # variable penalty
  a3<-Table_ss[,25] # Number of sensors removed
  
  mae.200<-cbind( sum(abs(a1-3))/length(a1), sum(abs(a2-3))/length(a2) ,sum(abs(a3-3))/length(a3))

  
  
  ##Here: only when 3 PCA components from one sensor is removed, the sensor is removed
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:12, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:12, sensor_val = 1:6, col_index = 22)$result_vec)
  
  
  
  sensor.rm.mean.200 <- c(m1,m2,m3)

  
  
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<- mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 22)$False_vec)
  
  
  sensor.rm.f.mean.200<-c(m1,m2,m3)

  
  # Number of variable removed
  
  Table_ss <- Table_N50[1:n,]
  var.rm.mean.50<-colMeans(Table_ss[,c(17,18,19)])
  var.rm.sd.50<- sapply(Table_ss[,c(17,18,19)], sd)
  
  
  # To calculate MAE for K
  a1<-Table_ss[,23] # individual penalty
  a2<-Table_ss[,24] # variable penalty
  a3<-Table_ss[,25] # Number of sensors removed
  
  mae.50<-cbind( sum(abs(a1-3))/length(a1), sum(abs(a2-3))/length(a2) ,sum(abs(a3-3))/length(a3))
  colnames(mae.50)<-c("individual","variable","group")
  
  
  # Number of sensors removed
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 22)$result_vec)
  
  s1<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 20)$result_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 21)$result_vec)
  s3<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sensor.rm.mean.50 <- c(m1,m2,m3)

  # Number of falsely removed sensors
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 22)$False_vec)
  
  s1<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 20)$False_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 21)$False_vec)
  s3<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:12, sensor_val = 1:6, col_index = 22)$False_vec)
  
  sensor.rm.f.mean.50<-c(m1,m2,m3)

  mat.no.500<-cbind(c(mae.500),c(var.rm.mean.500),c(sensor.rm.mean.500),
                    c(sensor.rm.f.mean.500))
  mat.no.350<-cbind(c(mae.350),c(var.rm.mean.350),c(sensor.rm.mean.350),
                    c(sensor.rm.f.mean.350))
  mat.no.200<-cbind(c(mae.200),c(var.rm.mean.200),c(sensor.rm.mean.200),
                    c(sensor.rm.f.mean.200))
  mat.no.50<-cbind(c(mae.50),c(var.rm.mean.50),c(sensor.rm.mean.50),
                   c(sensor.rm.f.mean.50))
  
  mat.no<-rbind(mat.no.50,mat.no.200,mat.no.350,mat.no.500)
  #############################################################################
  
  mat.res<-cbind(mat,mat.no[,1])
  
  ###############################################################################
  
  M=round(mat.res,2)
  m=nrow(M)
  n=ncol(M)
  cat("\\begin{pmatrix}",fill=T)
  for(i in 1:m)
  {
    
    if(i==1) {cat("\\multirow{3}{*}{50} &",fill=T)}else if(i==4){
      cat("\\multirow{3}{*}{200} &",fill=T)} else if (i==7) {
        cat("\\multirow{3}{*}{350} &",fill=T)
      }else if(i==10) {cat("\\multirow{3}{*}{500} &",fill=T)}else{cat("&",fill=T)}
    
    if(i %in% c(1,4,7,10)) cat("Individual &",fill=T)
    if(i %in% c(2,5,8,11)) cat("Variable &",fill=T)
    if(i %in% c(3,6,9,12)) cat("Group &",fill=T)
    
    
    tt=M[i,1]
    if (n>1)
    {
      for(j in 2:n)
      {
        tt=paste(tt,"&",M[i,j])
      }
    }
    cat(tt,"\\\\",fill=T)
    
    if(i==3) cat("\\hline\\hline",fill=T)
    if(i==6) cat("\\hline\\hline",fill=T)
    if(i==9) cat("\\hline\\hline",fill=T)
  }
  
  cat("\\end{pmatrix}",fill=T)
}
  

if(PLOT==TRUE){
  

  File_names_ss_50 <- list.files(path = path4, pattern = NULL,
                                 all.files = FALSE,
                                 full.names = TRUE, recursive = FALSE,
                                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_ss_50 <- combine_simu_data_cluster(File_names = File_names_ss_50)
  
  

  File_names_ss_200 <- list.files(path = path2, pattern = NULL,
                                  all.files = FALSE,
                                  full.names = TRUE, recursive = FALSE,
                                  ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_ss_200 <- combine_simu_data_cluster(File_names = File_names_ss_200)
  
  File_names_ss_350 <- list.files(path = path3, pattern = NULL,
                                  all.files = FALSE,
                                  full.names = TRUE, recursive = FALSE,
                                  ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_ss_350 <- combine_simu_data_cluster(File_names = File_names_ss_350)
  

  File_names_ss_500 <- list.files(path = path1, pattern = NULL,
                                  all.files = FALSE,
                                  full.names = TRUE, recursive = FALSE,
                                  ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_ss_500 <- combine_simu_data_cluster(File_names = File_names_ss_500)
  

  
  File_N500 <- list.files(path = nppath1, pattern = NULL,
                          all.files = FALSE,
                          full.names = TRUE, recursive = FALSE,
                          ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_ss_500.no  <- combine_simu_data_cluster(File_names = File_N500 )
  
  
  File_N200 <- list.files(path = nppath3, pattern = NULL,
                          all.files = FALSE,
                          full.names = TRUE, recursive = FALSE,
                          ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_ss_200.no  <- combine_simu_data_cluster(File_names = File_N200 )
  
  File_names_ss_350 <- list.files(path = nppath2, pattern = NULL,
                                  all.files = FALSE,
                                  full.names = TRUE, recursive = FALSE,
                                  ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_ss_350.no <- combine_simu_data_cluster(File_names = File_names_ss_350)
  
  
  
  File_N50 <- list.files(path = nppath4, pattern = NULL,
                         all.files = FALSE,
                         full.names = TRUE, recursive = FALSE,
                         ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_ss_50.no  <- combine_simu_data_cluster(File_names = File_N50 )
  
  

  
  #########################
  n<-200
  ARI_val_50_mat <- cbind(Table_ss_50[1:n,c(8,9,10)],Table_ss_50.no[1:n,9])
  ARI_val_50_mat <- cbind.data.frame(ARI_val_50_mat, rep("50", dim(ARI_val_50_mat)[1]) )
  ARI_val_200_mat <- cbind(Table_ss_200[1:n,c(8,9,10)],Table_ss_200.no[1:n,9])
  ARI_val_200_mat <- cbind.data.frame(ARI_val_200_mat, rep("200", dim(ARI_val_200_mat)[1]) )
  ARI_val_350_mat <- cbind(Table_ss_350[1:n,c(8,9,10)],Table_ss_350.no[1:n,9])
  ARI_val_350_mat <- cbind.data.frame(ARI_val_350_mat, rep("350", dim(ARI_val_350_mat)[1]) )
  
  ARI_val_500_mat <- cbind(Table_ss_500[1:n,c(8,9,10)],Table_ss_500.no[1:n,9])
  ARI_val_500_mat <- cbind.data.frame(ARI_val_500_mat, rep("500", dim(ARI_val_500_mat)[1]) )
  
  colnames(ARI_val_50_mat) <- c("Individual Penalty","Variable Penalty", "Group Penalty","No Penalty","label")
  colnames(ARI_val_200_mat) <- c("Individual Penalty","Variable Penalty", "Group Penalty","No Penalty","label")
  colnames(ARI_val_350_mat) <- c("Individual Penalty","Variable Penalty", "Group Penalty","No Penalty","label")
  colnames(ARI_val_500_mat) <- c("Individual Penalty","Variable Penalty", "Group Penalty","No Penalty","label")
  
  ####
  
  ARI_val_mat_ss<-rbind.data.frame(ARI_val_50_mat,ARI_val_200_mat,ARI_val_350_mat,ARI_val_500_mat)
  ################################

  df2_long <- melt(ARI_val_mat_ss, id.vars=c("label"))
  df2_long$ngroup <- factor(df2_long$label,                
                            levels = c("50", "200","350","500")) 
  
  p3 <- ggplot(df2_long, aes(x = "", y=value, fill=ngroup )) +
    geom_boxplot(lwd=0.2,outlier.size = 0.1) + 
    labs(title=" ",x = "", y = "ARI", fill = "Sample Size") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.margin = margin(0,0,0,0),legend.box.margin = margin(-30,-30,0,0),legend.position = "bottom") + ylim(0,1)+
    facet_wrap(~variable)  #+scale_fill_grey(start=1, end=0.3) 
  p3
  
  
  
  
}
  
}


simu.collect.Noisy<-function(path1="64_200",
                         path2="32_200",
                         path3="16_200",
                         path4="8_200",
                         nppath1="../pgseq1.5_nopen/64_200",
                         nppath2="../pgseq1.5_nopen/32_200",
                         nppath3="../pgseq1.5_nopen/16_200",
                         nppath4="../pgseq1.5_nopen/8_200",
                         TABLE=TRUE,
                         PLOT=FALSE
){
  
if(TABLE==TRUE){
  File_names_ss_64 <- list.files(path = path1, pattern = NULL,
                                 all.files = FALSE,
                                 full.names = TRUE, recursive = FALSE,
                                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_ss_64 <- combine_simu_data_cluster(File_names = File_names_ss_64)
  
  
  File_names_ss_32 <- list.files(path = path2, pattern = NULL,
                                 all.files = FALSE,
                                 full.names = TRUE, recursive = FALSE,
                                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_ss_32 <- combine_simu_data_cluster(File_names = File_names_ss_32)
  
  
  
  File_names_ss_16 <- list.files(path = path3, pattern = NULL,
                                 all.files = FALSE,
                                 full.names = TRUE, recursive = FALSE,
                                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_ss_16 <- combine_simu_data_cluster(File_names = File_names_ss_16)
  
  
  
  
  File_names_ss_8 <- list.files(path = path4, pattern = NULL,
                                all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_ss_8 <- combine_simu_data_cluster(File_names = File_names_ss_8)
  
  
  
  n<-200
  
  Table_ss <- Table_ss_64[1:n,]
  
  
  var.rm.mean.64<-colMeans(Table_ss[,c(17,18,19)])
  var.rm.sd.64<- sapply(Table_ss[,c(17,18,19)], sd)
  
  
  
  a1<-Table_ss[,23] # individual penalty
  a2<-Table_ss[,24] # variable penalty
  a3<-Table_ss[,25] # Number of sensors removed
  
  mae.64<-cbind( sum(abs(a1-3))/nrow(Table_ss), sum(abs(a2-3))/nrow(Table_ss) ,sum(abs(a3-3))/nrow(Table_ss))
  colnames(mae.64)<-c("individual","variable","group")
  
  
  
  ##Here: only when 3 PCA components from one sensor is removed, the sensor is removed
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:198, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:198, sensor_val = 1:6, col_index = 22)$result_vec)
  
  s1<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 20)$result_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 21)$result_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 22)$result_vec)
  
  
  sensor.rm.mean.64 <- c(m1,m2,m3)
  sensor.rm.sd.64 <-c(s1,s2,s3)
  
  
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<- mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 22)$False_vec)
  
  s1<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 20)$False_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 21)$False_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 22)$False_vec)
  
  
  sensor.rm.f.mean.64<-c(m1,m2,m3)
  sensor.rm.f.sd.64 <- c(s1,s2,s3)
  
  mat.64<-cbind(c(mae.64),c(var.rm.mean.64),c(var.rm.sd.64),c(sensor.rm.mean.64),
                c(sensor.rm.sd.64),c(sensor.rm.f.mean.64))
  
  
  
  ####32
  Table_ss <- Table_ss_32[1:n,]
  
  
  
  a1<-colMeans(Table_ss[,c(17,18,19)])
  a2<- sapply(Table_ss[,c(17,18,19)], sd)
  
  var.rm.mean.32<-colMeans(Table_ss[,c(17,18,19)])
  var.rm.sd.32<- sapply(Table_ss[,c(17,18,19)], sd)
  
  
  
  
  a1<-Table_ss[,23] # individual penalty
  a2<-Table_ss[,24] # variable penalty
  a3<-Table_ss[,25] # Number of sensors removed
  
  mae.32<-cbind( sum(abs(a1-3))/nrow(Table_ss), sum(abs(a2-3))/nrow(Table_ss) ,sum(abs(a3-3))/nrow(Table_ss))
  colnames(mae.32)<-c("individual","variable","group")
  
  
  
  ##Here: only when 3 PCA components from one sensor is removed, the sensor is removed
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:102, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:102, sensor_val = 1:6, col_index = 22)$result_vec)
  
  s1<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 20)$result_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 21)$result_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sensor.rm.mean.32 <- c(m1,m2,m3)
  sensor.rm.sd.32 <-c(s1,s2,s3)
  
  
  
  #fasly remove
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<- mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 22)$False_vec)
  
  s1<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 20)$False_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 21)$False_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 22)$False_vec)
  
  sensor.rm.f.mean.32<-c(m1,m2,m3)
  sensor.rm.f.sd.32 <- c(s1,s2,s3)
  
  
  
  
  
  
  ####16
  Table_ss <- Table_ss_16[1:n,]
  
  
  
  a1<-colMeans(Table_ss[,c(17,18,19)])
  a2<- sapply(Table_ss[,c(17,18,19)], sd)
  
  var.rm.mean.16<-colMeans(Table_ss[,c(17,18,19)])
  var.rm.sd.16<- sapply(Table_ss[,c(17,18,19)], sd)
  
  
  
  
  a1<-Table_ss[,23] # individual penalty
  a2<-Table_ss[,24] # variable penalty
  a3<-Table_ss[,25] # Number of sensors removed
  
  mae.16<-cbind( sum(abs(a1-3))/nrow(Table_ss), sum(abs(a2-3))/nrow(Table_ss) ,sum(abs(a3-3))/nrow(Table_ss))
  colnames(mae.16)<-c("individual","variable","group")
  
  
  
  ##Here: only when 3 PCA components from one sensor is removed, the sensor is removed
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  s1<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sensor.rm.mean.16 <- c(m1,m2,m3)
  sensor.rm.sd.16 <-c(s1,s2,s3)
  
  
  
  #fasly remove
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<- mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  s1<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  sensor.rm.f.mean.16<-c(m1,m2,m3)
  sensor.rm.f.sd.16 <- c(s1,s2,s3)
  
  
  
  ##8
  Table_ss <- Table_ss_8[1:n,]
  
  
  
  a1<-colMeans(Table_ss[,c(17,18,19)])
  a2<- sapply(Table_ss[,c(17,18,19)], sd)
  
  var.rm.mean.8<-colMeans(Table_ss[,c(17,18,19)])
  var.rm.sd.8<- sapply(Table_ss[,c(17,18,19)], sd)
  
  
  
  
  a1<-Table_ss[,23] # individual penalty
  a2<-Table_ss[,24] # variable penalty
  a3<-Table_ss[,25] # Number of sensors removed
  
  mae.8<-cbind( sum(abs(a1-3))/nrow(Table_ss), sum(abs(a2-3))/nrow(Table_ss) ,sum(abs(a3-3))/nrow(Table_ss))
  colnames(mae.8)<-c("individual","variable","group")
  
  
  
  ##Here: only when 3 PCA components from one sensor is removed, the sensor is removed
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:30, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:30, sensor_val = 1:6, col_index = 22)$result_vec)
  
  s1<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 20)$result_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 21)$result_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sensor.rm.mean.8 <- c(m1,m2,m3)
  sensor.rm.sd.8 <-c(s1,s2,s3)
  
  
  
  #fasly remove
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<- mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 22)$False_vec)
  
  s1<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 20)$False_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 21)$False_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 22)$False_vec)
  
  sensor.rm.f.mean.8<-c(m1,m2,m3)
  sensor.rm.f.sd.8 <- c(s1,s2,s3)
  
  
  
  mat.64<-cbind(c(mae.64),c(var.rm.mean.64),c(sensor.rm.mean.64),
                c(sensor.rm.f.mean.64))
  
  mat.32<-cbind(c(mae.32),c(var.rm.mean.32),c(sensor.rm.mean.32),
                c(sensor.rm.f.mean.32))
  
  mat.16<-cbind(c(mae.16),c(var.rm.mean.16),c(sensor.rm.mean.16),
                c(sensor.rm.f.mean.16))
  
  mat.8<-cbind(c(mae.8),c(var.rm.mean.8),c(sensor.rm.mean.8),
               c(sensor.rm.f.mean.8))
  # mat.4<-cbind(c(mae.4),c(var.rm.mean.4),c(sensor.rm.mean.4),
  #             c(sensor.rm.f.mean.4))
  # mat.2<-cbind(c(mae.2),c(var.rm.mean.2),c(sensor.rm.mean.2),
  #              c(sensor.rm.f.mean.2))
  
  # mat<-rbind(mat.2,mat.4,mat.8)#,mat.64
  mat<-rbind(mat.8,mat.16,mat.32,mat.64)
  ##################No penalty#############################
  
  File_names_ss_64 <- list.files(path = nppath1, pattern = NULL,
                                 all.files = FALSE,
                                 full.names = TRUE, recursive = FALSE,
                                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_ss_64 <- combine_simu_data_cluster(File_names = File_names_ss_64)
  
  
  File_names_ss_32 <- list.files(path = nppath2, pattern = NULL,
                                 all.files = FALSE,
                                 full.names = TRUE, recursive = FALSE,
                                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_ss_32 <- combine_simu_data_cluster(File_names = File_names_ss_32)
  

  
  File_names_ss_16 <- list.files(path = nppath3, pattern = NULL,
                                 all.files = FALSE,
                                 full.names = TRUE, recursive = FALSE,
                                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_ss_16 <- combine_simu_data_cluster(File_names = File_names_ss_16)
  
  
  
  
  
  
  File_names_ss_8 <- list.files(path = nppath4, pattern = NULL,
                                all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_ss_8 <- combine_simu_data_cluster(File_names = File_names_ss_8)
  
  
  
  
  
  n<-200
  
  Table_ss <- Table_ss_64[1:n,]
  
  
  var.rm.mean.64<-colMeans(Table_ss[,c(17,18,19)])
  var.rm.sd.64<- sapply(Table_ss[,c(17,18,19)], sd)
  
  
  
  a1<-Table_ss[,23] # individual penalty
  a2<-Table_ss[,24] # variable penalty
  a3<-Table_ss[,25] # Number of sensors removed
  
  mae.64<-cbind( sum(abs(a1-3))/nrow(Table_ss), sum(abs(a2-3))/nrow(Table_ss) ,sum(abs(a3-3))/nrow(Table_ss))
  colnames(mae.64)<-c("individual","variable","group")
  
  
  
  ##Here: only when 3 PCA components from one sensor is removed, the sensor is removed
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:198, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:198, sensor_val = 1:6, col_index = 22)$result_vec)
  
  s1<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 20)$result_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 21)$result_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 22)$result_vec)
  
  
  sensor.rm.mean.64 <- c(m1,m2,m3)
  sensor.rm.sd.64 <-c(s1,s2,s3)
  
  
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<- mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 22)$False_vec)
  
  s1<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 20)$False_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 21)$False_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:198, sensor_val = 1:6, col_index = 22)$False_vec)
  
  
  sensor.rm.f.mean.64<-c(m1,m2,m3)
  sensor.rm.f.sd.64 <- c(s1,s2,s3)
  
  mat.64<-cbind(c(mae.64),c(var.rm.mean.64),c(var.rm.sd.64),c(sensor.rm.mean.64),
                c(sensor.rm.sd.64),c(sensor.rm.f.mean.64))
  
  
  ##32
  
  
  
  ####32
  Table_ss <- Table_ss_32[1:n,]
  
  
  
  a1<-colMeans(Table_ss[,c(17,18,19)])
  a2<- sapply(Table_ss[,c(17,18,19)], sd)
  
  var.rm.mean.32<-colMeans(Table_ss[,c(17,18,19)])
  var.rm.sd.32<- sapply(Table_ss[,c(17,18,19)], sd)
  
  
  
  
  a1<-Table_ss[,23] # individual penalty
  a2<-Table_ss[,24] # variable penalty
  a3<-Table_ss[,25] # Number of sensors removed
  
  mae.32<-cbind( sum(abs(a1-3))/nrow(Table_ss), sum(abs(a2-3))/nrow(Table_ss) ,sum(abs(a3-3))/nrow(Table_ss))
  colnames(mae.32)<-c("individual","variable","group")
  
  
  
  ##Here: only when 3 PCA components from one sensor is removed, the sensor is removed
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:102, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:102, sensor_val = 1:6, col_index = 22)$result_vec)
  
  s1<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 20)$result_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 21)$result_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sensor.rm.mean.32 <- c(m1,m2,m3)
  sensor.rm.sd.32 <-c(s1,s2,s3)
  
  
  
  #fasly remove
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<- mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 22)$False_vec)
  
  s1<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 20)$False_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 21)$False_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:102, sensor_val = 1:6, col_index = 22)$False_vec)
  
  sensor.rm.f.mean.32<-c(m1,m2,m3)
  sensor.rm.f.sd.32 <- c(s1,s2,s3)
  
  
  
  
  
  ##16
  Table_ss <- Table_ss_16[1:n,]
  
  
  
  a1<-colMeans(Table_ss[,c(17,18,19)])
  a2<- sapply(Table_ss[,c(17,18,19)], sd)
  
  var.rm.mean.16<-colMeans(Table_ss[,c(17,18,19)])
  var.rm.sd.16<- sapply(Table_ss[,c(17,18,19)], sd)
  
  
  
  
  a1<-Table_ss[,23] # individual penalty
  a2<-Table_ss[,24] # variable penalty
  a3<-Table_ss[,25] # Number of sensors removed
  
  mae.16<-cbind( sum(abs(a1-3))/nrow(Table_ss), sum(abs(a2-3))/nrow(Table_ss) ,sum(abs(a3-3))/nrow(Table_ss))
  colnames(mae.16)<-c("individual","variable","group")
  
  
  
  ##Here: only when 3 PCA components from one sensor is removed, the sensor is removed
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  s1<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sensor.rm.mean.16 <- c(m1,m2,m3)
  sensor.rm.sd.16 <-c(s1,s2,s3)
  
  
  
  #fasly remove
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<- mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  s1<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  sensor.rm.f.mean.16<-c(m1,m2,m3)
  sensor.rm.f.sd.16 <- c(s1,s2,s3)
  
  
  ##8
  Table_ss <- Table_ss_8[1:n,]
  
  
  
  a1<-colMeans(Table_ss[,c(17,18,19)])
  a2<- sapply(Table_ss[,c(17,18,19)], sd)
  
  var.rm.mean.8<-colMeans(Table_ss[,c(17,18,19)])
  var.rm.sd.8<- sapply(Table_ss[,c(17,18,19)], sd)
  
  
  
  
  a1<-Table_ss[,23] # individual penalty
  a2<-Table_ss[,24] # variable penalty
  a3<-Table_ss[,25] # Number of sensors removed
  
  mae.8<-cbind( sum(abs(a1-3))/nrow(Table_ss), sum(abs(a2-3))/nrow(Table_ss) ,sum(abs(a3-3))/nrow(Table_ss))
  colnames(mae.8)<-c("individual","variable","group")
  
  
  
  ##Here: only when 3 PCA components from one sensor is removed, the sensor is removed
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:30, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table_ss, km = 3, true_val =  7:30, sensor_val = 1:6, col_index = 22)$result_vec)
  
  s1<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 20)$result_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 21)$result_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sensor.rm.mean.8 <- c(m1,m2,m3)
  sensor.rm.sd.8 <-c(s1,s2,s3)
  
  
  
  #fasly remove
  m1<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<- mean(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 22)$False_vec)
  
  s1<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 20)$False_vec)
  s2<-sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 21)$False_vec)
  s3<- sd(sensor_remove_check(Table_ss, km = 3, true_val = 7:30, sensor_val = 1:6, col_index = 22)$False_vec)
  
  sensor.rm.f.mean.8<-c(m1,m2,m3)
  sensor.rm.f.sd.8 <- c(s1,s2,s3)
  
 
  
  mat.no.64<-cbind(c(mae.64),c(var.rm.mean.64),c(sensor.rm.mean.64),
                   c(sensor.rm.f.mean.64))
  
  
  mat.no.32<-cbind(c(mae.32),c(var.rm.mean.32),c(sensor.rm.mean.32),
                   c(sensor.rm.f.mean.32))
  
  mat.no.16<-cbind(c(mae.16),c(var.rm.mean.16),c(sensor.rm.mean.16),
                   c(sensor.rm.f.mean.16))
  
  mat.no.8<-cbind(c(mae.8),c(var.rm.mean.8),c(sensor.rm.mean.8),
                  c(sensor.rm.f.mean.8))

  mat.no<-rbind(mat.no.8,mat.no.16,mat.no.32,mat.no.64)
  
  ##################################
  mat.res<-cbind(mat,mat.no[,1])
  ################################
  
  M=round(mat.res,2)
  m=nrow(M)
  n=ncol(M)
  cat("\\begin{pmatrix}",fill=T)
  for(i in 1:m)
  {
    
    if(i==1) {cat("\\multirow{3}{*}{8} &",fill=T)}else if(i==4){
      cat("\\multirow{3}{*}{16} &",fill=T)} else if (i==7) {
        cat("\\multirow{3}{*}{32} &",fill=T)
      }else if(i==10) {cat("\\multirow{3}{*}{64} &",fill=T)}else{cat("&",fill=T)}
    
    if(i %in% c(1,4,7,10)) cat("Individual &",fill=T)
    if(i %in% c(2,5,8,11)) cat("Variable &",fill=T)
    if(i %in% c(3,6,9,12)) cat("Group &",fill=T)
    
    
    tt=M[i,1]
    if (n>1)
    {
      for(j in 2:n)
      {
        tt=paste(tt,"&",M[i,j])
      }
    }
    cat(tt,"\\\\",fill=T)
    
    if(i==3) cat("\\hline\\hline",fill=T)
    if(i==6) cat("\\hline\\hline",fill=T)
    if(i==9) cat("\\hline\\hline",fill=T)
  }
  
  cat("\\end{pmatrix}",fill=T)
}
  
if(PLOT==TRUE){
  
  

  n<-200
  
  

  File_names_N8 <- list.files(path = path4, pattern = NULL,
                              all.files = FALSE,
                              full.names = TRUE, recursive = FALSE,
                              ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_N_8 <- combine_simu_data_cluster(File_names = File_names_N8)
  
  File_names_N16 <- list.files(path = path3, pattern = NULL,
                               all.files = FALSE,
                               full.names = TRUE, recursive = FALSE,
                               ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_N_16 <- combine_simu_data_cluster(File_names = File_names_N16)
  
  File_names_N32 <- list.files(path = path2, pattern = NULL,
                               all.files = FALSE,
                               full.names = TRUE, recursive = FALSE,
                               ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_N_32 <- combine_simu_data_cluster(File_names = File_names_N32)
  

  File_names_N64 <- list.files(path = path1, pattern = NULL,
                               all.files = FALSE,
                               full.names = TRUE, recursive = FALSE,
                               ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_N_64 <- combine_simu_data_cluster(File_names = File_names_N64)
  
  
  
  
  ####################No penal###############################
  File_names_ss_64 <- list.files(path = nppath1, pattern = NULL,
                                 all.files = FALSE,
                                 full.names = TRUE, recursive = FALSE,
                                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_N_64.no <- combine_simu_data_cluster(File_names = File_names_ss_64)
  
  
  
  
  File_names_N32 <- list.files(path = nppath2, pattern = NULL,
                               all.files = FALSE,
                               full.names = TRUE, recursive = FALSE,
                               ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_N_32.no <- combine_simu_data_cluster(File_names = File_names_N32)
  
  
  
  File_names_ss_16 <- list.files(path = nppath3, pattern = NULL,
                                 all.files = FALSE,
                                 full.names = TRUE, recursive = FALSE,
                                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_N_16.no <- combine_simu_data_cluster(File_names = File_names_ss_16)
  
  
  
  
  
  File_names_ss_8 <- list.files(path =nppath4, pattern = NULL,
                                all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_N_8.no <- combine_simu_data_cluster(File_names = File_names_ss_8)
  

  ARI_val8_mat <- cbind(Table_N_8[1:n,c(8,9,10)],Table_N_8.no[1:n,9])
  ARI_val8_mat <- cbind.data.frame(ARI_val8_mat, rep("8", dim(ARI_val8_mat)[1]) )
  ARI_val16_mat <- cbind(Table_N_16[1:n,c(8,9,10)],Table_N_16.no[1:n,9])
  ARI_val16_mat <- cbind.data.frame(ARI_val16_mat, rep("16", dim(ARI_val16_mat)[1]) )
  
  ARI_val32_mat <- cbind(Table_N_32[1:n,c(8,9,10)],Table_N_32.no[1:n,9])
  ARI_val32_mat <- cbind.data.frame(ARI_val32_mat, rep("32", dim(ARI_val32_mat)[1]) )
  ARI_val64_mat <- cbind(Table_N_64[1:n,c(8,9,10)],Table_N_64.no[1:n,9])
  ARI_val64_mat <- cbind.data.frame(ARI_val64_mat, rep("64", dim(ARI_val64_mat)[1]) )
  
  colnames(ARI_val8_mat) <- c("Individual Penalty","Variable Penalty", "Group Penalty","No Penalty","label")
  colnames(ARI_val16_mat) <- c("Individual Penalty","Variable Penalty", "Group Penalty","No Penalty","label")
  colnames(ARI_val32_mat) <- c("Individual Penalty","Variable Penalty", "Group Penalty","No Penalty","label")
  colnames(ARI_val64_mat) <- c("Individual Penalty","Variable Penalty", "Group Penalty","No Penalty","label")
  
 
  ARI_val_mat <- rbind.data.frame( ARI_val8_mat,ARI_val16_mat, ARI_val32_mat,ARI_val64_mat)
  
 
  df1_long <- melt(ARI_val_mat, id.vars=c("label"))
  df1_long$ngroup <- factor(df1_long$label,                 # Relevel group factor
                            levels = c("8", "16","32", "64"))#, "8","64"
  

  p2 <- ggplot(df1_long, aes(x="",y=value,fill=factor(ngroup)))+ 
    geom_boxplot(lwd=0.2,outlier.size=0.1) +facet_wrap(~variable)+ labs(title="",x = "", y = "ARI", fill = "Noisy Sensor") + 
    theme( plot.title = element_text(hjust = 0.5),legend.position = "bottom",
           legend.margin = margin(0,0,0,0),legend.box.margin = margin(-30,-30,0,0)) #+scale_fill_grey(start=1, end=0.3)
  p2

}
  
  
  
  
}




simu.collect.SigS<-function(path1="16_var1",
                             path2="16_200",
                             path3="16_var2",
                             path4="16_var2.5",
                             nppath1="../pgseq1.5_nopen/16_var1",
                             nppath2="../pgseq1.5_nopen/16_200",
                             nppath3="../pgseq1.5_nopen/16_var2",
                             nppath4="../pgseq1.5_nopen/16_var2.5",
                             TABLE=TRUE,
                             PLOT=FALSE
){
  
if(TABLE==TRUE){  
  
  n<-200
  
  
  File_names_var1 <- list.files(path = path1, pattern = NULL,
                                all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_var_1 <- combine_simu_data_cluster(File_names = File_names_var1)
  
  
  
  
  
  
  File_names_var2 <- list.files(path = path2, pattern = NULL,
                                all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_var_2 <- combine_simu_data_cluster(File_names = File_names_var2)
  
  
  
  
  File_names_var1 <- list.files(path = path3, pattern = NULL,
                                all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_var_25 <- combine_simu_data_cluster(File_names = File_names_var1)
  
  
  File_names_var2.5 <- list.files(path = path4, pattern = NULL,
                                  all.files = FALSE,
                                  full.names = TRUE, recursive = FALSE,
                                  ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_var_2.5 <- combine_simu_data_cluster(File_names = File_names_var2.5)
  
  
  
  #####var1##########
  Table <- Table_var_1[1:n,]
  
  
  
  # Number of variable removed
  var.rm.m.4<-colMeans(Table[,c(17,18,19)])
  
  
  # To calculate MAE for K
  a1<-Table[,23] # individual penalty
  a2<-Table[,24] # variable penalty
  a3<-Table[,25] # Number of sensors removed
  
  mae.4<-cbind( sum(abs(a1-3))/length(a1), sum(abs(a2-3))/length(a2) ,sum(abs(a3-3))/length(a3))
  colnames(mae.4)<-c("individual","variable","group")
  
  # Number of sensors removed
  m1<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6 ,col_index = 21)$result_vec)
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  
  sensor.rm.mean.4 <- c(m1,m2,m3)
  # Number of falsely removed sensors
  m1<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  sensor.rm.mean.f.4 <- c(m1,m2,m3)
  
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6,col_index = 22)$False_vec)
  
  ########var1.5#############
  
  Table <- Table_var_2[1:n,]
  
  
  
  
  var.rm.m.5<-colMeans(Table[,c(17,18,19)])
  
  
  # To calculate MAE for K
  a1<-Table[,23] # individual penalty
  a2<-Table[,24] # variable penalty
  a3<-Table[,25] # Number of sensors removed
  
  mae.5<-cbind( sum(abs(a1-3))/length(a1), sum(abs(a2-3))/length(a2) ,sum(abs(a3-3))/length(a3))
  colnames(mae.5)<-c("individual","variable","group")
  
  # Number of sensors removed
  m1<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6 ,col_index = 21)$result_vec)
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  
  sensor.rm.mean.5 <- c(m1,m2,m3)
  # Number of falsely removed sensors
  m1<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  sensor.rm.mean.f.5 <- c(m1,m2,m3)
  
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6,col_index = 22)$False_vec)
  
  
  
  
  #####var2##########
  Table <- Table_var_25[1:n,]
  
  
  # Number of variable removed
  var.rm.m.6<-colMeans(Table[,c(17,18,19)])
  
  
  # To calculate MAE for K
  a1<-Table[,23] # individual penalty
  a2<-Table[,24] # variable penalty
  a3<-Table[,25] # Number of sensors removed
  
  mae.6<-cbind( sum(abs(a1-3))/length(a1), sum(abs(a2-3))/length(a2) ,sum(abs(a3-3))/length(a3))
  colnames(mae.6)<-c("individual","variable","group")
  
  # Number of sensors removed
  m1<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6 ,col_index = 21)$result_vec)
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  
  sensor.rm.mean.6 <- c(m1,m2,m3)
  # Number of falsely removed sensors
  m1<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  sensor.rm.mean.f.6 <- c(m1,m2,m3)
  
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6,col_index = 22)$False_vec)
  
  
  
  
  
  
  ####var2.5#####
  Table <- Table_var_2.5[1:n,]
  
  
  
  # Number of variable removed
  var.rm.m.7<-colMeans(Table[,c(17,18,19)])
  
  
  # To calculate MAE for K
  a1<-Table[,23] # individual penalty
  a2<-Table[,24] # variable penalty
  a3<-Table[,25] # Number of sensors removed
  
  mae.7<-cbind( sum(abs(a1-3))/length(a1), sum(abs(a2-3))/length(a2) ,sum(abs(a3-3))/length(a3))
  colnames(mae.7)<-c("individual","variable","group")
  
  # Number of sensors removed
  m1<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6 ,col_index = 21)$result_vec)
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  
  sensor.rm.mean.7 <- c(m1,m2,m3)
  # Number of falsely removed sensors
  m1<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  sensor.rm.mean.f.7 <- c(m1,m2,m3)
  
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6,col_index = 22)$False_vec)
  
  
  
  mat.4<-cbind(c(mae.4),c(var.rm.m.4),c(sensor.rm.mean.4),
               c(sensor.rm.mean.f.4))
  mat.5<-cbind(c(mae.5),c(var.rm.m.5),c(sensor.rm.mean.5),
               c(sensor.rm.mean.f.5))
  mat.6<-cbind(c(mae.6),c(var.rm.m.6),c(sensor.rm.mean.6),
               c(sensor.rm.mean.f.6))
  mat.7<-cbind(c(mae.7),c(var.rm.m.7),c(sensor.rm.mean.7),
               c(sensor.rm.mean.f.7))
  
  mat<-rbind(mat.4,mat.5,mat.6,mat.7)
  ####################################No penalty########################
  
  
  
  
  n<-200
  
  
  File_names_var1 <- list.files(path = nppath1, pattern = NULL,
                                all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_var_1 <- combine_simu_data_cluster(File_names = File_names_var1)
  
  
  
  File_names_var2 <- list.files(path = nppath2, pattern = NULL,
                                all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_var_2 <- combine_simu_data_cluster(File_names = File_names_var2)
  
  
  
  
  File_names_var1 <- list.files(path = nppath3, pattern = NULL,
                                all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_var_25 <- combine_simu_data_cluster(File_names = File_names_var1)
  
  
  
  File_names_var2.5 <- list.files(path = nppath4, pattern = NULL,
                                  all.files = FALSE,
                                  full.names = TRUE, recursive = FALSE,
                                  ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_var_2.5 <- combine_simu_data_cluster(File_names = File_names_var2.5)
  
  
  ##var1###########
  Table <- Table_var_1[1:n,]
  
  
  
  # Number of variable removed
  var.rm.m.4<-colMeans(Table[,c(17,18,19)])
  
  
  # To calculate MAE for K
  a1<-Table[,23] # individual penalty
  a2<-Table[,24] # variable penalty
  a3<-Table[,25] # Number of sensors removed
  
  mae.4<-cbind( sum(abs(a1-3))/length(a1), sum(abs(a2-3))/length(a2) ,sum(abs(a3-3))/length(a3))
  colnames(mae.4)<-c("individual","variable","group")
  
  # Number of sensors removed
  m1<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6 ,col_index = 21)$result_vec)
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  
  sensor.rm.mean.4 <- c(m1,m2,m3)
  # Number of falsely removed sensors
  m1<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  sensor.rm.mean.f.4 <- c(m1,m2,m3)
  
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6,col_index = 22)$False_vec)
  
  
  ###var1.5#########
  Table <- Table_var_2[1:n,]
  
  
  
  # Number of variable removed
  var.rm.m.5<-colMeans(Table[,c(17,18,19)])
  
  
  # To calculate MAE for K
  a1<-Table[,23] # individual penalty
  a2<-Table[,24] # variable penalty
  a3<-Table[,25] # Number of sensors removed
  
  mae.5<-cbind( sum(abs(a1-3))/length(a1), sum(abs(a2-3))/length(a2) ,sum(abs(a3-3))/length(a3))
  colnames(mae.5)<-c("individual","variable","group")
  
  # Number of sensors removed
  m1<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6 ,col_index = 21)$result_vec)
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  
  sensor.rm.mean.5 <- c(m1,m2,m3)
  # Number of falsely removed sensors
  m1<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  sensor.rm.mean.f.5 <- c(m1,m2,m3)
  
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6,col_index = 22)$False_vec)
  
  
  
  
  ###############var2####
  Table <- Table_var_25[1:n,]
  
  
  # Number of variable removed
  var.rm.m.6<-colMeans(Table[,c(17,18,19)])
  
  
  # To calculate MAE for K
  a1<-Table[,23] # individual penalty
  a2<-Table[,24] # variable penalty
  a3<-Table[,25] # Number of sensors removed
  
  mae.6<-cbind( sum(abs(a1-3))/length(a1), sum(abs(a2-3))/length(a2) ,sum(abs(a3-3))/length(a3))
  colnames(mae.6)<-c("individual","variable","group")
  
  # Number of sensors removed
  m1<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6 ,col_index = 21)$result_vec)
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  
  sensor.rm.mean.6 <- c(m1,m2,m3)
  # Number of falsely removed sensors
  m1<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  sensor.rm.mean.f.6 <- c(m1,m2,m3)
  
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  sd(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6,col_index = 22)$False_vec)
  
  
  ####var2.5#####
  Table <- Table_var_2.5[1:n,]
  
  
  
  # Number of variable removed
  var.rm.m.7<-colMeans(Table[,c(17,18,19)])
  
  
  # To calculate MAE for K
  a1<-Table[,23] # individual penalty
  a2<-Table[,24] # variable penalty
  a3<-Table[,25] # Number of sensors removed
  
  mae.7<-cbind( sum(abs(a1-3))/length(a1), sum(abs(a2-3))/length(a2) ,sum(abs(a3-3))/length(a3))
  colnames(mae.7)<-c("individual","variable","group")
  
  # Number of sensors removed
  m1<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  m2<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$result_vec)
  m3<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6, col_index = 20)$result_vec)
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6 ,col_index = 21)$result_vec)
  sd(sensor_remove_check(Table, km = 2, true_val = 7:54, sensor_val = 1:6, col_index = 22)$result_vec)
  
  
  sensor.rm.mean.7 <- c(m1,m2,m3)
  # Number of falsely removed sensors
  m1<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 20)$False_vec)
  m2<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 21)$False_vec)
  m3<-mean(sensor_remove_check(Table, km = 3, true_val = 7:54, sensor_val = 1:6, col_index = 22)$False_vec)
  
  
  sensor.rm.mean.f.7 <- c(m1,m2,m3)
  
  
  ##########
  
  mat.4.no<-cbind(c(mae.4),c(var.rm.m.4),c(sensor.rm.mean.4),
                  c(sensor.rm.mean.f.4))
  mat.5.no<-cbind(c(mae.5),c(var.rm.m.5),c(sensor.rm.mean.5),
                  c(sensor.rm.mean.f.5))
  mat.6.no<-cbind(c(mae.6),c(var.rm.m.6),c(sensor.rm.mean.6),
                  c(sensor.rm.mean.f.6))
  mat.7.no<-cbind(c(mae.7),c(var.rm.m.7),c(sensor.rm.mean.7),
                  c(sensor.rm.mean.f.7))
  
  mat.no<-rbind(mat.4.no,mat.5.no,mat.6.no,mat.7.no)
  
  
  
  
  mat.res<-cbind(mat, mat.no[,1])
  
  
  ######################################
  
  M=round(mat.res,2)
  m=nrow(M)
  n=ncol(M)
  cat("\\begin{pmatrix}",fill=T)
  for(i in 1:m)
  {
    
    if(i==1) {cat("\\multirow{3}{*}{1} &",fill=T)}else if(i==4){
      cat("\\multirow{3}{*}{1.5} &",fill=T)} else if (i==7) {
        cat("\\multirow{3}{*}{2} &",fill=T)
      }else if(i==10) {cat("\\multirow{3}{*}{2.5} &",fill=T)}else{cat("&",fill=T)}
    
    if(i %in% c(1,4,7,10)) cat("Individual &",fill=T)
    if(i %in% c(2,5,8,11)) cat("Variable &",fill=T)
    if(i %in% c(3,6,9,12)) cat("Group &",fill=T)
    
    
    tt=M[i,1]
    if (n>1)
    {
      for(j in 2:n)
      {
        tt=paste(tt,"&",M[i,j])
      }
    }
    cat(tt,"\\\\",fill=T)
    
    if(i==3) cat("\\hline\\hline",fill=T)
    if(i==6) cat("\\hline\\hline",fill=T)
    if(i==9) cat("\\hline\\hline",fill=T)
  }
  
  cat("\\end{pmatrix}",fill=T)
}
  
if(PLOT==TRUE){
  
  
  
  File_names_var1 <- list.files(path = path1, pattern = NULL,
                                all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_var_1 <- combine_simu_data_cluster(File_names = File_names_var1)
  
  
  
  File_names_var2 <- list.files(path = path2, pattern = NULL,
                                all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_var_2 <- combine_simu_data_cluster(File_names = File_names_var2)
  
  
  File_names_var1 <- list.files(path = path3, pattern = NULL,
                                all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_var_25 <- combine_simu_data_cluster(File_names = File_names_var1)
  
  File_names_var2.5 <- list.files(path = path4, pattern = NULL,
                                  all.files = FALSE,
                                  full.names = TRUE, recursive = FALSE,
                                  ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_var_2.5 <- combine_simu_data_cluster(File_names = File_names_var2.5)
  
  
  
  ####No penl##############
  File_names_var1 <- list.files(path = nppath1, pattern = NULL,
                                all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_var_1.no <- combine_simu_data_cluster(File_names = File_names_var1)
  
  
  
  File_names_var2 <- list.files(path =nppath2, pattern = NULL,
                                all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_var_2.no <- combine_simu_data_cluster(File_names = File_names_var2)
  
  
  File_names_var1 <- list.files(path =nppath3, pattern = NULL,
                                all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_var_25.no <- combine_simu_data_cluster(File_names = File_names_var1)
  
  File_names_var2.5 <- list.files(path = nppath4, pattern = NULL,
                                  all.files = FALSE,
                                  full.names = TRUE, recursive = FALSE,
                                  ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Table_var_2.5.no <- combine_simu_data_cluster(File_names = File_names_var2.5)
  
  
  #############
  
  # 
  n<-200
  ARI_var1_mat <- cbind(Table_var_1[1:n,c(8,9,10)],Table_var_1.no[1:n,9])
  ARI_var1_mat <- cbind.data.frame(ARI_var1_mat, rep("1", dim(ARI_var1_mat)[1]) )
  ARI_var2_mat <- cbind(Table_var_2[1:n,c(8,9,10)],Table_var_2.no[1:n,9])
  ARI_var2_mat <- cbind.data.frame(ARI_var2_mat, rep("1.5", dim(ARI_var2_mat)[1]) )
  ARI_var25_mat <- cbind(Table_var_25[1:n,c(8,9,10)],Table_var_25.no[1:n,9])
  ARI_var25_mat <- cbind.data.frame(ARI_var25_mat, rep("2", dim(ARI_var25_mat)[1]) )
  
  ARI_var2.5_mat <- cbind(Table_var_2.5[1:n,c(8,9,10)],Table_var_2.5.no[1:n,9])
  ARI_var2.5_mat <- cbind.data.frame(ARI_var2.5_mat, rep("2.5", dim(ARI_var2.5_mat)[1]) )

  
  colnames(ARI_var1_mat) <- c("Individual Penalty","Variable Penalty", "Group Penalty","No Penalty","label")
  colnames(ARI_var2_mat) <- c("Individual Penalty","Variable Penalty", "Group Penalty","No Penalty","label")
  colnames(ARI_var25_mat) <- c("Individual Penalty","Variable Penalty", "Group Penalty","No Penalty","label")
  colnames(ARI_var2.5_mat) <- c("Individual Penalty","Variable Penalty", "Group Penalty","No Penalty","label")
  
  
  ARI_var_mat <- rbind.data.frame(ARI_var1_mat, ARI_var2_mat,ARI_var25_mat,ARI_var2.5_mat)
  
  
  # Extend to long data 
  df3_long <- melt(ARI_var_mat, id.vars=c("label"))
  df3_long$ngroup <- factor(df3_long$label,                 # Relevel group factor
                            levels = c("1", "1.5",  "2", "2.5"))
  
  
  p4 <- ggplot(df3_long, aes(x="",y=value,fill=factor(label)))+ 
    geom_boxplot(lwd=0.2,outlier.size = 0.1) + labs(title="",x = "", y = "ARI", fill = expression(~delta)) + 
    theme( plot.title = element_text(hjust = 0.5),
           legend.margin = margin(0,0,0,0),legend.box.margin = margin(-30,-30,0,0),legend.position = "bottom") + ylim(0, 1) + 
    facet_wrap(~variable)
  p4

 
  
}
  
  
  
}
