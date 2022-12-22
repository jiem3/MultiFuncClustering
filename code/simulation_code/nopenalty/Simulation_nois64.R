####################################
# Simulation Study: change number of noisy sensors
####################################
  


source("../../funcs.R")


  SEED <- 0

##read in the command line arguments
##run with: R CMD BATCH "--args seed=0 reps=5"
args <- commandArgs(TRUE)
if(length(args) > 0) 
    for(i in 1:length(args)) 
        eval(parse(text=args[[i]]))


  cat("seed is ", SEED, "\n", sep="")

  set.seed(SEED)


  
  N <- 200 ##sample size
  num_eval <- 31 
  
 
  mu1 <- c(c(14,14,rep(14,9),0), c(11,20,22,rep(23,5),seq(23,18,length.out=4)) ,rep(rep(0,12),64)) 
  mu2 <- c(c(11.5,12.5,13,seq(13,10,length.out=4),seq(9,4,length.out=4),0),
           c(10,15,13,seq(13,19,length.out=3),19,seq(20,23,length.out=3),20,18),
           rep(rep(0,12),64))
  mu3 <- c(c(5,seq(3,6.5,length.out=2),seq(6.5,3,length.out=2),seq(3,6.5,length.out=2),seq(6.5,3,length.out=2),seq(4,0,length.out=3)),
           c(rep(18,12)),
           rep(rep(0,12),64))

  
  pgseq<-1.5
  
  
  sigma1 <- diag(c(c(20,10,8,rep(2.8,6),rep(2.8,2),c(0.2)), 
                     c(rep(16,11),1), rep(15,12*64) ) )/pgseq  

  sigma2 <-diag(c(c(15,15,rep(10,5),seq(5,0.1,length.out=5)),
                    c(rep(30,5),rep(20,6),1),rep(15,12*64) ))/pgseq

  sigma3 <- diag( c(c(15,rep(12,2),rep(12,2),rep(12,2),rep(12,2),seq(3,0.2,length.out=3)),
                      c(rep(8,12)),rep(15,12*64)) )/pgseq  
 
  MU <- list(mu1, mu2, mu3)
  SIGMA <- list(sigma1, sigma2, sigma3)
  
  
  Data_result <- FDATA_GENERATOR(N = N, K =3, P = 66, PC = 12, order= 3, evalRange = c(0,30), MU = MU, SIGMA = SIGMA, prob = c(1/3,1/3,1/3) )
  
  Data_rescale<-cbind(apply(Data_result$Data[,1:66],2,scale),Data_result$Data[,67])

  FPCA_expansion <- MFPCA_sensory(Data_rescale, M = 3) 
  
  Lambda_set<-c(0)*(N^(1/3))
  gamma_set<-c(0)


  hyper_individual <- Hyper_Parameter_Selection(Lambda_set, gamma_set,Data = FPCA_expansion, K  = c(1,2,3,4), 
                                                penalty_method = "individual_equal", km = 3, 
                                                min_percentage = 0.03, plot = FALSE)
  
  hyper_variable <- Hyper_Parameter_Selection(Lambda_set, gamma_set,Data = FPCA_expansion, K = c(1,2,3,4), 
                                              penalty_method = "variable", km = 3, 
                                              min_percentage = 0.03, plot = FALSE)
  
  
  hyper_grouped <- Hyper_Parameter_Selection(Lambda_set, gamma_set,Data = FPCA_expansion, K = c(1,2,3,4), 
                                             penalty_method = "grouped", km = 3, 
                                             min_percentage = 0.03, plot = FALSE)

  
  optimal_individual <- EM_adaptive(MU_tilda = hyper_individual$MU_tilda, FPCA_expansion, K=hyper_individual$K, Lambda= hyper_individual$lambda ,
                                     Lambda2 = 10, sigma_method = 'equal', method = "individual_equal",
                                     M=200,gamma = hyper_individual$gamma, km = 3,
                                     min_percentage = 0.05)
  
  
  
  optimal_variable<- EM_adaptive(MU_tilda = hyper_variable$MU_tilda, FPCA_expansion, K=hyper_variable$K, 
                                          Lambda= hyper_variable$lambda, Lambda2 = 10, sigma_method = 'equal', 
                                          method = "variable", M=200, gamma = hyper_variable$gamma, km = 3,
                                          min_percentage = 0.05)
  optimal_grouped <- EM_adaptive(MU_tilda = hyper_grouped$MU_tilda, FPCA_expansion, K=hyper_grouped$K, Lambda= hyper_grouped$lambda, Lambda2 = 10, 
                                          sigma_method = 'equal', method = "grouped", M=200,
                                          gamma = hyper_grouped$gamma, km = 3, min_percentage = 0.05 )
    

  Rand_individual <- adjustedRand(optimal_individual$Cluster, Data_result$Cluster)[1]
  Rand_variable <- adjustedRand(optimal_variable$Cluster, Data_result$Cluster)[1]
  Rand_grouped <- adjustedRand(optimal_grouped$Cluster, Data_result$Cluster)[1]
  
  ARI_individual <- adjustedRand(optimal_individual$Cluster, Data_result$Cluster)[2]
  ARI_variable <- adjustedRand(optimal_variable$Cluster, Data_result$Cluster)[2]
  ARI_grouped <- adjustedRand(optimal_grouped$Cluster, Data_result$Cluster)[2]
  
  Jac_individual <- adjustedRand(optimal_individual$Cluster, Data_result$Cluster)[5]
  Jac_variable <- adjustedRand(optimal_variable$Cluster, Data_result$Cluster)[5]
  Jac_grouped <- adjustedRand(optimal_grouped$Cluster, Data_result$Cluster)[5]
  
  Remove_individual <- length(optimal_individual$VAR_removed)
  Remove_variable <- length(optimal_variable$VAR_removed)
  Remove_grouped <- length(optimal_grouped$VAR_removed)
  
  
  
  Remove_individual_s <- paste(optimal_individual$VAR_removed, sep = ",", collapse = ',')
  Remove_variable_s <- paste(optimal_variable$VAR_removed, sep = ",", collapse = ",")
  Remove_group_s <- paste(optimal_grouped$VAR_removed, sep = ",", collapse = ",")
  
  if(Remove_individual_s==''){Remove_individual_s <- "Not Removed"}
  if(Remove_variable_s==''){Remove_variable_s <- "Not Removed"}
  if(Remove_group_s==''){Remove_group_s <- "Not Removed"}
  
  result_data <- data.frame(matrix(c(hyper_individual$lambda, hyper_individual$gamma, 
                                     hyper_variable$lambda, hyper_variable$gamma,
                                     hyper_grouped$lambda, hyper_grouped$gamma,
                                     ARI_individual, ARI_variable, ARI_grouped,
                                     Rand_individual, Rand_variable, Rand_grouped,
                                     Jac_individual, Jac_variable, Jac_grouped,
                                     Remove_individual, Remove_variable, Remove_grouped,
                                     as.character(Remove_individual_s), as.character(Remove_variable_s), as.character(Remove_group_s),
                                     hyper_individual$K, hyper_variable$K, hyper_grouped$K), nrow = 1))
  
  

  
  names(result_data) <- c("lambda_individual", "gamma_individual",
                          "lambda_variable", "gamma_variable",
                          "lambda_group", "gamma_group",
                          "ARI_individual", "ARI_variable", "ARI_group",
                          "Rand_individual", "Rand_variable", "Rand_group",
                          "Jac_individual", "Jac_variable", "Jac_group",
                          "Var_remove_individual", "Var_remove_vaiable", "Var_remove_group",
                          "Var_remove_individual_d", "Var_remove_vaiable_d", "Var_remove_group_d",
                          "K_individual", "K_variable", "K_grouped")
  
  

 
 

  print(paste("Iteration ",SEED," Finished",sep=""))
  write.csv(result_data, file = paste("result_MB_FDA_nois64_nopara",SEED,".csv", sep = ""))
  




 