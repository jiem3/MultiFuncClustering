suppressPackageStartupMessages(library(mvtnorm))
suppressPackageStartupMessages(library(xts))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(fpca))
suppressPackageStartupMessages(library(MFPCA))
suppressPackageStartupMessages(library(Rmpfr))  


suppressPackageStartupMessages(library(fda))
suppressPackageStartupMessages(library(funData))
suppressPackageStartupMessages(library(clues))



FDATA_GENERATOR <- function(N, K, P, PC, order, evalRange, MU, SIGMA, prob)
{

  num_cluster <- K
  num_sensor <- P
  num_pc <- PC
  num_eval <- length(seq(evalRange[1], evalRange[2], 1))
  num_order <- order
  
  ##prob: proportion for mixed normals
  labels_cluster <- label_generator(N, K, prob) ## Generate labels using multinomial dist'n with K catagories for N obs's
  
  
  for(k in 1:num_cluster)
  {
    assign(paste('mu_', k, sep = ""),  MU[[k]] )
    assign(paste('sigma_', k, sep = ""), SIGMA[[k]]) ##Give names to MU??
  }
  
  # Create basis function
  ## generate spline basis, evalRange: time interval
  ## num_pc: number of spline basis 
  ## num_order: order of B-spline bsis
  basis_fn <- create.bspline.basis(rangeval = evalRange, nbasis = num_pc, norder = num_order)
 # plot(basis_fn)
  
  # Coefficient Matrix
  ## rows: obs
  ## cols: 3 coeff each spline basis( 3 spline bases for one sensor), 4 sensors in total
  coef_matrix <- matrix( rep(NA, N*num_sensor * num_pc), ncol = num_sensor * num_pc)
  
  # Data Matrix
  ## sensor+1 bc/ in num_eval: eg, 0-30 has 31 positions
  ##rows: 50 obs with 31 time evaluations per obs: 50x31
  ##cols: 4 sensors (4 cols bc/ adding up 3 spline basis together for each sensor )
  ## The last col is used to indicate the obs index (1:N)
  Eval_matrix <- matrix( rep(NA, N* num_eval*(num_sensor+1)), ncol = num_sensor+1)
  
  
  for(i in 1:N)
  {
    
    cluster <- labels_cluster[i] #Get label for i'th obs
    mu <- get(paste('mu_', cluster, sep = ""))
    sigma <- get(paste('sigma_', cluster, sep = ""))
    coef_matrix[i,] <- rmvnorm(1, mean = mu, sigma = sigma)
    
    
    # Create functional data with the basis for each sensor
    current_coef <- coef_matrix[i,] ##coeff for 3 basis of one sensor x 4 sensors=12 coeffs
    
    #browser()
    for(p in 1:num_sensor)
    {
      col_index <- ( (p-1)*num_pc +1 ) : ( (p-1)*num_pc + num_pc ) ##cols for the sepecific sensor
      row_index <- ( (i-1)*num_eval +1 ) : ( (i-1)*num_eval +num_eval ) ## rows for the evaluation of 0-30 ( 31 positions)
      
      coef_val <- current_coef[col_index]  #coeff for 4 sensors, 3 spline bases each
      fun_data <- fd(coef_val, basis_fn)
      
      # Obtain the distinct value
      fd_matrix <- eval.fd(seq(evalRange[1], evalRange[2], 1),fun_data) ##evaluate b-spline basis at the evaluation range and get the vector
      
      # Save functional evaluation in the matrix
      Eval_matrix[row_index,p] <- as.numeric( fd_matrix )
      
    }
    Eval_matrix[row_index,p+1] <- rep(i, num_eval)
    
  }
  
  return(list( Data = Eval_matrix, Coef = coef_matrix, Cluster = labels_cluster))
  
}



######################Data preperation#######################################

## Trim data on both tails 
# Calculating the variance after the trim
# data: vector of data 
# prob: vector of length 2 containing head and toe probabilities

trim_data <- function(data, prob=c(0.2,0.8))
{
  a <- data 
  result <- var( a[ which( a > quantile(a, prob, type =3 )[1] & a < quantile(a, prob, type =3)[2] ) ] )
  return(result)
}


# Calculate the Initial Value in E-M algorithm, assuming uniform Mu and Dirichlet of Alpha
# DATA: dataset with each row being an observation
# K: Number of groups users would like to generate
Initial_Values <- function(DATA, K)
{
  #browser()
  N <- dim(DATA)[1]
  p <- dim(DATA)[2]

  set.seed(Sys.time())
  kmeans_result <- kmeans(DATA, centers =K, iter.max=30)
  
  Mu <- kmeans_result$centers
  
  
  
  
  # calculate the initial proportions pi from Kmean results
  Alpha <- as.numeric(table(kmeans_result$cluster))/sum(as.numeric(table(kmeans_result$cluster)))

  Sigma <- matrix(rep(0,p*K), ncol = p)
  

  B <- rep(NA, K)
  C <- matrix(rep(NA, K*p), ncol = p )
  
  for(i in 1:K)
  {
    B[i] <- N*Alpha[i]/2
    
    for(j in 1:p)
    {
      C[i,j] <-  sum(Alpha[i]*(DATA[,j] - Mu[i,j])^2 )/2
      Sigma[i,j] <-  C[i,j]/B[i] 
    }
    
  }
  

  
  colnames(Mu) <- paste("Mu", 1:p, sep="")
  colnames(Sigma) <- paste("Sigma", 1:p, sep = "")
  
  return(list(Mu= Mu , Sigma =  Sigma, Alpha = Alpha) )
}

Initial_Values_old <- function(DATA, K)
{

  N <- dim(DATA)[1]
  p <- dim(DATA)[2]

  kmeans_result <- kmeans(DATA, centers = K, algorithm = "Lloyd", nstart = 500, iter.max=200)
  
  Mu <- kmeans_result$centers
  
  
  
  
  # calculate the initial proportions pi from Kmean results
  Alpha <- as.numeric(table(kmeans_result$cluster))/sum(as.numeric(table(kmeans_result$cluster)))

  Sigma <- matrix( apply(DATA,2,var) , ncol = p)

  
  colnames(Mu) <- paste("Mu", 1:p, sep="")
  colnames(Sigma) <- paste("Sigma", 1:p, sep = "")
  
  return(list(Mu= Mu , Sigma =  Sigma, Alpha = Alpha) )
}


PI_cal_old <- function(DATA, PI, sigma, mu, k)
{
  K <- length(mu)
  N <- dim(DATA)[1]
  DELTA <- data.frame()  # probability for each observation/proportion (ALL probability for one step of E-M)
  
  
  Pi <- matrix(rep(NA, N*K), ncol = K) # Create empty matrix to store delta
  
  for(i in 1:N)
  {
    
    X <- DATA[i,] 
    
    for(j in 1:K)
    {
      Mu_k <- mu[[j]][k,]
      Pi[i,j] <-  dmvnorm(x=X, mean=Mu_k, sigma = diag(sqrt(sigma[[j]][k,])),log = TRUE)
    }
    
    
    
    Pi[i,] <-  Pi[i,] - max(Pi[i,])  # add the same positive value to avoid 0
    Pi[i,] <- as.numeric( mpfr( exp(Pi[i,]), 1000 )  )# Transform probabilities back to regular scale rather than log

    Sum_Pi <-  sum(Pi[i,]*PI[k,])
    DELTA <- rbind(DELTA, Pi[i,]*PI[k,]/Sum_Pi) # Combine with the prior PI from last iteration
  }
  
  
  
  PI_update <-  colMeans(DELTA) # need to pass to the next iteration k+1
  
  return(list(DELTA = DELTA, Pi = PI_update))
}


PI_cal <- function(DATA, PI, sigma, mu, k)
{
  K <- length(mu)
  N <- dim(DATA)[1]
  DELTA <- data.frame()  # probability for each observation/proportion (ALL probability for one step of E-M)
  
  
  Pi <- matrix(rep(NA, N*K), ncol = K) # Create empty matrix to store delta
  
  for(i in 1:N)
  {
    
    X <- DATA[i,]  
    for(j in 1:K)
    {
      Mu_k <- mu[[j]][k,]
      Pi[i,j] <- dmvnorm_adj(X, Mu_k, sigma[[j]][k,])
    }
    

    if(min(Pi[i] < -745))
    {
      print(paste("less than 745",i))
      Pi[i,] <-  Pi[i,] - max(Pi[i,])  # add the same positive value to avoid 0
      Pi[i,] <- as.numeric( mpfr( exp(Pi[i,]), 1000 )  )# Transform probabilities back to regular scale rather than log
    }
    else
    {
      Pi[i,] <- as.numeric( mpfr( exp(Pi[i,]), 1000 )  )# Transform probabilities back to regular scale rather than log
    }

    
    # if all probabilities are 0, then assign probabilities equally. 
    if(all(Pi[i,]==0))
    {
      Pi[i,] <- 1/length(Pi[i,])
      print("changes Pi")
    }
    else{}

    Sum_Pi <- as.numeric( mpfr( sum(Pi[i,]*PI[k,]) , 1000 ) )
    DELTA <- rbind(DELTA, as.numeric( mpfr( Pi[i,]*PI[k,]/Sum_Pi, 1000)) )
  }
  
  
  
  PI_update <-  colMeans(DELTA) # need to pass to the next iteration k+1
  
  return(list(DELTA = DELTA, Pi = PI_update))
}

# Caculate the multivariate normal probability with un-correlated variance matrix
# We try to use more digit to avoid the exponential numerical issue. 
# X: one observation from a multivariate Gaussian distribution
# mu: the assumed mean vector for the Gassian distribution
# sigma: the assumed variance for the Gaussian distribution


dmvnorm_adj <- function(X, mu, sigma)
{
  K <- length(X)
  absV <- abs(prod(sigma))
  result <- - (absV)/2 - (1/2)*sum( (1/sigma)*((X-mu)^2) )

  return(result)
}


MU_cal_noEM <- function(DATA, DELTA, mu, sigma, Lambda, k, method = "individual", Adaptive = FALSE, gamma = 1, km = km)
{
  
    if(method == "individual_equal")
    {
      K <- length(mu)
      N <- dim(DATA)[1]
      p <- dim(DATA)[2]
      
      
      
      for(q in 1:K)
      {  
        for(r in 1:p)
        {
          if( Lambda < abs(sum(DATA[,r]))/(sigma[[q]][k,r]) ) 
          {
                  mu[[q]][k+1,r] <- sign( sum(DATA[,r]) ) * ( abs(sum(DATA[,r]))/N -(Lambda*sigma[[q]][k,r])/N   )
          }
          else
          {
            mu[[q]][k+1,r] <- 0
          }
        }
      }
    }
    
 
    
    if(method =="grouped")
    {
      K <- length(mu)
      N <- dim(DATA)[1]
      p <- dim(DATA)[2]
   
      M_M <- p/km # dimensions for each group
      
      if(M_M%%1==0) 
      {
        for(q in 1:K)
        {  
          #browser()
          for(r in 1:M_M)
          {
            
            index <- (km*(r-1)+1):(km*r) 
            # Check if the MU for this group has already been forced to be 0, if so, keep them 0
            if(all(mu[[q]][k,index]== rep(0, km)))
            {
              mu[[q]][k+1,index] <- 0
            }
            else
            {  
             
              X <- DATA[,index]
              sign_value <- (1- (Lambda*sqrt(km))/sqrt(sum( ( apply( (as.matrix(X)%*%as.matrix(solve(diag(sigma[[q]][k,index])))) , 2, sum ) )^2 )) )
              
              if( sign_value  > 0 )
              {
                coef <- (Lambda*sqrt(km)) / (N*(sqrt(sum((mu[[q]][k,index])^2)))) # sigma depends on mu
                tilda_mu <- apply(X, 2, sum)/N
                mu[[q]][k+1,index] <- solve( diag(km) + coef*diag(sigma[[q]][k,index]) )%*%( tilda_mu )
              }
              else
              {
                mu[[q]][k+1,index] <- 0
              }
            }
          }
        }
        
      }
      
      else
      {
        print("Length of Group is not integer, change the number of groups!")
        break
      }

    }
    
    if(method == "variable")
    {
    
      K <- length(mu)
      N <- dim(DATA)[1]
      p <- dim(DATA)[2]
      
      
      for(r in 1:p)
      {  
        max_compare <- rep(NA, K)
        for(q in 1:K)
        {
          
         mu[[q]][k+1,r] <-  sum(DATA[,r])/N 
         max_compare[q] <- mu[[q]][k+1,r]
        }
        
        #browser()
        # Check the cluster with the largest absolute mean component
        max_index <- which(abs(max_compare) == max(abs(max_compare)))
        
        if(length(max_index)!= 1 )
        {
          max_index <- max_index[1]
        }
        
        # penalize the largest the mean component
        if(  Lambda < abs(sum(DATA[,r]))/(sigma[[max_index]][k,r]) )
        {
          mu[[max_index]][k+1,r] <- sign( sum(DATA[,r]) ) * ( abs(sum(DATA[,r]))/N -(Lambda*sigma[[max_index]][k,r])/N   )
        }
        else
        {
          for(v in 1:K)
          {
            mu[[v]][k+1,r] <- 0
          }
          
          
        }
        
      }
      
      
    }
    
  

  return(list(MU = mu))
}



MU_cal_noEM_adaptive<- function(MU_tilda,DATA, DELTA, mu, sigma, Lambda, k, method = "individual", Adaptive = FALSE, gamma = 1, km = km)
{
  
  if(method == "individual_equal")
  {
    K <- length(mu)
    N <- dim(DATA)[1]
    p <- dim(DATA)[2]
    
    
    
    for(q in 1:K)
    {  
      for(r in 1:p)
      {
        
        
        
        if(abs(MU_tilda[[q]][r])==0){weight_para<-Inf}else{
          
          weight_para <- 1/ abs(MU_tilda[[q]][r])^gamma
        }
        
        
        # Deal with the infinity weight
        if(weight_para == Inf)
        {
          weight_para <- 1e238
        } 
        
        
        if( Lambda*weight_para < abs(sum(DATA[,r]))/(sigma[[q]][k,r]) )
        {
        
          mu[[q]][k+1,r] <- sign( sum(DATA[,r]) ) * ( abs(sum(DATA[,r]))/N -(Lambda*weight_para*sigma[[q]][k,r])/N   )
          
          }
        else
        {
          mu[[q]][k+1,r] <- 0
        }
      }
    }
  }
  
  
  
  if(method =="grouped")
  {
    K <- length(mu)
    N <- dim(DATA)[1]
    p <- dim(DATA)[2]
    
 
    M_M <- p/km # dimensions for each group
    
    if(M_M%%1==0) 
    {
      for(q in 1:K)
      {  
     
        for(r in 1:M_M)
        {
          
 
          index <- (km*(r-1)+1):(km*r) 
 
          if(all(mu[[q]][k,index]== rep(0, km)))
          {
            mu[[q]][k+1,index] <- 0
          }
          else
          {  
           
            X <- DATA[,index]
            
            weight_para <- 1/ (abs( base::norm(MU_tilda[[q]][index] ,type="2") ))^gamma

            if(weight_para == Inf)
            {
              weight_para <- 1e238
            }
            
           
            sign_value <- (1- (Lambda*weight_para*sqrt(km))/sqrt(sum( ( apply( (as.matrix(X)%*%as.matrix(solve(diag(sigma[[q]][k,index])))) , 2, sum ) )^2 )) )
            
            if( sign_value  > 0 )
            {
              coef <- (Lambda*weight_para*sqrt(km)) / (N*(sqrt(sum((mu[[q]][k,index])^2)))) # sigma depends on mu
              tilda_mu <- apply(X, 2, sum)/N
              mu[[q]][k+1,index] <- solve( diag(km) + coef*diag(sigma[[q]][k,index]) )%*%( tilda_mu )
            }
            else
            {
              mu[[q]][k+1,index] <- 0
            }
          }
        }
      }
      
    }
    
    else
    {
      print("Length of Group is not integer, change the number of groups!")
      break
    }

  }
  
  if(method == "variable")
  {
  
    
    K <- length(mu)
    N <- dim(DATA)[1]
    p <- dim(DATA)[2]
    
    
    for(r in 1:p)
    {  
      
    
      
      weight_para <- rep(NA, K)
      max_compare <- rep(NA, K)
      for(q in 1:K)
      {
        
          weight_para[q] <- 1/ (abs( MU_tilda[[q]][r]))^gamma 
    
        
        # Deal with the infinity weight
        if(weight_para[q] == Inf)
        {
          weight_para[q] <- 1e238
        }
          
        mu[[q]][k+1,r] <-  sum(DATA[,r])/N 
        max_compare[q] <- mu[[q]][k+1,r]
      }
      
      
      
      
      max_index <- which(abs(max_compare) == max(abs(max_compare)))
      
      if(length(max_index)!= 1 )
      {
        max_index <- max_index[1]
      }
      
   
      if( Lambda*weight_para [max_index]< abs(sum(DATA[,r]))/(sigma[[max_index]][k,r])  )
      {
           mu[[max_index]][k+1,r] <- sign( sum(DATA[,r]) ) * ( abs(sum(DATA[,r]))/N -(Lambda*weight_para[max_index]*sigma[[max_index]][k,r])/N   )
        
        }
      else
      {

        for(v in 1:K)
        {
          mu[[v]][k+1,r] <- 0
        }
        
        
      }
      
    }
    
    
  }
  
  
  
  return(list(MU = mu))
}

# Calculate the mean vector
# DATA: dataset, each row is an observation
# sigma: sigma for each dimension, This sigma stores results from all iterations
# mu: mu for each dimension, This sigma stores results from all iterations
# Lambda: tuning parameter
# k: number of iteration
# DELTA: probality of being in one of the cluster for each observation


MU_cal <- function(DATA, DELTA, mu, sigma, Lambda, k, method = "individual", Adaptive = FALSE, gamma = 1, km = km)
{
 

    if(method == "individual_equal")
    {
      K <- length(mu)
      N <- dim(DATA)[1]
      p <- dim(DATA)[2]
      
      
      
      for(q in 1:K)
      {  
        for(r in 1:p)
        {
          if( Lambda < abs(sum(DELTA[,q]*DATA[,r]))/(sigma[[q]][k,r]) ) 
          {
            mu[[q]][k+1,r] <- ( sum(DELTA[,q]*DATA[,r])/sum(DELTA[,q]) ) * ( (1-(Lambda*sigma[[q]][k,r])/(abs(sum(DELTA[,q]*DATA[,r])) ) )  )
            if( is.na( mu[[q]][k+1,r]) == TRUE){
              print("sum DELTA*Data == 0")
              print("Lambda")
              print(Lambda)
              print("K")
              print(q)
              mu[[q]][k+1,r] <- 1e238
              
            }
            }
          else
          {
            mu[[q]][k+1,r] <- 0
          }
        }
      }
    }
    
    if(method == "individual")
    {
      K <- length(mu)
      N <- dim(DATA)[1]
      p <- dim(DATA)[2]
      
      
      
      for(q in 1:K)
      {  
        for(r in 1:p)
        {
          if( (1-(Lambda*sigma[[q]][k,r])/abs(sum(DELTA[,q]*DATA[,r])) ) > 0 )
          {
            mu[[q]][k+1,r] <- ( sum(DELTA[,q]*DATA[,r])/sum(DELTA[,q]) ) * ( (1-(Lambda*sigma[[q]][k,r])/(abs(sum(DELTA[,q]*DATA[,r])) ) )  )
          }
          else
          {
            mu[[q]][k+1,r] <- 0
          }
        }
      }
    }
    
    if(method =="grouped")
    {
      K <- length(mu)
      N <- dim(DATA)[1]
      p <- dim(DATA)[2]
      
      #browser()
      M_M <- p/km # dimensions for each group
      
      if(M_M%%1==0) 
      {
        for(q in 1:K)
        {  
          #browser()
          for(r in 1:M_M)
          {
            
            #browser()
            index <- (km*(r-1)+1):(km*r) 
            # Check if the MU for this group has already been forced to be 0, if so, keep them 0
            if(all(mu[[q]][k,index]== rep(0, km)))
            {
              mu[[q]][k+1,index] <- 0
            }
            else
            {  
              #browser()
              X <- DATA[,index]
             sign_value <- (1- (Lambda*sqrt(km))/sqrt(sum( ( apply( DELTA[,q]*(as.matrix(X)%*%as.matrix(solve(diag(sigma[[q]][k,index])))) , 2, sum ) )^2 )) )
              
              if( sign_value  > 0 )
              {
                coef <- (Lambda*sqrt(km)) / sum(DELTA[,q]*(sqrt(sum((mu[[q]][k,index])^2)))) # sigma depends on mu
                tilda_mu <- apply(DELTA[,q]*X, 2, sum)/sum(DELTA[,q])
                mu[[q]][k+1,index] <- solve( diag(km) + coef*diag(sigma[[q]][k,index]) )%*%( tilda_mu )
              }
              else
              {
                mu[[q]][k+1,index] <- 0
              }
            }
          }
        }
        
      }
      
      else
      {
        print("Length of Group is not integer, change the number of groups!")
        break
      }
 
    }
    
    if(method == "variable")
    {
     
      
      K <- length(mu)
      N <- dim(DATA)[1]
      p <- dim(DATA)[2]
      
      
      for(r in 1:p)
      {  
        max_compare <- rep(NA, K)
        for(q in 1:K)
        {
          if(sum(DELTA[,q])==0){ mu[[q]][k+1,r] <- 1e238 # should be inf?
          }else{
          mu[[q]][k+1,r] <- ( sum(DELTA[,q]*DATA[,r])/sum(DELTA[,q]) )}
          max_compare[q] <- mu[[q]][k+1,r]
        }
        
       
        max_index <- which(abs(max_compare) == max(abs(max_compare)))
        
        if(length(max_index)!= 1 )
        {
          max_index <- max_index[1]
        }
        
        # penalize the largest the mean component
        if( (1-(Lambda*sigma[[max_index]][k,r])/abs(sum(DELTA[,max_index]*DATA[,r])) ) > 0 )
        {
          mu[[max_index]][k+1,r] <- ( sum(DELTA[,max_index]*DATA[,r])/sum(DELTA[,max_index]) ) * ( (1-(Lambda*sigma[[max_index]][k,r])/(abs(sum(DELTA[,max_index]*DATA[,r])) ) )  )
        }
        else
        {
         
          for(v in 1:K)
          {
            mu[[v]][k+1,r] <- 0
          }
          
        }
        
      }
      
    }
    
  


    

  return(list(MU = mu))
}

MU_cal_adaptive<- function(MU_tilda,DATA, DELTA, mu, sigma, Lambda, k, method = "individual", gamma = 1, km = km)
{
  

    if(method == "individual_equal")
    {
      K <- length(MU_tilda)
      N <- dim(DATA)[1]
      p <- dim(DATA)[2]
      
      for(q in 1:K)
      {  
        for(r in 1:p)
        {
          
          if(abs(MU_tilda[[q]][r])==0){weight_para<-Inf}else{
          
            weight_para <- 1/ abs(MU_tilda[[q]][r])^gamma
          }
          
        
      
          if(weight_para == Inf)
          {
            weight_para <- 1e238
    
          } 
          

          
         if(  Lambda < abs(sum(DELTA[,q]*DATA[,r]))/(sigma[[q]][k,r]* weight_para )  )
          {
            mu[[q]][k+1,r] <- ( sum(DELTA[,q]*DATA[,r])/sum(DELTA[,q]) ) * ( (1-(Lambda* weight_para *sigma[[q]][k,r])/(abs(sum(DELTA[,q]*DATA[,r])) ) )  )
            if(is.na(mu[[q]][k+1,r])==TRUE){
              
              print("In adaptive step, sum Delat*data==0,individual inf mu")
              print("Lambda")
              print(Lambda)
              print("Gamma")
              print(gamma)
              print("K")
              print(q)
              
              mu[[q]][k+1,r] <- 1e238 #Inf
              
              
            }
            
          }
          else
          {
            mu[[q]][k+1,r] <- 0
          }
        }
      }
    }
    
    if(method == "individual")
    {
      K <- length(mu)
      N <- dim(DATA)[1]
      p <- dim(DATA)[2]
      
      
      
      for(q in 1:K)
      {  
        for(r in 1:p)
        {
          weight_para <- 1/ (abs(MU_tilda[[q]][r]))^gamma
          
          # Deal with the infinity weight
          if(weight_para == Inf)
          {
            weight_para <- 1e238
          }
          
          if( (1-(Lambda* weight_para* sigma[[q]][k,r])/abs(sum(DELTA[,q]*DATA[,r])) ) > 0 )
          {
            mu[[q]][k+1,r] <- ( sum(DELTA[,q]*DATA[,r])/sum(DELTA[,q]) ) * ( (1-(Lambda* weight_para* sigma[[q]][k,r])/(abs(sum(DELTA[,q]*DATA[,r])) ) )  )
          }
          else
          {
            mu[[q]][k+1,r] <- 0
          }
        }
      }
    }
    
    if(method =="grouped")
    {
      K <- length(MU_tilda)
      N <- dim(DATA)[1]
      p <- dim(DATA)[2]

      M_M <- p/km # dimensions for each group
      
      if(M_M%%1==0) 
      {
        for(q in 1:K)
        {  

          for(r in 1:M_M)
          {
            
    
            index <- (km*(r-1)+1):(km*r) 
    
            if(all(mu[[q]][k,index]== rep(0, km)))
            {
              mu[[q]][k+1,index] <- 0
            }
            else
            {  
              
              X <- DATA[,index]
              
              weight_para <- 1/ (abs( base::norm(MU_tilda[[q]][index] ,type="2") ))^gamma
      
              if(weight_para == Inf)
              {
                weight_para <- 1e238
              }
              
           
              sign_value <- (1- (Lambda* weight_para* sqrt(km))/sqrt(sum( ( apply( DELTA[,q]*(as.matrix(X)%*%as.matrix(solve(diag(sigma[[q]][k,index])))) , 2, sum ) )^2 )) )
              
              if( sign_value  > 0 )
              {
                coef <- (Lambda* weight_para* sqrt(km)) / sum(DELTA[,q]*(sqrt(sum((mu[[q]][k,index])^2)))) # sigma depends on mu
                tilda_mu <- apply(DELTA[,q]*X, 2, sum)/sum(DELTA[,q])
                mu[[q]][k+1,index] <- solve( diag(km) + coef*diag(sigma[[q]][k,index]) )%*%( tilda_mu )
              }
              else
              {
                mu[[q]][k+1,index] <- 0
              }
            }
          }
        }
        
      }
      
      else
      {
        print("Length of Group is not integer, change the number of groups!")
        break
      }
     
    }
    
    if(method == "variable")
    {
      
      K <- length(mu)
      N <- dim(DATA)[1]
      p <- dim(DATA)[2]
      
      
      for(r in 1:p)
      {  
      
        weight_para <- rep(NA, K)
        max_compare <- rep(NA, K)
        for(q in 1:K)
        {

            weight_para[q] <- 1/ (abs( MU_tilda[[q]][r]))^gamma 
          
          # Deal with the infinity weight
          if(weight_para[q] == Inf)
          {
            weight_para[q] <- 1e238
          }
            
            
         mu[[q]][k+1,r] <- ( sum(DELTA[,q]*DATA[,r])/sum(DELTA[,q]) )
         max_compare[q] <- mu[[q]][k+1,r]
        }
        
        max_index <- which.max(abs(max_compare))
        
        
        if(length(max_index)!= 1 )
        {
          max_index <- max_index[1]
        }
        
        # penalize the largest the mean component
        if( (1-(Lambda*weight_para[max_index] * sigma[[max_index]][k,r])/abs(sum(DELTA[,max_index]*DATA[,r])) ) > 0 )
        {
          mu[[max_index]][k+1,r] <- ( sum(DELTA[,max_index]*DATA[,r])/sum(DELTA[,max_index]) ) * ( (1-(Lambda*weight_para[max_index] *sigma[[max_index]][k,r])/(abs(sum(DELTA[,max_index]*DATA[,r])) ) )  )
        }else{

          for(v in 1:K)
          {
            mu[[v]][k+1,r] <- 0
          }
          
          
        }
        
      
      
      

      
      
    }
    
  }
  

  return(list(MU = mu))
}



# Calculate the sigma vector
# DATA: dataset, each row is an observation
# DELTA: probality of being in one of the cluster for each observation
# sigma: sigma for each dimension, This sigma stores results from all iterations
# mu: for each dimension, This sigma stores results from all iterations
# k: number of iteration
# Lambda2: needed for the penalized variance


SIGMA_cal_noEM <- function(DATA, DELTA, mu, sigma, Lambda2 = 5, k, sigma_method = "unequal")
{

  K <- length(mu)
  N <- dim(DATA)[1]
  p <- dim(DATA)[2]
  
  

    for(l in 1:p)
    {
      # calculated the MLE for SIGMA for the first cluster
      
      # Initiate the simga
      sigma_accumulate <- 0
      
      for(r in 1:K)
      {
        
        sigma_accumulate <- sigma_accumulate + sum((DATA[,l] - mu[[r]][k,l])^2 )/N  # need more look
      }
      for(q in 1:K)
      {
        sigma[[q]][k+1,l] <- sigma_accumulate
      }
      
    }
  
  
  return(list(SIGMA = sigma))
  
}


SIGMA_cal <- function(DATA, DELTA, mu, sigma, Lambda2 = 5, k, sigma_method = "unequal")
{
  
  K <- length(mu)
  N <- dim(DATA)[1]
  p <- dim(DATA)[2]
  
  if(sigma_method == "unequal")
  {
    # For unequal variance, we define b and c vectors in Pan and Shen's paper
    B <- rep(NA, K)
    C <- matrix(rep(NA, K*p), ncol = p )
    
    for(i in 1:K)
    {
      B[i] <- sum(DELTA[,i]/2)
      
      for(j in 1:p)
      {
        C[i,j] <-  sum(DELTA[,i]*(DATA[,j] - mu[[i]][k,j])^2 )/2
      }
    }
    
    for(i in 1:K)
    {
      for(j in 1:p)
      {
        if( sign(abs(B[i]-C[i,j])-Lambda2 ) > 0   )
        {
          sigma[[i]][k+1,j] <-  ( ( C[i,j]/B[i] )/(1+ Lambda2*sign(C[i,j]-B[i])/B[i]) - 1 ) +1
        }
        else 
        {
          sigma[[i]][k+1,j] <- 1
        }
      }
    }
    
  }
  

  
  if(sigma_method == "equal")
  {
    for(l in 1:p)
    {
      # calculated the MLE for SIGMA for the first cluster
      
      # Initiate the simga
      sigma_accumulate <- 0
      
      for(r in 1:K)
      {
       
        sigma_accumulate <- sigma_accumulate + sum(DELTA[,r]*(DATA[,l] - mu[[r]][k,l])^2 )/N  # need more look
      }
      for(q in 1:K)
      {
        sigma[[q]][k+1,l] <- sigma_accumulate
      }
      
    }
  }
  
  return(list(SIGMA = sigma))
  
}





########################## E-M algorithm####################################
# 
# DATA: dataset for E-M Gaussian Mixture model
# K: Number of clusters 
# M: Maximum number of iterations(for high dimensional, could be time consumed)
# method: "individual", "grouped" for different penalty
# km: dimension for each group if using "grouped" method, in our case, it should equal to the number of component 
#     we choose in the MFPCA



# method is one of the 'individual', 'variable' and 'grouped', representing three penalty terms respectively
# sigma_method is one of the 'equal' and 'unequal', standing for the variances across each cluster.
# Lambda2 is only used in sigma_method == 'unequal' case, to penalize variance towards 1. 
# M is the maximum iteration
# km is only used in method == 'grouped' case, representing the length of each mean component group. 

EM_simple_first_step_zhongnan <- function(DATA, K, Lambda, Lambda2, sigma_method = "unequal", method = "individual", Adaptive = FALSE, gamma = 1, M, 
                                 km = 3, min_percentage = 0.15,PCA_exp=3,
                                 loglikn = FALSE)
{
  p <- dim(DATA)[2]
  N <- dim(DATA)[1]
  # Set initial values for Parameters
  PI <- matrix(rep(NA, K*M), ncol=K) # probability of observations from different group from each E-M process.
  MU <- rep(list(matrix(rep(NA, p*M),ncol=p )),K) # mean vector for each group from each E-M process. 
  SIGMA <- rep(list(matrix(rep(0, p*M),ncol=p )),K)  # sigma vector for each group. (has to be 0, since we use it as an initial to sum)
  IND_matrix <- matrix(rep(FALSE, p*K), ncol=p)
  IND <- rep(FALSE, p) # vector of length p, if FALSE, means the variabel can not be eliminated
  
  
  for( k in 1:(M-1))
  {
    
    if(k==1)
    {
  
      
      repeat{
        PI[k, ] <- Initial_Values_old(DATA, K=K)$Alpha
        mu_initial <- Initial_Values_old(DATA, K=K)$Mu # generate initial value containing K rows 
        if(min(PI[k, ]) > min_percentage)
        {break}
      }
      

      for(i in 1:K)
      {
        MU[[i]][k,] <- as.numeric(mu_initial[i,])
        SIGMA[[i]][k,] <- as.numeric(Initial_Values_old(DATA, K=K)$Sigma)  # assume same sigma for all cluster
      }
      

      DELTA <- PI_cal_old(DATA, PI= PI,mu=MU,sigma = SIGMA, k)$DELTA
      PI[k+1,] <- PI_cal_old(DATA, PI,SIGMA, MU, k)$Pi
      
      
      MU <- MU_cal(DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda = Lambda, k = k, method = method, Adaptive = Adaptive, gamma = gamma, km = km)$MU
     SIGMA <- SIGMA_cal(DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda2 = Lambda2, k = k, sigma_method = sigma_method)$SIGMA
      
      
    }
    
    else
    {
      
      DELTA <- PI_cal_old(DATA, PI, SIGMA, MU, k)$DELTA
      PI[k+1,] <- PI_cal_old(DATA, PI, SIGMA, MU, k)$Pi
      
      MU <- MU_cal(DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda = Lambda, k = k, method = method, Adaptive = Adaptive, gamma = gamma, km = km)$MU
      SIGMA <- SIGMA_cal(DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda2 = Lambda2, k = k, sigma_method = sigma_method)$SIGMA
      
      
    }
    
    if(max(abs(PI[k+1, ] - PI[k,])) < 5e-5 )
    {
      criteria <- rep(FALSE, K)
      for(i in 1:K)
      {
        if(max(abs(MU[[i]][k+1, ] - MU[[i]][k, ] )) < 5e-5){criteria[i] = TRUE}
      }  
      
      if(all(criteria))
      {
        break
      }
    }else{}
    
  }
  
  #################################
  # Calculate the cluster result
  cluster_result <- rep(NA, dim(DELTA)[1])
  for(i in 1:dim(DELTA)[1])
  {
    cluster_result[i] <- which(DELTA[ i,]==max(DELTA[i,]) ) # For each oberservation 
  }
  
 
  #################################
  # Calcualte the likelihood 

  Likelihood <- 0
  logLik <- 0
  
  
  zero_prob_ind <- matrix(rep(NA, N*K) , ncol=K)
  
  for(i in 1:N)
  {
    
    for(j in 1:K)
    {
      prob <- exp( mpfr(dmvnorm(x=DATA[i,], mean=MU[[j]][k,], sigma = diag(sqrt(SIGMA[[j]][k,])), log = T) , 300) )
      lik <- mpfr(PI[k,j],300)*prob
      Likelihood <- Likelihood + lik
      if(lik==0){zero_prob_ind[i,j] = 1}else{zero_prob_ind[i,j]=0}
    }
    

    logLik <- logLik + log(Likelihood)
  }
 
  
  if(loglikn==TRUE){
    logLikBIC <- logLik/sqrt(N) 
  }else if(loglikn==FALSE){
    logLikBIC <- logLik
  }
  

  if(method == "individual_equal")
  {
    # Initial penalty term
    penalty <- 0
    for(i in 1:p)
    {
      for(j in 1:K)
      {
        penalty <- penalty + Lambda*abs(MU[[j]][k,i])
      }
    }
  }

  if(method == "individual")
  {
    penalty <- 0
    for(i in 1:p)
    {
      for(j in 1:K)
      {
        penalty <- penalty + Lambda*abs(MU[[j]][k,i]) + Lambda2*abs(log(SIGMA[[j]][k,i]))
      }
    }
  }

  if(method == "grouped")
  {

    penalty <- 0
    M_M <- p/km 
    for(q in 1:K)
    {
 
      for(r in 1:M_M)
      {

        index <- (km*(r-1)+1):(km*r)
        penalty <- penalty + Lambda*sqrt(km)*sqrt( sum( (MU[[q]][k,index])^2 ) )
      }
    }

  }

  if(method == "variable")
  {
    # Initial penalty term
    penalty <- 0
    for(i in 1:p)
    {
      mu_vec_compare <- rep(NA, K)
      for(j in 1:K)
      {
        mu_vec_compare[j] <- abs(MU[[j]][k,i])
      }

      penalty <- penalty + Lambda*max(abs(mu_vec_compare))
    }

  }
 
  
  # Calculate the BIC
  Num_zero_mean <- 0
  for(m in 1:K){Num_zero_mean <- Num_zero_mean+ length( which(MU[[m]][k,]==0)) }
  BIC <- -2*logLikBIC + log(N)*( p+K+p*K-1- Num_zero_mean)
  

  ALL_ZERO <- rep(FALSE, K)
  for(q in 1:K)
  {
    if(all(MU[[q]][k,] == rep(0,p)))
    {
      ALL_ZERO[q] <- TRUE
    }
  }
  
  if(all(ALL_ZERO))
  {
    print("All mean elements are 0!")
    all_zero = TRUE
  }
  else
  {
    all_zero = FALSE
  }
  

  for(i in 1:K)
  {
    IND_matrix[i,which(MU[[i]][k+1,]==0) ] <- TRUE
  }
  
  if(method == "variable" || method =="individual"|| method =="individual_equal"){
  for(j in 1:p)
  {
    if(all( IND_matrix[,j] ) )
    {
      IND[j] <- TRUE
    }
  }
  }else if( method == "grouped"){
    
    ind<-seq(1,p,by=PCA_exp)
    
    for(j in ind){
      if(all( IND_matrix[,j:(j+(PCA_exp-1))] ) )
      {
        IND[j:(j+(PCA_exp-1))] <- TRUE
      }
    }
    
      }
    
    
  
  
  
  return(list(PI = PI, MU = MU, SIGMA = SIGMA, DELTA = DELTA, Stop_id = k, Cluster = cluster_result, 
              logLik = logLik, penalty = penalty, BIC =BIC,  
              all_zero = all_zero, zero_prob_ind = zero_prob_ind, VAR_removed = which(IND) ))
}


EM_simple_first_step <- function(DATA, K, Lambda, Lambda2, sigma_method = "unequal", method = "individual", M, 
                                          km = 3, min_percentage = 0.15,PCA_exp=3)
{
  p <- dim(DATA)[2]
  N <- dim(DATA)[1]
  # Set initial values for Parameters
  PI <- matrix(rep(NA, K*M), ncol=K) # probability of observations from different group from each E-M process.
  MU <- rep(list(matrix(rep(NA, p*M),ncol=p )),K) # mean vector for each group from each E-M process. 
  SIGMA <- rep(list(matrix(rep(0, p*M),ncol=p )),K)  # sigma vector for each group. (has to be 0, since we use it as an initial to sum)
  IND_matrix <- matrix(rep(FALSE, p*K), ncol=p)
  IND <- rep(FALSE, p) # vector of length p, if FALSE, means the variabel can not be eliminated
  
  
  for( k in 1:(M-1))
  {
 
    if(k==1)
    {
      
      count_ini<-0
      repeat{
      
        PI[k, ] <- Initial_Values_old(DATA, K=K)$Alpha
        mu_initial <- Initial_Values_old(DATA, K=K)$Mu # generate initial value containing K rows 
        if(min(PI[k, ]) > min_percentage){break}
        
        count_ini<-count_ini+1
        


        if(count_ini>3){
          print("initial repeated >3.,min_percentage halfed")
          min_percentage=min_percentage/2}
        
        if(count_ini>5){
          print("initial repeated >5.")
          print(min_percentage)
          print(K)
          print(Lambda)
        
          break
          
        }
      
        
        
          }
             

      for(i in 1:K)
      {
        MU[[i]][k,] <- as.numeric(mu_initial[i,])
        SIGMA[[i]][k,] <- as.numeric(Initial_Values_old(DATA, K=K)$Sigma)  # assume same sigma for all cluster
      }
      
      
      DELTA <- PI_cal_old(DATA, PI= PI,mu=MU,sigma = SIGMA, k)$DELTA
      PI[k+1,] <- PI_cal_old(DATA, PI,SIGMA, MU, k)$Pi
      
      if(K==1){
      MU <- MU_cal_noEM(DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda = Lambda, k = k, method = method, Adaptive = FALSE, gamma = gamma, km = km)$MU
      SIGMA <- SIGMA_cal_noEM(DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda2 = Lambda2, k = k, sigma_method = sigma_method)$SIGMA
      }else{
        
        MU <- MU_cal(DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda = Lambda, k = k, method = method, Adaptive = FALSE, gamma = gamma, km = km)$MU
        SIGMA <- SIGMA_cal(DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda2 = Lambda2, k = k, sigma_method = sigma_method)$SIGMA
        
      }
    }
    
    else
    {
      
      
      
      if(K==1){
        MU <- MU_cal_noEM(DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda = Lambda, k = k, method = method, Adaptive = FALSE, gamma = gamma, km = km)$MU
        SIGMA <- SIGMA_cal_noEM(DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda2 = Lambda2, k = k, sigma_method = sigma_method)$SIGMA
        
        if(max(abs(MU[[K]][k+1, ] - MU[[K]][k, ] )) < 5e-5){ break }

          

        
        
        
      }else if(K>1){
        
        
      DELTA <- PI_cal_old(DATA, PI, SIGMA, MU, k)$DELTA
      PI[k+1,] <- PI_cal_old(DATA, PI, SIGMA, MU, k)$Pi
      MU <- MU_cal(DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda = Lambda, k = k, method = method, Adaptive = FALSE, gamma = gamma, km = km)$MU
      SIGMA <- SIGMA_cal(DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda2 = Lambda2, k = k, sigma_method = sigma_method)$SIGMA
      
      
      if(max(abs(PI[k+1, ] - PI[k,])) < 5e-5 )
      {
        criteria <- rep(FALSE, K)
        for(i in 1:K)
        {
          if(max(abs(MU[[i]][k+1, ] - MU[[i]][k, ] )) < 5e-5){criteria[i] = TRUE}
        }  
        
        if(all(criteria))
        {
 
          break
        }
      }else{}
      
      }
    }
    
   
    
  }
  
  #################################
  # Calculate the cluster result
  
  if(K==1){cluster_result<-rep(1,N)}else{
  cluster_result <- rep(NA, N)
  for(i in 1:N)
  {
    cluster_result[i] <- which(DELTA[ i,]==max(DELTA[i,]) ) # For each oberservation 
  }
  
  }
  #################################
  # Calcualte the likelihood 
  

  Likelihood <- 0
  logLik <- 0
  
  
  zero_prob_ind <- matrix(rep(NA, N*K) , ncol=K)
  
  for(i in 1:N)
  {
    
    for(j in 1:K)
    {
      prob <- exp( mpfr(dmvnorm(x=DATA[i,], mean=MU[[j]][k,], sigma = diag(sqrt(SIGMA[[j]][k,])), log = T) , 300) )
       lik <- mpfr(PI[k,j],300)*prob
       Likelihood <- Likelihood + lik
      if(lik==0){zero_prob_ind[i,j] = 1}else{zero_prob_ind[i,j]=0}
    }
    
    
    logLik <- logLik + log(Likelihood)
  }
  
  
  
  # Calculate the BIC
  Num_zero_mean <- 0
  for(m in 1:K){Num_zero_mean <- Num_zero_mean+ length( which(MU[[m]][k,]==0)) }
  BIC <- -2*logLik + log(N*p)*( p+K+p*K-1- Num_zero_mean)

  ALL_ZERO <- rep(FALSE, K)
  for(q in 1:K)
  {
    if(all(MU[[q]][k,] == rep(0,p)))
    {
      ALL_ZERO[q] <- TRUE
    }
  }
  
  if(all(ALL_ZERO))
  {
    print("All mean elements are 0!")
    all_zero = TRUE
  }
  else
  {
    all_zero = FALSE
  }
  

  for(i in 1:K)
  {
    IND_matrix[i,which(MU[[i]][k+1,]==0) ] <- TRUE
  }
  
  if(method == "variable" || method =="individual"|| method =="individual_equal"){
    for(j in 1:p)
    {
      if(all( IND_matrix[,j] ) )
      {
        IND[j] <- TRUE
      }
    }
  }else if( method == "grouped"){
    
    ind<-seq(1,p,by=PCA_exp)
    
    for(j in ind){
      if(all( IND_matrix[,j:(j+(PCA_exp-1))] ) )
      {
        IND[j:(j+(PCA_exp-1))] <- TRUE
      }
    }
    
  }
  

  
  
  return(list(PI = PI, MU = MU, SIGMA = SIGMA,  Stop_id = k,logLik = logLik,  BIC =BIC ))
}


EM_adaptive<-function(MU_tilda,DATA,K,Lambda,Lambda2,sigma_method="unequal",method,gamma, M, km=3,min_percentage,PCA_exp=3){
  
  
 
    p <- dim(DATA)[2]
    N <- dim(DATA)[1]
    # Set initial values for Parameters
    PI <- matrix(rep(NA, K*M), ncol=K) # probability of observations from different group from each E-M process.
    MU <- rep(list(matrix(rep(NA, p*M),ncol=p )),K) # mean vector for each group from each E-M process. 
    SIGMA <- rep(list(matrix(rep(0, p*M),ncol=p )),K)  # sigma vector for each group. (has to be 0, since we use it as an initial to sum)
    IND_matrix <- matrix(rep(FALSE, p*K), ncol=p)
    IND <- rep(FALSE, p) # vector of length p, if FALSE, means the variabel can not be eliminated
    
    
    for( k in 1:(M-1))
    {
      
      if(k==1)
      {

        count_ini <-0
        repeat{
          PI[k, ] <- Initial_Values_old(DATA, K=K)$Alpha
          mu_initial <- Initial_Values_old(DATA, K=K)$Mu # generate initial value containing K rows 
          if(min(PI[k, ]) > min_percentage)
          {break}
        

        count_ini<-count_ini+1

        if(count_ini>3){
          print("initial repeated >3.,min_percentage halfed")
          min_percentage=min_percentage/2}
        
        if(count_ini>5){
          print("initial repeated >5.")
          print(min_percentage)
          print(Lambda)
          print(gamma)
          break
          
        }
        }
        

        for(i in 1:K)
        {
          MU[[i]][k,] <- as.numeric(mu_initial[i,])
          SIGMA[[i]][k,] <- as.numeric(Initial_Values_old(DATA, K=K)$Sigma)  # assume same sigma for all cluster
        }
        
        
        DELTA <- PI_cal_old(DATA, PI= PI,mu=MU,sigma = SIGMA, k)$DELTA
        PI[k+1,] <- PI_cal_old(DATA, PI,SIGMA, MU, k)$Pi
        
        if(K==1){
        MU <- MU_cal_noEM_adaptive(MU_tilda, DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda = Lambda, k = k, method = method, gamma = gamma, km = km)$MU
        SIGMA <- SIGMA_cal(DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda2 = Lambda2, k = k, sigma_method = sigma_method)$SIGMA
        }else{
          
          MU <- MU_cal_adaptive(MU_tilda, DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda = Lambda, k = k, method = method, gamma = gamma, km = km)$MU
          SIGMA <- SIGMA_cal(DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda2 = Lambda2, k = k, sigma_method = sigma_method)$SIGMA
          
        }
        
      }
      
      else
      {  
        if(k==(M-1)){"EM adaptive reaches max iteration number"}
        
        if(K==1){
          MU <- MU_cal_noEM_adaptive(MU_tilda, DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda = Lambda, k = k, method = method, Adaptive = FALSE, gamma = gamma, km = km)$MU
          SIGMA <- SIGMA_cal_noEM(DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda2 = Lambda2, k = k, sigma_method = sigma_method)$SIGMA
          
          if(max(abs(MU[[K]][k+1, ] - MU[[K]][k, ] )) < 5e-5){ break }
    
        }else if(K>1){
        
        DELTA <- PI_cal_old(DATA, PI, SIGMA, MU, k)$DELTA
        PI[k+1,] <- PI_cal_old(DATA, PI, SIGMA, MU, k)$Pi
        
        MU <- MU_cal_adaptive(MU_tilda,DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda = Lambda, k = k, method = method, gamma = gamma, km = km)$MU
        SIGMA <- SIGMA_cal(DATA = DATA, DELTA = DELTA, mu = MU, sigma = SIGMA, Lambda2 = Lambda2, k = k, sigma_method = sigma_method)$SIGMA
        

      
      if(max(abs(PI[k+1, ] - PI[k,])) < 5e-5 )
      {
        criteria <- rep(FALSE, K)
        for(i in 1:K)
        {
          if(max(abs(MU[[i]][k+1, ] - MU[[i]][k, ] )) < 5e-5){criteria[i] = TRUE}
        }  
        
        if(all(criteria))
        {
  
          break
        }
      }else{}
      
        }
      }
    }
    
    #################################
    # Calculate the cluster result
    if(K==1){cluster_result<-rep(1,N)}else{
    cluster_result <- rep(NA, N)
    for(i in 1:N)
    {
      cluster_result[i] <- which(DELTA[ i,]==max(DELTA[i,]) ) # For each oberservation 
    }
    
    }
    #################################
    # Calcualte the likelihood 
    

    Likelihood <- 0
    logLik <- 0
    
    
    zero_prob_ind <- matrix(rep(NA, N*K) , ncol=K)
    
    for(i in 1:N)
    {
      
      for(j in 1:K)
      {
        prob <- exp( mpfr(dmvnorm(x=DATA[i,], mean=MU[[j]][k,], sigma = diag(sqrt(SIGMA[[j]][k,])), log = T) , 300) )
        lik <- mpfr(PI[k,j],300)*prob
        Likelihood <- Likelihood + lik
        if(lik==0){zero_prob_ind[i,j] = 1}else{zero_prob_ind[i,j]=0}
      }
      
      
      logLik <- logLik + log(Likelihood)
    }
    
    
    if(method == "individual_equal")
    {
      # Initial penalty term
      penalty <- 0
      for(i in 1:p)
      {
        for(j in 1:K)
        {
          penalty <- penalty + Lambda*abs(MU[[j]][k,i])
        }
      }

    }
    
    if(method == "individual")
    {
      # Initial penalty term
      penalty <- 0
      for(i in 1:p)
      {
        for(j in 1:K)
        {
          penalty <- penalty + Lambda*abs(MU[[j]][k,i]) + Lambda2*abs(log(SIGMA[[j]][k,i]))
        }
      }
    }
    
    if(method == "grouped")
    {

      penalty <- 0
      M_M <- p/km  # dimensions for each group
      for(q in 1:K)
      {
        for(r in 1:M_M)
        {

          index <- (km*(r-1)+1):(km*r)
          penalty <- penalty + Lambda*sqrt(km)*sqrt( sum( (MU[[q]][k,index])^2 ) )
        }
      }
    }
    
    if(method == "variable")
    {
      # Initial penalty term
      penalty <- 0
      for(i in 1:p)
      {
        mu_vec_compare <- rep(NA, K)
        for(j in 1:K)
        {
          mu_vec_compare[j] <- abs(MU[[j]][k,i])
        }
        
        penalty <- penalty + Lambda*max(abs(mu_vec_compare))
      }
    }

    
    # Calculate the BIC
    Num_zero_mean <- 0
    for(m in 1:K){Num_zero_mean <- Num_zero_mean+ length( which(MU[[m]][k,]==0)) }
    BIC <- -2*logLik + log(N*p)*( p+K+p*K-1- Num_zero_mean)
    

    ALL_ZERO <- rep(FALSE, K)
    for(q in 1:K)
    {
      if(all(MU[[q]][k,] == rep(0,p)))
      {
        ALL_ZERO[q] <- TRUE
      }
    }
    
    if(all(ALL_ZERO))
    {
      print("All mean elements are 0!")
      all_zero = TRUE
    }
    else
    {
      all_zero = FALSE
    }

    for(i in 1:K)
    {
      IND_matrix[i,which(MU[[i]][k+1,]==0) ] <- TRUE
    }
    
    if(method == "variable" || method =="individual"|| method =="individual_equal"){
      for(j in 1:p)
      {
        if(all( IND_matrix[,j] ) )
        {
          IND[j] <- TRUE
        }
      }
    }else if( method == "grouped"){
      
      ind<-seq(1,p,by=PCA_exp)
      
      for(j in ind){
        if(all( IND_matrix[,j:(j+(PCA_exp-1))] ) )
        {
          IND[j:(j+(PCA_exp-1))] <- TRUE
        }
      }
      
    }
    
    
    
    
    
    return(list(PI = PI, MU = MU, SIGMA = SIGMA, DELTA = DELTA, Stop_id = k, Cluster = cluster_result, 
                logLik = logLik, penalty = penalty, BIC =BIC,  
                all_zero = all_zero, zero_prob_ind = zero_prob_ind, VAR_removed = which(IND) ))
  
  
  
}





############################################################Analysis##################################################################

# Plot the EM result for clustering
#
# EM_results: results obtained from EM process (EM_simple function)
# pi : logical, options for plot the pi
# mu : logical, options for plot the mu
# sigma: logical, options for plot the sigma
# sensor: columns/dimensions of sensors users would like to plot for mu and sigma, 
#         can be a single value or a vector

PLOT_EM <- function(EM_results, pi=TRUE, mu = TRUE, sigma = TRUE, sensor = 1 )
{
  K <- dim(EM_results$PI)[2]
  if(pi == TRUE)
  {
    par(mfrow=c(ceiling(K/2), 2))
    for(i in 1:K)
    {  
      id <- EM_results$Stop_id
      plot(EM_results$PI[1:id,i], ylim = c(0,1), type = "l", 
           main = paste("porportion from cluster", i),
           ylab = paste("Pi", sensor[j] ), xlab = "Number of iterations")
    }  
  }
  if(mu == TRUE)
  {
    for(j in 1:length(sensor))
    {  
      dev.new()
      par(mfrow=c(ceiling(K/2), 2))
      for(i in 1:K)
      {  
        id <- EM_results$Stop_id
        plot(EM_results$MU[[i]][1:id,sensor[j]], ylim = c(min(EM_results$MU[[i]][1:id,sensor[j]]),max(EM_results$MU[[i]][1:id,sensor[j]])),
             type = "l", main = paste("Mean trace from cluster", i),
             ylab = paste("Mu", sensor[j]), xlab = "Number of iterations")
        abline(h = EM_results$MU[[i]][id-1,sensor[j]], col = 2, lty = 2)
      } 
    }
  }
  if(sigma == TRUE)
  {
    id <- EM_results$Stop_id
    for(j in 1:length(sensor))
    {  
      dev.new()
      par(mfrow=c(ceiling(K/2), 2))
      for(i in 1:K)
      {  
        plot(EM_results$SIGMA[[i]][1:id,sensor[j]], ylim = c(min(EM_results$SIGMA[[i]][1:id,sensor[j]]),max(EM_results$SIGMA[[i]][1:id,sensor[j]])),
             type = "l", main = paste("Sigma trace from cluster", i),
             ylab = paste("Sigma", sensor[j] ), xlab = "Number of iterations")
        abline(h = EM_results$SIGMA[[i]][id-1,sensor[j]], col = 2, lty = 2)
      }  
    }
  }
}

# Function: Cluster result in originial function format
#
# ORIG_DATA: Original Dataset with a functional formal(Defined in previous MFPCA function)
# EM_results: EM_result returned by EM_simple
# sensor: Indicate the number of sensor user would like to show cluster result from
# observation: 2 values indicate the the start and end index of the proportion of observations users would like to show
# len: Indicate the length of each functions data, returned from MFPCA functions

PLOT_CLUSTER <- function(ORIG_DATA, EM_results, sensor, observation = c(1,20), len)
{
  
  # Create cluster result vector
  cluster_result <- rep(NA, dim(EM_results$DELTA)[1])
  for(i in 1:dim(EM_results$DELTA)[1])
  {
    cluster_result[i] <- which(EM_results$DELTA[ i,]==max(EM_results$DELTA[i,]) ) # For each oberservation 
  }
  
  VAR <-sensor# which sensor/varible we want to plot
  L <- observation[1] # Starting point
  U <- observation[2] # End point
  
  for(j in 1:length(sensor))
  {
    dev.new()
    minY <- min(ORIG_DATA[,VAR[j]])
    maxY <- max(ORIG_DATA[,VAR[j]])
    
    
    plot(1:len,ORIG_DATA[ORIG_DATA[,67]==1,VAR[j]],xlab="grid",ylab="",type='l', ylim=c(minY, maxY), main=paste("Plot for sensor",VAR[j]))
    for(i in U:L){
      cur<-ORIG_DATA[ORIG_DATA[,67]==i,]
      lines(1:len,as.numeric(coredata(cur[,VAR[j]],type="l")), col = cluster_result[i]+1)
    }
  }
  
}


# Summary of the E-M algorithm
#
# EM_result: object returned by EM_simple function
# 
# Return estimated mean and variance vector

summary.EM <- function(EM_result,...)
{
  stop_id <- EM_result$Stop_id
  # Estimated mean vector
  L <- length(EM_result$MU)
  p <- length(EM_result$MU[[1]][stop_id,])
  N <- dim(EM_result$DELTA)[1]
  K <- dim(EM_result$DELTA)[2]
  
  Estimated_Mean <- matrix(rep(NA, p*K), ncol = p)
  colnames(Estimated_Mean) <- paste( "Mu_",1:p, sep="")
  rownames(Estimated_Mean) <- paste( "Cluster", 1:K, sep = "")
  
  for(i in 1:K)
  {
    Estimated_Mean[i,] <- EM_result$MU[[i]][stop_id, ] 
  }
  
  Estimated_Variance <- matrix(rep(NA, p*K), ncol = p)
  colnames(Estimated_Variance) <- paste( "Variance_",1:p, sep="")
  rownames(Estimated_Variance) <- paste( "Cluster", 1:K, sep = "")
  
  for(i in 1:K)
  {
    Estimated_Variance[i,] <- EM_result$SIGMA[[i]][stop_id, ] 
  }
  
  pi <- EM_result$PI[stop_id,]
  
  
  return(list(Mu = Estimated_Mean, Variance = Estimated_Variance, pi = pi, cluster_result = EM_result$Cluster  ))
  
}



# 
# K: number of true clusters
# p: number of dimensions/sensors
# N: number of observations
# pi: a vector with the same length of K, true proportions 
# mu: list of length K including mean vectors of length p for each list element
# sigma: list of length K including sigma vectors of length p for each list element
DATA_GENERATOR <- function(K, p, N, pi, mu, sigma)
{
  DATA <- matrix(rep(NA, p*N), ncol=p)
  if(K == length(pi))
  {
    INDEX <- sample(1:K, N, replace = T, prob = pi)
  }
  if(K != length(pi))
  {
    print("Length of probability is not equal to the number of clusters")
  }
  

  
  for(i in 1:N)
  {
    k <- INDEX[i]
    DATA[i,] <- rmvnorm(1, mean = mu[[k]], sigma = diag(sigma[[k]]) )
  }

  
  return(list(DATA = DATA, INDEX = INDEX, PI = pi ))
}



# Function to choose the best lambda and K
# SIMU_DATA: Multivariate data used for clustering
# MIN_LAMBDA: minimum tuning parameters users would like to try
# MAX_LAMBDA: maximum tuning parameters users would like to try
# STEP_LAMBDA: step length for tuning parameters change when finding the optimal lambda
# K: vector of Ks users would like to try
# M: maximum runing step in E-M algorithm
# method: Either "individual" or "grouped", same usage as EM_simple
# km: length for each group is method is "grouped"
# plot: Whether users choose to plot the result log-likelihood and BIC for all possible lambda and K combinations
# True_cluster: clustering label(INDEX) returned by DATA_GENERATOR function

Lambda_and_cluster <- function(SIMU_DATA, MIN_LAMBDA = 1, MAX_LAMBDA = 15, STEP_LAMBDA = 1, Lambda2 = 10,
                               K = c(2,3,4), M = 100, method = "individual", km = 2, plot = TRUE,
                               True_cluster = NA )
{
  K <- K
  N <- dim(SIMU_DATA)[1]
  p <- dim(SIMU_DATA)[2]
  Lambda <- seq(from= MIN_LAMBDA, to = MAX_LAMBDA, by = STEP_LAMBDA)
  Total_Iter <- length(K)*length(Lambda)
  
  Results <- matrix(rep(NA, Total_Iter), ncol = length(Lambda) )
  Converge_status <- matrix(rep(NA, Total_Iter), ncol = length(Lambda) )
  Num_zero_mean <- matrix(rep(0, Total_Iter), ncol = length(Lambda) )
  BIC <- matrix(rep(NA, Total_Iter), ncol = length(Lambda) )
  ARI <- matrix(rep(NA, Total_Iter), ncol = length(Lambda) )
  ZERO_IND <- matrix(rep(NA, Total_Iter), ncol = length(Lambda) )
  

  # Check true cluster 
  if(all(is.na(True_cluster)))
  {
    print("No true cluster label provided")
    True_cluster <- rep(1,N)
  }
  
  
  for(k in 1:length(K))
  {
    for(lamb in 1:length(Lambda))
    {

      EM_result <- EM_simple(SIMU_DATA, K= K[k], Lambda = Lambda[lamb], Lambda2 = Lambda2, M = M, method = method, km = km )           
      
      # See if all means are 0
      if(EM_result$all_zero)
      {
        Results[k, lamb] <- 0
        Converge_status[k, lamb] <- M - (EM_result$Stop_id + 1) 
        Num_zero_mean[k, lamb] <- p*K[k]
        BIC[k, lamb] <- 0
        ZERO_IND[k, lamb] <- 0
        ARI[k, lamb] <- mclust::adjustedRandIndex(EM_result$Cluster , True_cluster)

      }
      else 
      {

        Results[k, lamb] <- as.numeric(EM_result$logLik)
        Converge_status[k, lamb] <- M - (EM_result$Stop_id + 1) # if Converge status is 0 then means has not converged
        for(m in 1:K[k]){Num_zero_mean[k, lamb] <- Num_zero_mean[k, lamb]+ length( which(EM_result$MU[[m]][EM_result$Stop_id,]==0)) }
        BIC[k, lamb] <- -2*Results[k, lamb] + log(N)*( p+K[k]+p*K[k]-1- Num_zero_mean[k, lamb] )
        ARI[k, lamb] <- mclust::adjustedRandIndex(EM_result$Cluster , True_cluster)

        ZERO_IND[k, lamb] <- sum(EM_result$zero_prob_ind) # Total number of 0 prob for N*p combinations
      }
      
    }
  }
  
  rownames(Results) <- c(paste(K))
  colnames(Results) <- c(paste(Lambda))
  
  rownames(Converge_status) <- c(paste(K))
  colnames(Converge_status) <- c(paste(Lambda))
  
  rownames(BIC) <- c(paste(K))
  colnames(BIC) <- c(paste(Lambda))
  
  rownames(Num_zero_mean) <- c(paste(K))
  colnames(Num_zero_mean) <- c(paste(Lambda))
  
  rownames(ZERO_IND) <- c(paste(K))
  colnames(ZERO_IND) <- c(paste(Lambda))
  
  rownames(ARI) <- c(paste(K))
  colnames(ARI) <- c(paste(Lambda))
  
  
  
  par(mfrow=c(1,2))
  
  if(plot==TRUE)
  {

    
    MIN <- min(  Results[ intersect(which(Results!=0) , which(!is.infinite(Results)) ) ] )
    MAX <- max(  Results[ intersect(which(Results!=0) , which(!is.infinite(Results)) ) ] )

    Results[which(is.infinite(Results))] <- 0
    
    plot(Lambda , Results[1,] , main="Log-likelihood for K and Lambda", type = "l",xlab = "Lambda", ylab = "Log-Likelihood", ylim = c(MIN,MAX))  
    for(i in 2:length(K))
    {
      lines(Lambda,Results[i,], col = i, lty = i)
    }
    legend("bottomleft", legend = c(paste("K=",c(K))), lty = c(1:length(K)), col = c(1:length(K)), cex = 0.5)
  }
  if(plot==FALSE){}
  
  
  if(plot==TRUE)
  {

    
    MIN <- min(  BIC[ intersect(which(BIC!=0) , which(!is.infinite(BIC)) ) ] )
    MAX <- max(  BIC[ intersect(which(BIC!=0) , which(!is.infinite(BIC)) ) ] )
    
    BIC[which(is.infinite(BIC))] <- 0
    plot(Lambda , BIC[1,] , main="Adjusted BIC for K and Lambda", type = "l",xlab = "Lambda", ylab = "Adjusted BIC", ylim = c(MIN,MAX))  
    for(i in 2:length(K))
    {
      lines(Lambda,BIC[i,], col = i, lty = i)
    }
    legend("topleft", legend = c(paste("K=",c(K))), lty = c(1:length(K)), col = c(1:length(K)), cex = 0.5)
  }
  if(plot==FALSE){}
  
  
  
  min_bic <- which( (BIC==min(BIC[BIC!=0]) ), arr.ind = TRUE)
  best_lambda <- Lambda[min_bic[2]]
  best_K <- K[min_bic[1]]
  
  par(mfrow=c(1,1))
  
  
  return(list(Result = Results, Converge = Converge_status, BIC = BIC,Num_zero_mean = Num_zero_mean, 
              Best_Lambda = best_lambda, Best_K = best_K , ZERO_IND = ZERO_IND, ARI = ARI, min_loc = as.numeric(min_bic)))
}





# This function is designated to generate L dimensional functional data with N samples of length P
# each dimension is stored in each compoenent of the result
# 
# L: dimensions of mfdata, or number of variables/sensors in the data
# P: numebr of measure points in the functional data
# N: numerb of samples that users would like to generate
# nbasis: number of basis used in function 'create.bspline.basis'
# norder: number of order used in function 'create.bspline.basis'
# correlations: correlations among the functional data, this should be a choose(L, 2) vector, if not specified, rep(0,choose(L,2)) will be used
# Coef_Matrix: Coefficients that are used when generating the functional data from B-spline
#              dimension of this matrix should be L*nbasis, each row represent the coefficient for each centerline.


Gaussian_mfdata <- function(L, P, N, nbasis = 20, norder = 4, a = 0.01, b = 10,
                            Coef_Matrix = NA,correlations = NA, plot = TRUE)
{

  alpha <- rep(a, L)
  beta <- rep(b, L)
  

  t <- 1:P
  

  splines <- create.bspline.basis(rangeval=c(1, max(t)), nbasis = nbasis ,norder=norder)
  
  
  time_grid = seq( min(t), max(t), length.out = P )
  for(i in 1:L)
  {
    assign(paste('C','_',i, sep="") ,exp_cov_function( time_grid, alpha = alpha[i], beta = beta[i] ) )
  }
  

  
  if(all(is.na(Coef_Matrix)))
  {  
    for(i in 1:L)
    {
      assign( paste("key",i, sep=""), fd(coef = as.vector( rmvnorm(1, mean = rep(0,nbasis), sigma = 2*diag(nbasis)) ) , basisobj = splines) )   
      assign(   paste("key",i,"_","centerline", sep="") , eval.fd(seq(1,P, 1), get(paste("key",i, sep=""))  )    )
    
    }
  }
  else
  {
    for(i in 1:L)
    {
      assign( paste("key",i, sep=""), fd(coef = Coef_Matrix[i,] , basisobj = splines) )   
      assign(   paste("key",i,"_","centerline", sep="") , eval.fd(seq(1,P, 1), get(paste("key",i, sep=""))  )    )
    }
  }
  centerline <- do.call(cbind,mget(paste("key",1:L,"_","centerline",sep="")))
  colnames(centerline) <- paste("key",1:L,"_","centerline",sep="")
  centerline <- t(centerline)
  

  if(is.na(correlations))
  {
    correlations <- rep(0, choose(L,2))
    
    
    mfdata <- generate_gauss_mfdata( N, L, centerline,
                                     correlations = correlations,
                                     listCov = mget(paste("C","_",1:L,sep="")) )
  }
  else
  {
    mfdata <- generate_gauss_mfdata( N, L, centerline,
                                     correlations = correlations,
                                     listCov = mget(paste("C","_",1:L,sep="")) )
  }
  
  if(plot==TRUE)
  {
    row <- ceiling(L/2)
    par(mfrow=c(row,2))
    for(i in 1:L)
    {
      matplot(t(mfdata[[i]]), lty=1,type="l",col=1, cex = 0.01, xlab = "X", ylab = "Y")
    }
    par(mfrow=c(1,1))
  }
  else{}
  
  return(mfdata = mfdata)
  
}




# This function is designate to convert multivariate functional data into matrix form which is compatible for our EM method
#
# 
# value: This function will return a matrix with sensory data observation with ID index. 
#        basically, it will return N submatrix stacked together the last column is ID that indicates the observation
#        it comes from. 

mfdata_to_matrix <- function(mfdata)
{
  N <- dim(mfdata[[1]])[1]
  P <- dim(mfdata[[1]])[2]
  L <- length(mfdata)
  
  # Create a matrix to restore the result
  Mul_Fun_Matrix <- matrix( rep(NA, N*P*(L+1)) , ncol = (L+1) )
  
  for(i in 1:N)
  {
    index <- ((i-1)*P+1) : (i*P)
    for(j in 1:L)
    {
      Mul_Fun_Matrix[index, j] <- mfdata[[j]][i,]
    }
    Mul_Fun_Matrix[index, L+1] <- rep(i, P)
  }
  
  colnames(Mul_Fun_Matrix) <- c(paste("Cencor", 1:L) , "ID")
  
  return(Mul_Fun_Matrix)
}




MFPCA_sensory <- function(Mul_Fun_Matrix , M)
{

  L <- dim(Mul_Fun_Matrix)[2]-1 # Number of sensors
  TT <- dim(Mul_Fun_Matrix[Mul_Fun_Matrix[,L+1]==1,])[1] 
  N <- length(unique(Mul_Fun_Matrix[,L+1])) # number of observations
  p <- L  # number of sensors
  t=1:TT
  
  splines <- create.bspline.basis(rangeval=c(1, max(t)), nbasis = 20,norder=3) 

  MFPCA_results <- matrix(rep(NA, N*p*M), ncol = p*M)
  
  #M eigenvalues for each obs
  EigenValues <- matrix(rep(NA, N*M), ncol = M)
  
  Eigenfunctions <- list()

  
  
  for( i in 1: L)
  {
    
    pca_data <- matrix(Mul_Fun_Matrix[,i], ncol = N)

    pca_data_fd <- Data2fd(y = pca_data, argvals=t, basisobj=splines) 

    

    mfpca<- pca.fd(pca_data_fd, nharm = M, centerfns = TRUE) 
    

    index <- ((i-1)*M+1 ) : (i*M)
    
    MFPCA_results[, index] <- mfpca$scores
  }
  
  colnames(MFPCA_results) <- rep( paste("sensor", 1:L) , each = M)
  
  
  return(MFPCA_results)
  
}



# Function Cluster_coef_mat_generator
# This function generates coefficients from MFPCA 
# L: number of variables/sensors 
# P: numebr of measure points
# N: number of samples for each cluster
# nbasis: number of basis used in B-spline generator
# norder: number of orders used in B-spline generator
# M: number of principle component uses in MFPCA
# ncluster: numebr of clusters uses woudl like to generate
# Coef_Matrix_list: a list of length 'ncluster', with each component to be a matrix as 'Coef_Matrix' define
#                   in the function 'Gaussian_mfdata'


Cluster_coef_mat_generator <- function(L, P, N, nbasis, norder, M, ncluster, plot = FALSE,
                                       a = 0.01, b = 10, Coef_Matrix_list = rep(list(NA),L) )
{
  MFPCA <- data.frame()
  mfdata_total <- data.frame()
  
  if(all(is.na(Coef_Matrix_list)))
  {
    Coef_Matrix_list <- rep(list(NA),ncluster)
  }
  else
  {}
  
  #browser()
  for(i in 1:ncluster)
  {
    mfdata <- Gaussian_mfdata(L, P, N, nbasis = nbasis, norder = norder, a = a, b = b, Coef_Matrix =Coef_Matrix_list[[i]], plot = plot)
    Mul_Fun_Matrix <- mfdata_to_matrix(mfdata)
    MFPCA <- rbind(MFPCA, cbind(MFPCA_sensory(Mul_Fun_Matrix, M) , rep(i, N)) )
    mfdata_total <- rbind.data.frame(mfdata_total, Mul_Fun_Matrix)
  }
  return(list(MFPCA = MFPCA, mfdata = mfdata_total) )
  
}





# Function to generate labels
# K: number of clusters
# N: number of observaitons
# prob: probability vector

label_generator <- function(N, K, prob)
{
  if(K != length(prob))
  {
    print("Error: Number of probabilities does not match the number of clusters!")
    break
  }
  label_result <- vector()
  labels_cluster <-  t( rmultinom(N, 1, prob) ) ##multinomial with sample size N and prob for K catagories
  for(i in 1:N)
  {
    label_result <- c( label_result, which(labels_cluster[i,]==1) )##transform from hot coding back to labels
  }
  
  return(label_result)
}



# Function to generate functional data
# Input
# N: number of observations
# K: number of clusters
# P: number of variables
# PC: number of basis we use
# order: order used in bassis creation
# evalRange: range of functional data
# MU: A list of length K, include mean vectors for each cluster, each mean vector is of length P*PC. 
# SIGMA: A list of length K, include the covariance matrice for each cluster, each mean vector is of length P*PC. 
# prob: probability for 
FDATA_GENERATOR <- function(N, K, P, PC, order, evalRange, MU, SIGMA, prob)
{
  #browser()
  num_cluster <- K
  num_sensor <- P
  num_pc <- PC
  num_eval <- length(seq(evalRange[1], evalRange[2], 1))
  num_order <- order
  
  labels_cluster <- label_generator(N, K, prob)
  
  
  for(k in 1:num_cluster)
  {
    assign(paste('mu_', k, sep = ""),  MU[[k]] )
    assign(paste('sigma_', k, sep = ""), SIGMA[[k]])
  }
  
  # Create basis function
  basis_fn <- create.bspline.basis(rangeval = evalRange, nbasis = num_pc, norder = num_order)

  # Coefficient Matrix
  coef_matrix <- matrix( rep(NA, N*num_sensor * num_pc), ncol = num_sensor * num_pc)
  
  # Data Matrix
  Eval_matrix <- matrix( rep(NA, N* num_eval*(num_sensor+1)), ncol = num_sensor+1)
  
  
  for(i in 1:N)
  {

    
    cluster <- labels_cluster[i]
    mu <- get(paste('mu_', cluster, sep = ""))
    sigma <- get(paste('sigma_', cluster, sep = ""))
    coef_matrix[i,] <- rmvnorm(1, mean = mu, sigma = sigma)
    
    
    # Create functional data with the basis for each sensor
    current_coef <- coef_matrix[i,]
    
    #browser()
    for(p in 1:num_sensor)
    {
      col_index <- ( (p-1)*num_pc +1 ) : ( (p-1)*num_pc + num_pc )
      row_index <- ( (i-1)*num_eval +1 ) : ( (i-1)*num_eval +num_eval )
      
      coef_val <- current_coef[col_index]
      fun_data <- fd(coef_val, basis_fn)
      
      # Obtain the distinct value
      fd_matrix <- eval.fd(seq(evalRange[1], evalRange[2], 1),fun_data)
      
      # Save functional evaluation in the matrix
      Eval_matrix[row_index,p] <- as.numeric( fd_matrix )
      
    }
    Eval_matrix[row_index,p+1] <- rep(i, num_eval)
    
  }

  
  return(list( Data = Eval_matrix, Coef = coef_matrix, Cluster = labels_cluster))
  
}

# Function: to return the optimal hyperparameters including lambda and gamma 
# Input
# Lambda_set: vector with lambda parameters to choose from
# gamma_set: vector with gamma parameters to choose from
# Data: Data to conduct clustering with
# K: number of clusters we use
# penalty_method: one of 'individual', 'variable' and 'grouped'.
# km: Only used in 'grouped' penalty, indicating the length of variable groups.
# min_percentage: Minimum clustering percentage for initialzation
# plot: if TRUE, plot the contour plot

Hyper_Parameter_Selection <- function(Lambda_set, gamma_set,Data, K, penalty_method = "individual", km, min_percentage, PCA_exp=3,plot = FALSE )
{

  HP_result_matrix_BIC <- matrix(rep(NA, length(Lambda_set) * length(gamma_set)), ncol = length(gamma_set) )
  HP_result_matrix_Lik <- matrix(rep(NA, length(Lambda_set) * length(gamma_set)), ncol = length(gamma_set) )
  

  result_list_BIC <- rep(list(HP_result_matrix_BIC),length(K))
  result_list_Lik <- rep(list(HP_result_matrix_Lik),length(K))
  
  MU_tilda_list<-vector(mode="list",length=length(K))
  total_iteration <- length(Lambda_set)*length(gamma_set)*length(K)
  
  for(k in 1:length(K)){
    
    tmpBIC<-rep(NA,length(Lambda_set))
  
    
    for( i in 1:length(Lambda_set)){
      Count <- 0
      
      repeat{
      
      tmp <- EM_simple_first_step(DATA=Data, K[k], Lambda=Lambda_set[i], Lambda2=10, sigma_method = "equal", method = penalty_method,  M=150,
                                  km = km, min_percentage = 0.1,PCA_exp=PCA_exp)
      
      
      if(min(tmp$PI[tmp$Stop_id,]) > 0)
      {
        break()
      }
      
      if(Count>0){print ("restarted ")}
      
      if(Count > 5){ 
        print("restarted 5 times, break")
        paste("K",K[k],sep="")
        print("Lambda",Lambda_set[i],sep="")
        break()
      }
      
      Count <- Count +1
      
    
      
      }

      
      if(Count >5){
        tmpBIC[i]<-Inf
      }else{
        
        tmpBIC[i]<-as.numeric(tmp$BIC)
      }
      
      }
      
    
    
    idbic<-which.min((tmpBIC))
    
    
    tmp <- EM_simple_first_step(DATA=Data, K[k], Lambda_set[idbic], Lambda2=10, sigma_method = "equal", method = penalty_method,  M=200,
                                km = km, min_percentage = 0.1,PCA_exp=PCA_exp)
    tmpstopid <- tmp$Stop_id
    MU_tilda <- vector(mode="list",length=K[k])
    
    
    for( i in 1:K[k]){
      
      MU_tilda[[i]]<-tmp$MU[[i]][tmpstopid,]
    }
    MU_tilda_list[[k]]<-MU_tilda

    for(j in 1:length(gamma_set)){
      
      
      for(i in 1:length(Lambda_set))
      {
        
        Count <- 0
        
        repeat{
          result <- EM_adaptive(MU_tilda=MU_tilda, Data, K=K[k], Lambda=Lambda_set[i], Lambda2 = 10,
                                         sigma_method = 'equal', method = penalty_method,
                                         M=100, gamma = gamma_set[j], km = km, min_percentage = min_percentage,PCA_exp=PCA_exp)
       
          
          if(min(result$PI[result$Stop_id,]) > 0)
          {
            break()
          }
          

          
          if(Count>0){print ("restarted ")}
          
          if(Count > 5 ){ 
            print("restarted 5 times, break")
            print (Lambda_set[i])
            print (gamma_set[j])
            print("K")
            print(K[k])
            
            break()
          }
          
          Count <- Count +1
          

          
          
          
        }

        if(Count >5 ){  #|| checkpoint==1
          result_list_BIC[[k]][i,j] <- Inf
          result_list_Lik[[k]][i,j] <- as.numeric(result$logLik)
          
        }else{
        
        result_list_BIC[[k]][i,j] <- as.numeric(result$BIC)
        result_list_Lik[[k]][i,j] <- as.numeric(result$logLik)
        }
      }
      
    }
  }

  MIN_BIC_vec <- unlist(lapply(result_list_BIC, min))
  min_index <- which(unlist(lapply(result_list_BIC, min)) == min(unlist(lapply(result_list_BIC, min))))
  
  best_MU_tilda<-MU_tilda_list[[min_index]]
  

  HP_result_matrix_BIC <- result_list_BIC[[min_index]]
  HP_result_matrix_Lik <- result_list_Lik[[min_index]]
  
  colnames(HP_result_matrix_BIC) <- gamma_set
  rownames(HP_result_matrix_BIC) <- Lambda_set
  
  colnames(HP_result_matrix_Lik) <- gamma_set
  rownames(HP_result_matrix_Lik) <- Lambda_set
  
  if(plot == TRUE)
  {
    par(mfrow = c(1,1))
    contour(Lambda_set, gamma_set, HP_result_matrix_BIC, xlab = expression(lambda), ylab = expression(gamma), main = "BIC on Engineer Data")
    contour(Lambda_set, gamma_set, HP_result_matrix_Lik, xlab = expression(lambda), ylab = expression(gamma), main = "Log-Likelihood on Engineer Data")
  }
  

  row <- dim(which(HP_result_matrix_BIC == min(HP_result_matrix_BIC), arr.in=TRUE))[1]
  
  row_index <- which(HP_result_matrix_BIC == min(HP_result_matrix_BIC), arr.in=TRUE)[row,1]
  col_index <- which(HP_result_matrix_BIC == min(HP_result_matrix_BIC), arr.in=TRUE)[row,2]
  
  best_lambda <- Lambda_set[row_index]
  best_gamma <- gamma_set[col_index]
  
  Min_BIC <- min(HP_result_matrix_BIC)
  

  
  return(list(MU_tilda = best_MU_tilda, lambda = best_lambda, gamma = best_gamma, K = K[min_index],
              min_BIC_vec = MIN_BIC_vec,
              BIC_matrix  = result_list_BIC,
              loglike_matrix  = result_list_Lik))
  
}





# Function that combines the simulation results
combine_simu_data_cluster <- function(File_names)
{

  Result_Data <- data.frame()
  for(i in 1:length(File_names) )
  {
    print(i)
    new_data <- read.table(File_names[i], header = T, sep = ",", stringsAsFactors = F)
    

    if(all(!is.na(new_data)))
    {

      Result_Data <- rbind.data.frame(Result_Data, new_data)
      
      
    }
    else{
      print(i)
      print(File_names[i])
      next
    }

  }
  
  return(Result_Data)
}


MSE_cal <- function(Tables = Tables, LABELS = LABELS, var = "Mu1_c", True_val, abs.val = FALSE)
{
  if(length(True_val)==1)
  {
    True_val <- rep(True_val, length(LABELS))
  }
  

  result_vec <- rep(NA, length(LABELS))
  for(i in 1:length(LABELS))
  {
    index <- which(Tables[,c("label")]== LABELS[i] )
    result_vec[i] <- mean((Tables[index,var]- True_val[i] )^2)
  }
  
  if(abs.val == TRUE)
  {
    for(i in 1:length(LABELS))
    {
      index <- which(Tables[,c("label")]== LABELS[i] )
      result_vec[i] <- mean((abs(Tables[index,var]- True_val )^2))
    }
  }
  return(result_vec)
}


# Function to check the number of sensors removed by our models 
sensor_remove_check <- function(Table_result, km, true_val, sensor_val, col_index)
{

  result_vec <- rep(NA, dim(Table_result)[1])
  False_vec <- rep(NA, dim(Table_result)[1])
  
  for(i in 1:dim(Table_result)[1])
  {
    if(is.na(Table_result[i,col_index]))
    {
      result_vec[i] <- 0
      False_vec[i] <- 0
      next
    }
    if(as.character(Table_result[i,col_index]) == "Not Removed")
    {
      result_vec[i] <- 0
      False_vec[i] <- 0
      next
    }
    var_removed <- as.numeric(unlist(strsplit(as.character(Table_result[i,col_index]), ",")))
    count_removed <- 0
    for(m in 1:(length(true_val)/km))
    {
 
      index <- ((m-1)*km +1) : ((m-1)*km +km)
      if(all(true_val[index] %in% var_removed))
      {
        count_removed <- count_removed + 1
      }
    }
    result_vec[i] <- count_removed
    
    
    count_removed_fp <- 0
    for(n in 1:(length(sensor_val)/km))
    {
      index <- ((n-1)*km +1) : ((n-1)*km +km)
      if(all(sensor_val[index] %in% var_removed))
      {
        count_removed_fp <- count_removed_fp + 1
      }
    }
    False_vec[i] <- count_removed_fp

  }
  return(list(result_vec = result_vec, False_vec = False_vec))
}



MFPCA_sensory_check <- function(Mul_Fun_Matrix , M)
{

  L <- dim(Mul_Fun_Matrix)[2]-1 # Number of sensors
  TT <- dim(Mul_Fun_Matrix[Mul_Fun_Matrix[,L+1]==1,])[1] 
  N <- length(unique(Mul_Fun_Matrix[,L+1])) # number of observations
  p <- L  # number of sensors
  t=1:TT
  
  splines <- create.bspline.basis(rangeval=c(1, max(t)), nbasis = 60,norder=4) 

  MFPCA_varprop <- matrix(rep(NA, L*M), ncol = M)
  
  
  for( i in 1: L)
  {
    
    pca_data <- matrix(Mul_Fun_Matrix[,i], ncol = N)

    pca_data_fd <- Data2fd(y = pca_data, argvals=t, basisobj=splines) 
    
    mfpca<- pca.fd(pca_data_fd, nharm = M, centerfns = TRUE) 
    
    MFPCA_varprop[i, 1:M] <- mfpca$varprop
  }
  
  rownames(MFPCA_varprop) <- paste("sensor", 1:L)
  
  
  return(MFPCA_varprop)
  
}

##########################REAL dat##################
real.onerun<-function(sf.id, data=FALSE, num_PCA=3){
  
  if(data==TRUE){

    b <- read.csv("System_B.csv", head = T, sep=",")
  

  
  b_rescale<-cbind(apply(b[,1:66],2,scale),b[,67])
  colnames(b_rescale)[67]<-"Stop_ID"
  

  b_var<-apply(b[,1:66],2,var)
  b_mean <- apply(b[,1:66],2,mean)
  
  idx<-which(b_var==0)
  for( i in 1:length(idx)){
    b_rescale[,idx[i]]<- b[,idx[i]]-b_mean[idx[i]]
  }
  
  
  
  data_pre <- MFPCA_sensory(b_rescale, M = 2) 
  
  rm.var.index<-which(apply(data_pre, 2, sd)==0) #29,30,31,32,49,50,97,98,99,100
  rm.var.index<-sort(c(rm.var.index,c(7,8,9,10,17,18,21,22,23,24,25,26,29,30,33,34,37,38,39,40,41,42,
                                      43,44,47,48,55,56,61,62,101,102, 107,108,109,110,111,112)))
  rm.sensor.index<-round(rm.var.index[seq(1,length(rm.var.index),by=2)+1]/2)
  
  b_rescale_new<-b_rescale[,-rm.sensor.index]
  
  
  car_data_new<-MFPCA_sensory(b_rescale_new,M=num_PCA)}else if(data==FALSE){
    
    car_data_new<-readRDS('sysBpca.rds')
    
  }

  
  
  
  ####################
  # Hyper parameter Selection
  ####################
  N<-nrow(car_data_new)
  Lambda_set<-c(0,0.5,1,1.5, 2,3,5,7,10, 15,20, 25)*(N^(1/3))
  gamma_set<-c(0.5,1,1.5,2)
  

  
  
  hyper_individual <- Hyper_Parameter_Selection(Lambda_set, gamma_set,Data = car_data_new, K = c(1,2,3,4), 
                                                penalty_method = "individual_equal", km = num_PCA,PCA_exp=num_PCA, 
                                                min_percentage = 0.03, plot = FALSE)
  
  write.csv(hyper_individual$lambda,paste(sf.id,"hyperindividual_lambda.csv",sep=""))
  write.csv(hyper_individual$gamma,paste(sf.id,"hyperindividual_gamma.csv",sep=""))
  write.csv(hyper_individual$K,paste(sf.id,"hyperindividual_K.csv",sep=""))
  
  hyper_variable <- Hyper_Parameter_Selection(Lambda_set, gamma_set,Data = car_data_new, K = c(1,2,3,4), 
                                              penalty_method = "variable", km = num_PCA, PCA_exp= num_PCA, 
                                              min_percentage = 0.03, plot = FALSE )
  
  
  write.csv(hyper_variable$lambda,paste(sf.id,"hypervariable_lambda.csv",sep=""))
  write.csv(hyper_variable$gamma,paste(sf.id,"hypervariable_gamma.csv",sep=""))
  write.csv(hyper_variable$K,paste(sf.id,"hypervariable_K.csv",sep=""))
  
  
  
  hyper_grouped <- Hyper_Parameter_Selection(Lambda_set, gamma_set,Data = car_data_new, K = c(1,2,3,4), 
                                             penalty_method = "grouped", km = num_PCA,PCA_exp=num_PCA, 
                                             min_percentage = 0.03, plot = FALSE )
  
  
  
  write.csv(hyper_grouped$lambda,paste(sf.id,"hypergrouped_lambda.csv",sep=""))
  write.csv(hyper_grouped$gamma,paste(sf.id,"hypergrouped_gamma.csv",sep=""))
  write.csv(hyper_grouped$K,paste(sf.id,"hypergrouped_K.csv",sep=""))
  
  
  
  
  EM_result_individual <- EM_adaptive(MU_tilda = hyper_individual$MU_tilda, DATA = car_data_new, 
                                      K=hyper_individual$K, Lambda= hyper_individual$lambda ,
                                      Lambda2 = 10, sigma_method = 'equal', method = "individual_equal",
                                      M=200,gamma = hyper_individual$gamma, km = num_PCA,PCA_exp=num_PCA,
                                      min_percentage = 0.05)
  
  
  EM_result_variable<- EM_adaptive(MU_tilda = hyper_variable$MU_tilda, car_data_new, K=hyper_variable$K, 
                                   Lambda= hyper_variable$lambda, Lambda2 = 10, sigma_method = 'equal', 
                                   method = "variable", M=200, gamma = hyper_variable$gamma, km = num_PCA,PCA_exp=num_PCA,
                                   min_percentage = 0.05)
  
  
  
  
  EM_result_grouped  <- EM_adaptive(MU_tilda = hyper_grouped$MU_tilda, car_data_new, 
                                    K=hyper_grouped$K, Lambda= hyper_grouped$lambda, Lambda2 = 10, 
                                    sigma_method = 'equal', method = "grouped", M=200,
                                    gamma = hyper_grouped$gamma, km = num_PCA, min_percentage = 0.05,PCA_exp=num_PCA )
  
  
  
  
  
  saveRDS(EM_result_individual, file = paste(sf.id,"EM_result_individual_pca.RDS",sep=""))
  saveRDS(EM_result_variable, file = paste(sf.id,"EM_result_variable_pca.RDS",sep=""))
  saveRDS(EM_result_grouped, file = paste(sf.id,"EM_result_grouped_pca.RDS",sep=""))
  
}
#########FIGURE2#################
figure2.plot<-function(x){
  N <- 200 
  colseq = rep(2:8,length.out = N)
  
  if(x=='a'){
    
    
    N <- 200 
    num_eval <- 30 
    
    
    mu1 <- c(c(14,14,rep(14,9),0), c(11,20,22,rep(23,5),seq(23,18,length.out=4)) ,rep(rep(0,12),2)) 
    mu2 <- c(c(11.5,12.5,13,seq(13,10,length.out=4),seq(9,4,length.out=4),0),
             c(10,15,13,seq(13,19,length.out=3),19,seq(20,23,length.out=3),20,18),
             rep(rep(0,12),2))
    mu3 <- c(c(5,seq(3,6.5,length.out=2),seq(6.5,3,length.out=2),seq(3,6.5,length.out=2),seq(6.5,3,length.out=2),seq(4,0,length.out=3)),
             c(rep(18,12)),
             rep(rep(0,12),2))
    
    pgseq<-2.5
    
    sigma1 <- diag(c(c(20,10,8,rep(2.8,6),rep(2.8,2),c(0.2)), 
                     c(rep(16,11),1), rep(15,12*2) ) )/pgseq  
    sigma2 <-diag(c(c(15,15,rep(10,5),seq(5,0.1,length.out=5)),
                    c(rep(30,5),rep(20,6),1),rep(15,12*2) ))     /pgseq
    sigma3 <- diag( c(c(15,rep(12,2),rep(12,2),rep(12,2),rep(12,2),seq(3,0.2,length.out=3)),
                      c(rep(8,12)),rep(15,12*2)) )
    
    MU <- list(mu1, mu2, mu3)
    SIGMA <- list(sigma1, sigma2, sigma3)
    
    
    Data_result <- FDATA_GENERATOR(N = N, K =3, P = 4, PC = 12, order= 2, evalRange = c(0,30), MU = MU, SIGMA = SIGMA, prob = c(1/3,1/3,1/3) )
    Data_rescale<-cbind(apply(Data_result$Data[,1:4],2,scale),Data_result$Data[,5])
    
    data_matrix<- Data_result$Data
    
    sensor1<-data_matrix
    sensor1.matrix<-data.frame(matrix(NA,ncol=31,nrow=N) )
    sensor1.matrix$cluster<-Data_result$Cluster
    for(i in 1:N){
      sensor1.matrix[i,1:31]<-sensor1[which(sensor1[,5]==i),1]
    }
    par(mai=c(0.85, 0.85, .05, .05))
    n<-nrow(t(sensor1.matrix[,1:31]))
    matplot(t(sensor1.matrix[,1:31]),type='l',xlab="Measure Points",ylab="Measurements",lty=1,col=colseq,cex.axis=1.5,cex.lab=1.3)#col=sensor1.matrix[,32],
    m1<-sensor1.matrix[which(sensor1.matrix[,32]==1),1:31]
    m1.mean<-apply(m1,2,mean)
    lines(1:31,m1.mean,lwd=2)
    m2<-sensor1.matrix[which(sensor1.matrix[,32]==2),1:31]
    m2.mean<-apply(m2,2,mean)
    lines(1:31,m2.mean,lwd=2)
    m3<-sensor1.matrix[which(sensor1.matrix[,32]==3),1:31]
    m3.mean<-apply(m3,2,mean)
    lines(1:31,m3.mean,lwd=2)
    
    
    
    
    
  }
  if(x=='b'){
    N <- 200 
    num_eval <- 30 
    
    mu1 <- c(c(14,14,rep(14,9),0), c(11,20,22,rep(23,5),seq(23,18,length.out=4)) ,rep(rep(0,12),2)) 
    mu2 <- c(c(11.5,12.5,13,seq(13,10,length.out=4),seq(9,4,length.out=4),0),
             c(10,15,13,seq(13,19,length.out=3),19,seq(20,23,length.out=3),20,18),
             rep(rep(0,12),2))
    mu3 <- c(c(5,seq(3,6.5,length.out=2),seq(6.5,3,length.out=2),seq(3,6.5,length.out=2),seq(6.5,3,length.out=2),seq(4,0,length.out=3)),
             c(rep(18,12)),
             rep(rep(0,12),2))
    
    
    
    pgseq<-1
    
    
    sigma1 <- diag(c(c(20,10,8,rep(2.8,6),rep(2.8,2),c(0.2)), 
                     c(rep(16,11),1), rep(15,12*2) ) )/pgseq  #Variance, corresponding
    
    sigma2 <-diag(c(c(15,15,rep(10,5),seq(5,0.1,length.out=5)),
                    c(rep(30,5),rep(20,6),1),rep(15,12*2) ))     /pgseq
    
    sigma3 <- diag( c(c(15,rep(12,2),rep(12,2),rep(12,2),rep(12,2),seq(3,0.2,length.out=3)),
                      c(rep(8,12)),rep(15,12*2)) )
    
    MU <- list(mu1, mu2, mu3)
    SIGMA <- list(sigma1, sigma2, sigma3)
    
    
    Data_result <- FDATA_GENERATOR(N = N, K =3, P = 4, PC = 12, order= 2, evalRange = c(0,30), MU = MU, SIGMA = SIGMA, prob = c(1/3,1/3,1/3) )
    Data_rescale<-cbind(apply(Data_result$Data[,1:4],2,scale),Data_result$Data[,5])
    
    data_matrix<- Data_result$Data
    
  
    sensor1<-data_matrix
    sensor1.matrix<-data.frame(matrix(NA,ncol=31,nrow=N) )
    sensor1.matrix$cluster<-Data_result$Cluster##in the explanation of real data, which means 50 stops, each obs period is 31
    for(i in 1:N){
      
      sensor1.matrix[i,1:31]<-sensor1[which(sensor1[,5]==i),1]
      
    }
    
    colseq = rep(2:8,length.out = N)
    par(mai=c(0.85, 0.85, .05, .05))
    n<-nrow(t(sensor1.matrix[,1:31]))
    matplot(t(sensor1.matrix[,1:31]),type='l',xlab="Measure Points",ylab="Measurements",lty=1,col=colseq,cex.axis = 1.5,cex.lab = 1.3) #,col=1
    m1<-sensor1.matrix[which(sensor1.matrix[,32]==1),1:31]
    m1.mean<-apply(m1,2,mean)
    lines(1:31,m1.mean,lwd=2)
    m2<-sensor1.matrix[which(sensor1.matrix[,32]==2),1:31]
    m2.mean<-apply(m2,2,mean)
    lines(1:31,m2.mean,lwd=2)
    m3<-sensor1.matrix[which(sensor1.matrix[,32]==3),1:31]
    m3.mean<-apply(m3,2,mean)
    lines(1:31,m3.mean,lwd=2)
    
  }
  
  if(x=='c'){
    
    N <- 200 
    num_eval <- 30 
    
   
    mu1 <- c(c(14,14,rep(14,9),0), c(11,20,22,rep(23,5),seq(23,18,length.out=4)) ,rep(rep(0,12),2)) 
    mu2 <- c(c(11.5,12.5,13,seq(13,10,length.out=4),seq(9,4,length.out=4),0),
             c(10,15,13,seq(13,19,length.out=3),19,seq(20,23,length.out=3),20,18),
             rep(rep(0,12),2))
    mu3 <- c(c(5,seq(3,6.5,length.out=2),seq(6.5,3,length.out=2),seq(3,6.5,length.out=2),seq(6.5,3,length.out=2),seq(4,0,length.out=3)),
             c(rep(18,12)),
             rep(rep(0,12),2))
    
    
    
    pgseq<-1
    
    
    sigma1 <- diag(c(c(20,10,8,rep(2.8,6),rep(2.8,2),c(0.2)), 
                     c(rep(16,11),1), rep(15,12*2) ) )/pgseq  #Variance, corresponding
    
    sigma2 <-diag(c(c(15,15,rep(10,5),seq(5,0.1,length.out=5)),
                    c(rep(30,5),rep(20,6),1),rep(15,12*2) ))     /pgseq
    
    sigma3 <- diag( c(c(15,rep(12,2),rep(12,2),rep(12,2),rep(12,2),seq(3,0.2,length.out=3)),
                      c(rep(8,12)),rep(15,12*2)) )
    
    MU <- list(mu1, mu2, mu3)
    SIGMA <- list(sigma1, sigma2, sigma3)
    
    
    Data_result <- FDATA_GENERATOR(N = N, K =3, P = 4, PC = 12, order= 2, evalRange = c(0,30), MU = MU, SIGMA = SIGMA, prob = c(1/3,1/3,1/3) )
    Data_rescale<-cbind(apply(Data_result$Data[,1:4],2,scale),Data_result$Data[,5])
    
    data_matrix<- Data_result$Data
    
    
    sensor1<-data_matrix
    sensor1.matrix<-data.frame(matrix(NA,ncol=31,nrow=N) )
    sensor1.matrix$cluster<-Data_result$Cluster##in the explanation of real data, which means 50 stops, each obs period is 31
    for(i in 1:N){
      
      sensor1.matrix[i,1:31]<-sensor1[which(sensor1[,5]==i),1]
      
    }
    
    sensor3.matrix<-data.frame(matrix(NA,ncol=31,nrow=N) )
    sensor3.matrix$cluster<-Data_result$Cluster
    for(i in 1:N){
      
      sensor3.matrix[i,1:31]<-sensor1[which(sensor1[,5]==i),3]
      
    }
    
    
    par(mai=c(0.85, 0.85, .05, .05))
    n<-nrow(t(sensor3.matrix[,1:31]))
    matplot(t(sensor3.matrix[,1:31]),type='l',xlab="Measure Points",ylab="Measurements",lty=1,col=colseq,cex.axis=1.5, cex.lab=1.3)#,col=sensor3.matrix[,32]
    m1<-sensor3.matrix[which(sensor3.matrix[,32]==1),1:31]
    m1.mean<-apply(m1,2,mean)
    lines(1:31,m1.mean,lwd=2)
    m2<-sensor3.matrix[which(sensor3.matrix[,32]==2),1:31]
    m2.mean<-apply(m2,2,mean)
    lines(1:31,m2.mean,lwd=2)
    m3<-sensor3.matrix[which(sensor3.matrix[,32]==3),1:31]
    m3.mean<-apply(m3,2,mean)
    lines(1:31,m3.mean,lwd=2)
    
    
  }
  
  if( !(x %in% c('a','b','c'))) "incorrect argument"
  

  
}

