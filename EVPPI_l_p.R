# REEEVR - Assessing reliability limits of regression and approximation methods for VoI estimation 
# Develop a generic Markov model that can have any number of states and input parameters with varying level of correlation. 
# Qian Xin October 2023

# MLMC level l estimator
EVPPI_l_p<-function(l,N, EVPPI_std_p = NULL)
{
  if(is.null(EVPPI_std_p)) return("EVPPI_std_p must be supplied")
  
  # N is the number of outer samples, M=2^l is the number of inner samples
  sum1 <- rep(0, 7) # 7 is the number of level? why 7?
  # Need to store results of different stats
  M = 2^(l + 1)
  Np = max(M, 128) # what's the meaning of Np?
  
  inputs = 1:ceiling(M * N / Np)
  
  # Iterate over the inputs
  # Can add parallel computation here
  results <- foreach(i=inputs,.export=c("n_samples", "n_tx",
                                        "generate_net_benefit", "generate_input_parameters"),.packages='MASS') %dorng%
    {
      NN=min(Np, N*M-(i-1)*Np)
      EVPPI_std_p(M,NN/M)
    }  
  #NetB <- EVPPI_x_std_p(M, N)
  
  # Convert results of foreach to NetB matrix
  # If not interested in parallelisation could merge the foreach and for loop to simplify
  NetB = matrix(NA, M*N, n_tx)
  for(i in inputs){
    NN <- min(Np, N*M-(i-1)*Np)
    nn = min(i*Np,N*M)
    NetB[(nn-NN+1):nn,] = matrix(unlist(results[i]),NN,n_tx)
  }
  # Net benefit based on partial perfect information for every sample
  NetB_max = apply(NetB,1,max)
  # Expected value of max over each set of inner samples
  NetB_max_sample = apply(matrix(NetB_max,N,M,byrow=TRUE),1,mean)
  # Matrices needed for antithetic variable variance reduction
  NetB_low_c_1 = matrix(NA,N,n_tx)
  NetB_low_c_2 = matrix(NA,N,n_tx)
  NetB_low_f = matrix(NA,N,n_tx)
  for(n in 1:n_tx){
    # Net benefIts for treatment n in matrix of size outer by inner samples
    temp = matrix(NetB[,n],N,M,byrow=TRUE)
    # Average net benefit over each set of inner samples
    NetB_low_f[,n] = apply(temp,1,mean)
    # Antithetic variable construction splits samples into first and second halves
    # Split formula into cases l=0 and l>0
    if(M==2){
      NetB_low_c_1[,n] = temp[,1]
      NetB_low_c_2[,n] = temp[,2]
    }else{
      if(N==1){
        NetB_low_c_1[,n] = sum(temp[,1:max(M/2,2)])/max(M/2,2)
        NetB_low_c_2[,n] = sum(temp[,min(M/2+1,M-1):M])/(M-min(M/2+1,M-1)+1)  
      }else{
        NetB_low_c_1[,n] = rowSums(temp[,1:max(M/2,2)])/max(M/2,2)
        NetB_low_c_2[,n] = rowSums(temp[,min(M/2+1,M-1):M])/(M-min(M/2+1,M-1)+1)        
      }
      
    }
  }
  NetB_low_f_sample = apply(NetB_low_f,1,max)
  # Put antithetic variable construction together
  NetB_low_c_sample = (apply(NetB_low_c_1,1,max)+apply(NetB_low_c_2,1,max))/2
  # Fine estimator (i.e. e_l^(n))
  Pf = NetB_max_sample - NetB_low_f_sample
  # Coarse estimator (i.e. e_(l-1)^(n))
  Pc = NetB_max_sample - NetB_low_c_sample
  
  # Sum the moments of the estimator
  # First is the mean of the difference estimator d_l^(n)
  # Summing these difference estimates gives an estimate of DIFF= EVPI-EVPPI
  sum1[1] = sum1[1] + sum(Pf-Pc);
  sum1[2] = sum1[2] + sum((Pf-Pc)^2);
  sum1[3] = sum1[3] + sum((Pf-Pc)^3);
  sum1[4] = sum1[4] + sum((Pf-Pc)^4);
  sum1[5] = sum1[5] + sum(Pf);
  sum1[6] = sum1[6] + sum(Pf^2);
  sum1[7] = sum1[7] + M*N;
  
  return(list(sums = sum1[1:6], cost = sum1[7]))
}
