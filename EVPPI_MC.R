# REEEVR - Assessing reliability limits of regression and approximation methods for VoI estimation 
# Develop a generic Markov model that can have any number of states and input parameters with varying level of correlation. 
# Qian Xin October 2023

library(foreach)
library(doParallel)
library(parallel)

EVPPI_MC <- function(M, N, treatment_costs_df, state_costs_df, state_utility_df, state_transition_matrices,
                     transition_names, state_names, tx_names, n_columns, 
                     uncertain_level = 100, n_tx, n_cycles, n_states,
                     hold_constant = c(), threshold = 25000) {
  
  # Register parallel backend
  num_cores <- detectCores() - 4  # Use all but one core
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  clusterExport(cl, c("rdirichlet", 'rdiric', 'bcea', 'createInputs', 'evppi'))
  clusterExport(cl, c('genCorData','melt','setkey','dcast'))
  
  # Use foreach to parallelize over N iterations
  NetB <- foreach(i = 1:N, .combine = rbind, .packages = c('MASS')) %dopar% {
    source("./generate_base_case_values.R")
    source("./generate_net_benefit.R")
    source("./generate_input_parameters.R")
    source("./generate_transition_matrix.R")
    
    # Generate input parameters
    input_parameters <- generate_input_parameters(
      n_samples = M, 
      n_states = n_states,
      n_tx =n_tx,
      treatment_costs_df, 
      state_costs_df, 
      state_utility_df,
      state_transition_matrices,
      transition_names,
      state_names,
      tx_names,
      n_columns,	
      hold_constant = hold_constant,
      uncertain_level = uncertain_level
    )
    
    # Generate net benefit
    results <- generate_net_benefit(
      input_parameters = input_parameters,
      n_samples = M,
      n_tx = n_tx,
      n_cycles = n_cycles,
      n_states = n_states,
      tx_names = tx_names,
      threshold = threshold
    )
    
    # Return mean net benefit for this iteration
    apply(results$net_benefit, 1, mean)
  }
  
  # Stop the parallel backend
  stopCluster(cl)
  
  # Compute the maximum net benefit for each sample and the maximum mean
  Net_B_max <- apply(NetB, 1, max)
  max_mean <- max(apply(NetB, 2, mean))
  
  # Calculate the difference and standard error
  Diff <- Net_B_max - max_mean
  MCSE <- sd(Diff) / sqrt(N)
  
  # Expected Value of Perfect Parameter Information (EVPPI)
  EVPPI <- mean(Diff)
  
  return(data.frame(EVPPI, MCSE))
}
