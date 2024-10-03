# REEEVR - Assessing reliability limits of regression and approximation methods for VoI estimation 
# Develop a generic Markov model that can have any number of states and input parameters with varying level of correlation. 
# Qian Xin October 2023

EVPPI_MC <- function(M, N,treatment_costs_df, state_costs_df, state_utility_df, state_transition_matrices, 
                     uncertain_level = 100,  n_tx, n_cycles, n_states,hold_constant=c(), threshold = 25000){
  NetB <- array(0, dim = c(N, n_tx),
                dimnames = list(NULL, tx_names))

  
  for (i in 1:N){
    input_parameters <-
      generate_input_parameters(
        n_samples = M, 
        treatment_costs_df, 
        state_costs_df, 
        state_utility_df,
        state_transition_matrices, 
        hold_constant = hold_constant, 
        uncertain_level = 100
      )
    results <-
      generate_net_benefit(
        input_parameters = input_parameters,
        n_samples = M,
        n_tx, 
        n_cycles, 
        n_states, 
        threshold = 25000
      )
    NetB[i, ] <- apply(results$net_benefit, 1, mean)
    # NetB2[i] <- t(apply(results$net_benefit, 2, max))
  }
  
  Net_B_max <- apply(NetB,1,max)
  max_mean <- max(apply(NetB,2,mean))
  
  Diff <- Net_B_max - max_mean
  MCSE <- sd(Diff)/sqrt(N)
  
  EVPPI <- mean(Diff)
  return(data.frame(EVPPI, MCSE))
}
