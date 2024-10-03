#

generate_net_benefit <- function(input_parameters, n_tx, n_samples, n_cycles, n_states, threshold = 25000){
  
  transition_matrices <- generate_transition_matrices(input_parameters, n_samples, n_tx, n_states)
  #transition_matrices[1,2,,]
  # test
  #matrix_to_test <- transition_matrices[1, 1, , ]
  #row_sums <- apply(matrix_to_test, 1, sum)
  #row_sums
  t1 = Sys.time()
  
  # generate the state costs for each sample, each state
  start_col <- n_tx*n_states*(n_states-1)+n_states
  end_col <-n_tx*n_states*(n_states-1)+2*n_states-2
  subset_df <- input_parameters[, start_col:end_col]
  column_name <- paste("costs_state_", n_states, sep = "") # last column/absorbing state
  subset_df[[column_name]] <- rep(0, nrow(subset_df))
  state_costs <- as.array(as.matrix(subset_df))
  t2 = Sys.time()
  
  # generate the treatment costs for each sample
  start_col <- n_tx*n_states*(n_states-1)+2*n_states-1
  end_col <-n_tx*n_states*(n_states-1)+2*n_states-1 +n_tx-1
  subset_df <- input_parameters[, start_col:end_col]
  tx_costs <- as.array(as.matrix(subset_df))
  # print(tx_costs[1,])
  t3 =Sys.time()
  # generate the state QALY for each sample, each state
  start_col <- n_tx*n_states*(n_states-1)+1
  end_col <-n_tx*n_states*(n_states-1)+n_states-1
  subset_df <- input_parameters[, start_col:end_col]
  column_name <- paste("costs_state_", n_states, sep = "")
  subset_df[[column_name]] <- rep(0, nrow(subset_df))
  state_qalys <- as.array(as.matrix(subset_df))
  t4 =Sys.time()
  # There is one cohort vector for each treatment, for each sample, for each cycle
  cohort_vectors <- array(dim = c(n_tx, n_samples, n_cycles, n_states),
                          dimnames = list(NULL, NULL, NULL, NULL))
  # Assume that everyone starts in state 1 no matter the treatment
  cohort_vectors[, , 1, ] <- 0
  cohort_vectors[, , 1, 1] <- 1
  
  #print(cohort_vectors[1, 1, , ]) #first tx, first sample, first cycle
  
  # Build an array to store the costs and QALYs accrued per cycle
  cycle_costs <- array(dim = c(n_tx, n_samples, n_cycles),
                       dimnames = list(tx_names, NULL, NULL))
  cycle_qalys <- array(dim = c(n_tx, n_samples, n_cycles),
                       dimnames = list(tx_names, NULL, NULL))
  
  cycle_tx_costs <- array(dim = c(n_tx, n_samples, n_cycles),
                          dimnames = list(tx_names, NULL, NULL))
  
  # Build arrays to store the total costs and total QALYs
  total_costs <- array(dim = c(n_tx, n_samples),
                       dimnames = list(tx_names, NULL))
  total_qalys <- array(dim = c(n_tx, n_samples),
                       dimnames = list(tx_names, NULL))
  net_benefit <- array(dim = c(n_tx, n_samples),
                       dimnames = list(tx_names, NULL))
  t5 =Sys.time()
  # Loop over the treatment options
  for (i_treatment in 1:n_tx)
  {
    # Loop over the PSA samples
    for (i_sample in 1:n_samples)
    {
      # Loop over the cycles
      # Cycle 1 is already defined so only need to update cycles 2:n_cycles
      for (i_cycle in 2:n_cycles)
      {
        # Markov update
        # Multiply previous cycle's cohort vector by transition matrix
        # i_e_ pi_j = pi_(j-1)*P
        cohort_vectors[i_treatment, i_sample, i_cycle, ] <-
          cohort_vectors[i_treatment, i_sample, i_cycle-1, ]%*%
          transition_matrices[i_treatment, i_sample, , ]
      }
      
      # total costs for each cycle
      cycle_costs[i_treatment, i_sample, ] <-
        cohort_vectors[i_treatment, i_sample, , ] %*% state_costs[i_sample, ]
      
      
      # And total QALYs for each cycle
      cycle_qalys[i_treatment, i_sample, ] <-
        cohort_vectors[i_treatment, i_sample, , ] %*% state_qalys[i_sample, ]
      
      # treatment costs for each cycle
      cycle_tx_costs[i_treatment, i_sample, ] <-
        (1.0-cohort_vectors[i_treatment, i_sample, ,n_states ]) * tx_costs[i_sample,i_treatment]
      
      # Combine the cycle_costs and treatment_costs to get total costs
      # Apply the discount factor
      # (1 in first year, 1_035 in second, 1_035^2 in third, and so on)
      total_costs[i_treatment, i_sample] <- (cycle_costs[i_treatment, i_sample, ]
                                             +cycle_tx_costs[i_treatment, i_sample, ]) %*%
        (1 / 1.035)^(0:(n_cycles - 1))
      
      # Combine the cycle_qalys to get total qalys
      # Apply the discount factor
      # (1 in first year, 1_035 in second, 1_035^2 in third, and so on)
      total_qalys[i_treatment, i_sample] <- cycle_qalys[i_treatment, i_sample, ]%*%
        (1 / 1.035)^(0:(n_cycles - 1))
      
      net_benefit[i_treatment, i_sample] <- total_qalys[i_treatment, i_sample]*threshold - total_costs[i_treatment, i_sample]
    }
  }
  t6 =Sys.time()
  
  #print(paste(t6-t5,t5-t4,t4-t3,t3-t2,t2-t1))
  
  return(list(
    total_costs = total_costs,
    total_qalys = total_qalys,
    net_benefit = net_benefit
  ))
  
}
