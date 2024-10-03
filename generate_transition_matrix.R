# REEEVR - Assessing reliability limits of regression and approximation methods for VoI estimation 
# Develop a generic Markov model that can have any number of states and input parameters with varying level of correlation. 
# Qian Xin October 2023

# generate the transition matrices for each sample, each tx
generate_transition_matrices <- function(input_parameters,n_samples,n_tx, n_states ) { # move this outside the function
  # Create a 3-dimensional array to store the transition matrices
  # Dimensions: treatments x samples x states x states
  transition_matrices <- array(0, dim = c(n_tx, n_samples, n_states, n_states),
                               dimnames = list(tx_names, NULL, NULL, NULL)) # Tx, sample, states, states
  
  for (tx in 1:n_tx) {
    for (sample in 1:n_samples) {
      for (state in 1:(n_states - 1)) {
        # Construct the column name pattern for the current treatment and state
        col_pattern <- paste0("tx_", tx, "_state_", state, "_to_state_")
        
        # Extract the columns corresponding to this pattern
        relevant_cols <- grep(col_pattern, colnames(input_parameters))
        
        # Assign the probabilities to the transition matrix
        transition_matrices[tx, sample, state, 1:n_states] <- as.numeric(input_parameters[sample, relevant_cols])
      }
      
      # For the absorbing state, it remains 1
      transition_matrices[tx, sample, n_states, n_states] <- 1
    }
  }
  
  return(transition_matrices)
}
