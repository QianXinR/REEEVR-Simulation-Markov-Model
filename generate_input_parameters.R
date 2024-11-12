# REEEVR - Assessing reliability limits of regression and approximation methods for VoI estimation 
# Develop a generic Markov model that can have any number of states and input parameters with varying level of correlation. 
# Qian Xin October 2023

# Generate input parameters data with n_samples
generate_input_parameters <- function(n_samples, n_states, n_tx, treatment_costs_df, state_costs_df, state_utility_df,
                                      state_transition_matrices,transition_names,state_names, tx_names,
                                      n_columns,hold_constant = c(), uncertain_level = 100)  {
  
  input_parameters <- as.data.frame(matrix(nrow = n_samples, ncol = n_columns))
  colnames(input_parameters) <- c(transition_names,
                                  paste0("utility_", state_names),
                                  paste0("cost_", state_names),
                                  paste0("cost_", tx_names))
  
  # treatment costs, state costs, and state utility are normally distributed
  # For treatment costs
  input_parameters <- populate_random_value_normal(input_parameters, n_samples, treatment_costs_df, "cost_", "Treatment", "Mean_Cost", "SD_Cost")
  
  # For state costs
  if(!is.element("state_cost", hold_constant)){
    input_parameters <- populate_random_value_normal(input_parameters, n_samples, state_costs_df, "cost_", "Health_State", "Mean_Cost", "SD_Cost")
  }else{
    input_parameters <- populate_random_value_normal(input_parameters, 1, state_costs_df, "cost_", "Health_State", "Mean_Cost", "SD_Cost")
  }
  
  
  # For state utility
  if (!is.element("state_utility", hold_constant)) {
    input_parameters <- populate_random_value_normal(input_parameters, n_samples, state_utility_df, "utility_", "Health_State", "Mean_Utility", "SD_Utility")
    # cap the values between 0 and 1 for all utility columns
    utility_columns <- grep("^utility_", colnames(input_parameters), value = TRUE)
    for (col in utility_columns) {
      input_parameters[[col]] <- pmin(pmax(input_parameters[[col]], 0), 1)
    }
    
  } else {
    input_parameters <- populate_random_value_normal(input_parameters, 1, state_utility_df, "utility_", "Health_State", "Mean_Utility", "SD_Utility")
    # cap the values between 0 and 1 for all utility columns
    utility_columns <- grep("^utility_", colnames(input_parameters), value = TRUE)
    for (col in utility_columns) {
      input_parameters[[col]] <- pmin(pmax(input_parameters[[col]], 0), 1)
  }
  }
  
  # Transition probability
  # Sample from the Dirichlet distribution for each row, except the absorbing state
  if(!is.element("transition_probability", hold_constant)){
    samples <- matrix(0, nrow = n_samples, ncol = n_tx * n_states * (n_states - 1))
    for (n in 1:n_tx){
      for (i in 1:(n_states - 1)) {
        matrix_name <- paste0("tx_", n, "_state_transition_matrix")
        # Calculate the start and end columns for the current iteration
        start_col = ((n - 1) * (n_states - 1) + (i - 1)) * n_states + 1
        end_col = ((n - 1) * (n_states - 1) + i) * n_states
        
        # Fill the appropriate columns of the matrix
        samples[, start_col:end_col] <- rdiric(n_samples, state_transition_matrices[[matrix_name]][i, ]*uncertain_level)
        input_parameters[, 1:(n_tx * n_states * (n_states - 1))] <- samples
      }
    }}else{
      samples <- matrix(0, nrow = 1, ncol = n_tx * n_states * (n_states - 1))
      for (n in 1:n_tx){
        for (i in 1:(n_states - 1)) {
          matrix_name <- paste0("tx_", n, "_state_transition_matrix")
          # Calculate the start and end columns for the current iteration
          start_col = ((n - 1) * (n_states - 1) + (i - 1)) * n_states + 1
          end_col = ((n - 1) * (n_states - 1) + i) * n_states
          
          # Fill the appropriate columns of the matrix
          samples[, start_col:end_col] <- rdiric(1, state_transition_matrices[[matrix_name]][i, ]*uncertain_level)
          replicated_samples <- matrix(rep(t(samples), n_samples), nrow = n_samples, byrow = TRUE)
          input_parameters[, 1:(n_tx * n_states * (n_states - 1))] <- replicated_samples
        
        }
      }
    }
  
  
  return(input_parameters)
  
}

# samples[, 4:6] <- rdiric(1, state_transition_matrices$tx_1_state_transition_matrix[2, ]*100)
# replicated_samples <- matrix(rep(t(samples), n_samples), nrow = n_samples, byrow = TRUE)
