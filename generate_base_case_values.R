# REEEVR - Assessing reliability limits of regression and approximation methods for VoI estimation 
# Develop a generic Markov model that can have any number of states and input parameters with varying level of correlation. 
# Qian Xin October 2023

# Generate column names
# Omitting the final state column, assumed to be an absorbing state with constant transition probabilities of 1.
generate_column_names <- function(n_states) {
  column_names <- character(0)  # Initialize an empty character vector
  for (i in 1:n_tx){
    for (n in 1:(n_states-1)) {
      for (m in 1:n_states) {
        transition_name <- paste0("tx_",i, "_state_", n, "_to_state_", m)
        column_names <- c(column_names, transition_name)
      }
    }
  }
  return(column_names)
}

# Generate random number
generate_state_data <- function(n_states, n_max, metric_type, sort_direction) {
  # Generate random state metrics and percentages, then calculate standard deviations
  state_metric <- runif(n_states - 1, min = 0, max = n_max)
  sd_state_metric <- runif(n_states - 1, min = 0.4, max = 0.5) * state_metric
  
  # Initialize the data frame
  state_data_df <- data.frame(
    Health_State = paste0("state_", 1:(n_states - 1)),
    Mean = state_metric,
    SD = sd_state_metric
  )
  
  # Sort by Mean and rearrange SD accordingly
  if (sort_direction == "ascending") {
    sorted_indices <- order(state_data_df$Mean)
  } else if (sort_direction == "descending") {
    sorted_indices <- order(-state_data_df$Mean)
  }
  
  state_data_df[c("Mean", "SD")] <- state_data_df[sorted_indices, c("Mean", "SD")]
  
  # Rename the columns based on metric_type
  new_colnames <- c(paste0("Mean_", metric_type), paste0("SD_", metric_type))
  colnames(state_data_df)[2:3] <- new_colnames
  
  return(state_data_df)
}

# Function to initialize mu for different health states
initialize_mu <- function(n_states) {
  # Always start with 1 and 0.4, 0.2, followed by (n-3) elements of 0.1
  if (n_states >= 3) {
    return(c(1, 0.4, 0.2, rep(0.1, n_states - 3)))
  } else {
    stop("n must be 3 or greater.")
  }
}

# Function to generate correlation matrix with 0.25 for adjacent states and 0.1 otherwise
generate_cor_matrix <- function(n_states) {
  cor_matrix <- diag(1, n_states)  # Start with identity matrix (diagonal elements = 1)
  
  # Fill the adjacent elements with 0.25
  for (i in 1:(n_states - 1)) {
    cor_matrix[i, i + 1] <- 0.25
    cor_matrix[i + 1, i] <- 0.25
  }
  
  # Fill the rest of the off-diagonal elements with 0.1
  cor_matrix[cor_matrix == 0] <- 0.1
  
  return(cor_matrix)
}

# Internal function to generate a random transition matrix for a given number of states
generate_matrix <- function(n_states, uncertain_level) {
  
  # Initialize mu based on n_states
  mu <- initialize_mu(n_states*n_states)
  # Initialize a correlation matrix based on n_states
  corMatrix <- generate_cor_matrix(n_states*n_states)
  # corMatrix_normalised <- t(apply(corMatrix, 1, function(row) row / sum(row)))
  
  # Initialize an empty n x n matrix filled with zeros
  mat <- matrix(0, n_states, n_states)
  
  # Generate n_states*n_states observations of n_states RVs from a multivariate normal 
  dt <- genCorData(1, mu = mu, sigma = 1, 
                   corMatrix = corMatrix, corstr = "cs" )
  
  # dt <- as.matrix(dt[, -1]) %*% corMatrix_normalised
  # dt <- as.data.table(cbind(id = 1:nrow(dt), dt))
  
  ### create a long version of the data set
  
  dtM <- melt(dt, id.vars = "id", variable.factor = TRUE, 
              value.name = "X", variable.name = "seq")
  setkey(dtM, "id")   # sort data by id
  dtM[, seqid := .I]  # add index for each record
  
  ### apply CDF to X to get uniform distribution
  
  dtM[, U := pnorm(X)]
  
  ### Generate correlated Poisson data with mean and variance 
  ### apply inverse CDF to U
  
  dtM[, Y_beta := qbeta(U, 4, 19), keyby = seqid]
  
  mat <- as.matrix(dcast(data = dtM, id~seq,value.var = "Y_beta"))
  mat <- as.matrix(mat[, -1])
  mat <- matrix(mat, nrow = n_states, ncol = n_states, byrow = TRUE)
  mat <- mat / rowSums(mat)

  # Update the last row to be all zeros, except the last element which will be uncertain_level
  mat[n_states, ] <- c(rep(0, n_states - 1), 1)
  mat <- mat * uncertain_level
  
  # Assign row and column names
  row_col_names <- paste0("state_", 1:n_states)
  rownames(mat) <- colnames(mat) <- row_col_names
  
  return(mat)
}

# Generate treatment specific transition matrices
generate_tx_transition_matrices <- function(n_tx, n_states, uncertain_level = 1) {
  
  # Initialize a list to store the matrices for each treatment
  tx_transition_matrices <- list()
  
  tx_transition_matrices[["tx_1_state_transition_matrix"]] <- generate_matrix(n_states, uncertain_level)
  
  alphas_dirichlet <- list()
  delta_values <- list()
  scales <- list()
  
  # Loop through each treatment
  for (n in 2:n_tx) {
    
    # tx_transition_matrices[[paste0("tx_", 1, "_state_transition_matrix")]] <- generate_matrix(n_states, uncertain_level)
    
    #get alphas to use for dirichlet distribution for the delta
    alphas_dirichlet[[n]] <- runif(n_states)^2
    
    delta_values[[n]]  <- rdirichlet(n_states, alphas_dirichlet[[n]])
    
    # get the minimum value for each row in initial transition matrix
    scales <- apply( tx_transition_matrices[["tx_1_state_transition_matrix"]],1,min)
    #print(scales)
    # center to sum =0 and scale by the minimum value of that row
    delta_values[[n]] <- diag(scales)%*% (delta_values[[n]]- 1/n_states) #sweep(delta_values[[n]]- 1/n_states, 1, scales, `*`)
    
    delta_values[[n]]  <- t(apply(delta_values[[n]],1,sort))
    
    tx_transition_matrices[[paste0("tx_", n, "_state_transition_matrix")]] <-
      tx_transition_matrices[["tx_1_state_transition_matrix"]] - delta_values[[n]]
    
    # paste the row 3 to n_states transition probabilities from tx1 to the rest of the txs
    tx_transition_matrices[[paste0("tx_", n, "_state_transition_matrix")]][3:n_states,] <-
      tx_transition_matrices[["tx_1_state_transition_matrix"]][3:n_states, ]
    
    # swap first column values with diagonal values for all matrices in the tx_transition_matrices list, except for the last row (death)
    # Loop through each matrix in tx_transition_matrices
    for (m in names(tx_transition_matrices)) {
      # Extract the current matrix
      current_matrix <- tx_transition_matrices[[m]]
      
      # Store the first column and diagonal values in temporary variables
      first_column <- current_matrix[1:(n_states - 1), 1]
      diagonal_values <- diag(current_matrix)[1:(n_states - 1)]
      
      # Swap the values
      diag(current_matrix)[1:(n_states - 1)] <- first_column
      current_matrix[1:(n_states - 1), 1] <- diagonal_values
      
      # Update the matrix in the list
      tx_transition_matrices[[m]] <- current_matrix
    }
    
    
  }
  
  # Return the list of generated matrices
  return(tx_transition_matrices)
}


# Generate function to calculate random values based on normal distribution
populate_random_value_normal <- function(input_parameters, n_samples, df, column_prefix, id_col, mean_col, sd_col) {
  for (i in 1:nrow(df)) {
    item_name <- df[i, id_col]
    mean_value <- df[i, mean_col]
    sd_value <- df[i, sd_col]
    
    # Generate the column name
    column_name <- paste0(column_prefix, item_name)
    
    # Populate the column in input_parameters
    input_parameters[[column_name]] <- rnorm(n_samples, mean = mean_value, sd = sd_value)
  }
  return(input_parameters)
}




#transition_matrices<- generate_transition_matrices(input_parameters)


