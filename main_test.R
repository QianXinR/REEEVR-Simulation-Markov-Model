# Set seed for random number generation
#set.seed(45843)
rm(list=ls())
# Load necessary libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl,BCEA,ggplot2,gtools,VGAM,mlmc,grid,Rcpp,doRNG,voi,simstudy, data.table, psych, digest, parallel)


#############################################################################
## Input parameters #########################################################
#############################################################################

# Define global parameters that are accessed by all functions

source("./generate_base_case_values.R")
source("./generate_net_benefit.R")
source("./generate_input_parameters.R")
source("./generate_transition_matrix.R")
source("./EVPPI_MC.R")
source("./calculate_EVPI.R")
source("./EVPPI_l_p.R")


# run_locally = FALSE
n_samples <- 1000
n_cycles <- 30
n_tx <- 3
n_states <- 3

n_columns <- n_tx * n_states * (n_states-1) + 2*n_states + n_tx -2
##run_evppi_simulation <- function(n_states){
# generate states and treatment names
state_names <- paste0("state_", 1:(n_states-1))
tx_names <- paste0("tx_", 1:n_tx)
transition_names <- generate_column_names(n_states)


# wrap evppi calculation in a function for later parallel 
run_evppi_simulation <- function(n_states, n_tx){
  
  # random states Utility
  n_states = n_states
  n_tx = n_tx
  state_utility_df <- generate_state_data(n_states, 1, "Utility", "descending")
  
  # random states costs
  state_costs_df <- generate_state_data(n_states, 10000,"Cost", "ascending" )
  
  # random drug costs per cycle
  # Initialize empty vectors to store the drug costs and standard deviations
  drug_costs <- runif(n_tx, min = 0, max = 1000)
  sd_drug_costs <- runif(n_tx, min = 0.1, max = 0.2) * drug_costs
  
  treatment_costs_df <- data.frame(
    Treatment = tx_names,
    Mean_Cost = drug_costs,
    SD_Cost = sd_drug_costs
  )
  
  treatment_costs_df <- treatment_costs_df[order(treatment_costs_df$Mean_Cost), ]
  treatment_costs_df$Treatment <- tx_names
  
  # random state transition matrix
  state_transition_matrices <- generate_tx_transition_matrices(n_tx,n_states)
  
  
  ## Inputs and results  ######################################################
  
  input_parameters <- generate_input_parameters(n_samples,treatment_costs_df, state_costs_df, state_utility_df,
                                                state_transition_matrices,hold_constant =c(), uncertain_level = 100)
  
  results <- generate_net_benefit(input_parameters, n_tx, n_samples, n_cycles,n_states, threshold = 25000)
  
  
  #############################################################################
  ## voi Analysis of results GP #################################################
  #############################################################################

  input_parameters2<-createInputs(input_parameters, print_is_linear_comb = TRUE)
  
  #utility
  evppi_utility <- voi::evppi(
    t(results$net_benefit),
    input_parameters2$mat,
    pars = paste0("utility_", state_names),
    se = TRUE,
    method="gp"
  )
  
  utility_evppi_gp = evppi_utility$evppi
  utility_evppi_gp_se = evppi_utility$se
  
  #cost
  evppi_cost <- voi::evppi(
    t(results$net_benefit),
    input_parameters2$mat,
    pars = paste0("cost_", state_names),
    se = TRUE,
    method="gp"
  )
  cost_evppi_gp = evppi_cost$evppi
  cost_evppi_gp_se = evppi_cost$se
  
  # probs
  pattern <- paste0("state_(1|2)_to_state_[1-", n_states - 1, "]")
  
  subset_list <- grep(pattern, input_parameters2$parameters, value = TRUE)
  
  evppi_prob <- voi::evppi(
    t(results$net_benefit),
    input_parameters2$mat,
    pars = subset_list,
    se = TRUE,
    method= "gp"
  )
  prob_evppi_gp = evppi_prob$evppi
  prob_evppi_gp_se = evppi_prob$se
  
  #############################################################################
  ## voi Analysis of results MARS #############################################
  #############################################################################
  
  #utility
  evppi_utility <- voi::evppi(
    t(results$net_benefit),
    input_parameters2$mat,
    pars = paste0("utility_", state_names),
    se = TRUE,
    method="earth"
  )
  
  utility_evppi_mars = evppi_utility$evppi
  utility_evppi_mars_se = evppi_utility$se
  
  #cost
  evppi_cost <- voi::evppi(
    t(results$net_benefit),
    input_parameters2$mat,
    pars = paste0("cost_", state_names),
    se = TRUE,
    method="earth"
  )
  cost_evppi_mars = evppi_cost$evppi
  cost_evppi_mars_se = evppi_cost$se

  # prob
  evppi_prob <- voi::evppi(
    t(results$net_benefit),
    input_parameters2$mat,
    pars = subset_list,
    se = TRUE,
    method= "earth"
  )
  prob_evppi_mars = evppi_prob$evppi
  prob_evppi_mars_se = evppi_prob$se
  
  #############################################################################
  ## voi Analysis of results BART #############################################
  #############################################################################
  
  #utility
  evppi_utility <- voi::evppi(
    t(results$net_benefit),
    input_parameters2$mat,
    pars = paste0("utility_", state_names),
    se = TRUE,
    method="bart"
  )
  
  utility_evppi_bart = evppi_utility$evppi
  utility_evppi_bart_se = evppi_utility$se
  
  #cost
  evppi_cost <- voi::evppi(
    t(results$net_benefit),
    input_parameters2$mat,
    pars = paste0("cost_", state_names),
    se = TRUE,
    method="bart"
  )
  cost_evppi_bart = evppi_cost$evppi
  cost_evppi_bart_se = evppi_cost$se
  
  # probs
  
  evppi_prob <- voi::evppi(
    t(results$net_benefit),
    input_parameters2$mat,
    pars = subset_list,
    se = TRUE,
    method= "bart"
  )
  prob_evppi_bart = evppi_prob$evppi
  prob_evppi_bart_se = evppi_prob$se
  
  #############################################################################
  ## Analysis of results using MC #############################################
  #############################################################################
  
  ## utility
  evppi_utility <- EVPPI_MC(100, 1000, treatment_costs_df, state_costs_df, state_utility_df, state_transition_matrices, 
                            uncertain_level = 100,  n_tx, n_cycles, n_states,hold_constant = "state_utility", threshold = 25000)
  
  utility_evppi_mc <- evppi_utility$EVPPI
  utility_evppi_mc_se <- evppi_utility$MCSE
  
  ## costs
  evppi_cost <- EVPPI_MC(100, 1000, treatment_costs_df, state_costs_df, state_utility_df, state_transition_matrices, 
                         uncertain_level = 100,  n_tx, n_cycles, n_states,hold_constant = "state_cost", threshold = 25000)
  cost_evppi_mc <- evppi_cost$EVPPI
  cost_evppi_mc_se <- evppi_cost$MCSE
  
  
  ## transition probs
  evppi_prob <- EVPPI_MC(100, 1000, treatment_costs_df, state_costs_df, state_utility_df, state_transition_matrices, 
                         uncertain_level = 100,  n_tx, n_cycles, n_states, hold_constant = "transition_probability", threshold = 25000)
  prob_evppi_mc <- evppi_prob$EVPPI
  prob_evppi_mc_se <- evppi_prob$MCSE
  
  #############################################################################
  ## Analysis of results using MLMC ###########################################
  #############################################################################
  
  # calculate evpi
  n_samples <- 100000
  
  input_parameters <- generate_input_parameters(n_samples, treatment_costs_df, state_costs_df, state_utility_df,
                                                state_transition_matrices, hold_constant = c(), uncertain_level = 100)
  
  results <- generate_net_benefit(input_parameters, n_tx, n_samples, n_cycles, n_states, threshold = 25000)
  
  evpi <- calculate_EVPI(results)
  
  #### MLMC
  
  ## utility
  EVPPI_utility_std_p <- function(M, N)
  {
    # N is the number of outer samples, M=2^l is the number of inner samples
    # Total number of samples is NN
    NN <- M * N
    
    input_parameters <-
      generate_input_parameters(
        n_samples = NN, 
        treatment_costs_df, 
        state_costs_df, 
        state_utility_df,
        state_transition_matrices, 
        hold_constant = "state_utility", 
        uncertain_level = 100
      )
    results <-
      generate_net_benefit(
        input_parameters = input_parameters,
        n_samples = NN,
        n_tx, 
        n_cycles, 
        n_states, 
        threshold = 25000
      )
    NetB <- t(results$net_benefit)
    
    return(NetB)
  }
  ## Costs
  EVPPI_cost_std_p <- function(M, N)
  {
    # N is the number of outer samples, M=2^l is the number of inner samples
    # Total number of samples is NN
    NN <- M * N
    
    input_parameters <-
      generate_input_parameters(
        n_samples = NN, 
        treatment_costs_df, 
        state_costs_df, 
        state_utility_df,
        state_transition_matrices, 
        hold_constant = "state_cost", 
        uncertain_level = 100
      )
    results <-
      generate_net_benefit(
        input_parameters = input_parameters,
        n_samples = NN,
        n_tx, 
        n_cycles, 
        n_states, 
        threshold = 25000
      )
    NetB <- t(results$net_benefit)
    
    return(NetB)
  }
  
  ## transition probabilities
  EVPPI_prob_std_p <- function(M, N)
  {
    # N is the number of outer samples, M=2^l is the number of inner samples
    # Total number of samples is NN
    NN <- M * N
    
    input_parameters <-
      generate_input_parameters(
        n_samples = NN, 
        treatment_costs_df, 
        state_costs_df, 
        state_utility_df,
        state_transition_matrices, 
        hold_constant = "transition_probability", 
        uncertain_level = 100
      )
    results <-
      generate_net_benefit(
        input_parameters = input_parameters,
        n_samples = NN,
        n_tx, 
        n_cycles, 
        n_states, 
        threshold = 25000
      )
    NetB <- t(results$net_benefit)
    
    return(NetB)
  }
  
  
  EVPPI_utility_l_p <- function(l = l,N = N) {
    return(EVPPI_l_p(l, N, EVPPI_std_p = EVPPI_utility_std_p))
  }
  
  
  EVPPI_cost_l_p <- function(l = l,N = N) {
    return(EVPPI_l_p(l, N, EVPPI_std_p = EVPPI_cost_std_p))
  }
  
  
  EVPPI_prob_l_p <- function(l = l,N = N) {
    return(EVPPI_l_p(l, N, EVPPI_std_p = EVPPI_prob_std_p))
  }
  
  #set.seed(33)
  tst_utility <- mlmc.test(EVPPI_utility_l_p, N=1000,
                           L=4, N0=1000,
                           eps.v=c(60, 30, 15,7),
                           Lmin=2, Lmax=10,
                           alpha = 1,
                           beta = 1.5)
  utility_evppi_mlmc <- evpi$EVPI - tst_utility$P[5]
  
  
  tst_cost <- mlmc.test(EVPPI_cost_l_p, N=1000,
                        L=4, N0=1000,
                        eps.v=c(60, 30, 15,7),
                        Lmin=2, Lmax=10,
                        alpha = 1,
                        beta = 1.5)
  cost_evppi_mlmc <- evpi$EVPI - tst_cost$P[5]
  
  
  tst_prob <- mlmc.test(EVPPI_prob_l_p, N=1000,
                        L=4, N0=1000,
                        eps.v=c(60, 30, 15,7),
                        Lmin=2, Lmax=10,
                        alpha = 1,
                        beta = 1.5)
  prob_evppi_mlmc <- evpi$EVPI - tst_prob$P[5]
  #############################################################################
  ## save the results #########################################################
  #############################################################################
  results_vector <- c(utility_evppi_gp, utility_evppi_gp_se, 
                      cost_evppi_gp, cost_evppi_gp_se,
                      prob_evppi_gp, prob_evppi_gp_se,
                      utility_evppi_bart, utility_evppi_bart_se, 
                      cost_evppi_bart, cost_evppi_bart_se,
                      prob_evppi_bart, prob_evppi_bart_se,
                      utility_evppi_mars, utility_evppi_mars_se, 
                      cost_evppi_mars, cost_evppi_mars_se,
                      prob_evppi_mars, prob_evppi_mars_se,
                      utility_evppi_mc,utility_evppi_mc_se,
                      cost_evppi_mc, cost_evppi_mc_se,
                      prob_evppi_mc, prob_evppi_mc_se,
                      utility_evppi_mlmc,
                      cost_evppi_mlmc,
                      prob_evppi_mlmc
  )
  
  results_names <- c('utility_evppi_gp', 'utility_evppi_gp_se', 
                     'cost_evppi_gp', 'cost_evppi_gp_se',
                     'prob_evppi_gp', 'prob_evppi_gp_se',
                     'utility_evppi_bart', 'utility_evppi_bart_se', 
                     'cost_evppi_bart', 'cost_evppi_bart_se',
                     'prob_evppi_bart', 'prob_evppi_bart_se',
                     'utility_evppi_mars', 'utility_evppi_mars_se', 
                     'cost_evppi_mars', 'cost_evppi_mars_se',
                     'prob_evppi_mars', 'prob_evppi_mars_se',
                     'utility_evppi_mc', 'utility_evppi_mc_se',
                     'cost_evppi_mc', 'cost_evppi_mc_se',
                     'prob_evppi_mc', 'prob_evppi_mc_se',
                     'utility_evppi_mlmc',
                     'cost_evppi_mlmc',
                     'prob_evppi_mlmc'
  )
  names(results_vector) <- results_names
  df <- as.data.frame(t(results_vector))
  return(df)
  
}

# a <- run_evppi_simulation(3)


if (run_locally == FALSE) {
  library(digest)
  script_run_hash = digest(format(Sys.time(), "%Y%m%d%H%M%S") , algo = 'md5')
  library(parallel)
  num_sims = 5 # can be changed
  save_window = 1
  input_values <- rep(c(n_states,n_tx), save_window)
  num_cores <- detectCores() / 2 + 10 #number of CPU cores available on the machine
  num_cores_to_use <- min(42, num_cores, detectCores()) # capped the number of cores at 42 to avoid overloading the system
  
  for (n_states in 3:10) {
    # Loop through n_states from 3 to 10
    for (n_tx in 2:10) {
      # Inner loop for n_tx from 2 to 10
      
      input_values <- rep(n_states, save_window)
      
      for (iteration in 1:ceiling(num_sims / save_window)) {
        print(
          paste(
            "Running for n_states =",
            n_states,
            "and n_tx =",
            n_tx,
            "using",
            num_cores_to_use,
            "cores"
          )
        )
        
        c1 <- makeCluster(num_cores_to_use) # creates a cluster cores, allowing parallel execution on multiple CPU cores
        
        # import my own libraries
        clusterEvalQ(c1, source("./generate_base_case_values.R"))
        clusterEvalQ(c1, source("./generate_net_benefit.R"))
        clusterEvalQ(c1, source("./generate_input_parameters.R"))
        clusterEvalQ(c1, source("./generate_transition_matrix.R"))
        clusterEvalQ(c1, source("./EVPPI_MC.R"))
        clusterEvalQ(c1, source("./calculate_EVPI.R"))
        clusterEvalQ(c1, source("./EVPPI_l_p.R"))
        
        # import external libraries
        
        #clusterEvalQ(c1, {library(MCMCpack)})
        clusterExport(
          c1,
          c(
            "n_samples",
            "n_cycles",
            "n_columns",
            "state_names",
            "tx_names",
            "transition_names"
          )
        )
        clusterExport(c1, c('n_states'))
        clusterExport(c1,
                      c(
                        "rdirichlet",
                        'rdiric',
                        'bcea',
                        'createInputs',
                        'evppi'
                      )) # n_states is defined as global var, not right, to be fixed
        clusterExport(c1, c("run_evppi_simulation", 'input_values'))
        results <- parLapply(c1, input_values, run_evppi_simulation)
        
        # bind the results together
        if (iteration == 1) {
          df_results <- do.call(rbind, results)
        } else {
          df_results_temp <- do.call(rbind, results)
          df_results <- rbind(df_results, df_results_temp)
        }
        
        print(df_results)
        
        # save results to a CSV file
        time_code <- format(Sys.time(), "%Y%m%d%H%M%S")
        file_path <- paste0(
          './results_folder/results_',
          script_run_hash,
          '_',
          n_states,
          '_',
          time_code,
          '.csv'
        )
        write.csv(df_results,
                  file = gsub("\\s+", "", file_path),
                  row.names = FALSE)
        
        message <- paste(
          iteration,
          '.',
          save_window * iteration,
          '/',
          num_sims,
          'jobs completed, job hash:',
          script_run_hash
        )
        print(message)
        
        # Stop the cluster
        stopCluster(c1)
        
      }
    }
  }
} else {
  for (n_states in 3:10) {
    for (n_tx in 2:10) {
      run_evppi_simulation(n_states, n_tx)  # run locally if run_locally is TRUE
    }
  }
}








