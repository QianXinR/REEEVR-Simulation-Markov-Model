# REEEVR - Assessing reliability limits of regression and approximation methods for VoI estimation 
# Develop a generic Markov model that can have any number of states and input parameters with varying level of correlation. 
# Qian Xin October 2023

# Set seed for random number generation
#set.seed(45843)

# Load necessary libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl,BCEA,ggplot2,gtools,VGAM,mlmc,grid,Rcpp,doRNG,voi,simstudy, data.table, psych)

# Challenges and assumptions
#' 1- The random numbers generate for drug costs are fairly close to each other, which might not be the case in real life
#'  1.1- Generate state costs values in increasing order, i.e., the more severer the health sate, the higher costs correlated.
#'        The final health state(death) is assumed to have no associated costs.
#'  1.2- Generate state utility values in decreasing order, i.e., the more severer the health sate, the lower utility correlated.
#'        The final health state(death) is assumed to have no associated utility.
#' 2- Assumed patient cohort all start from health state 1, and the the cohort changes with each model cycle
#' 3- Assumed the cycle length is one year, so no further calculation is conducted transferring utility to QALYs
#' 4- Assumed no treatment discontinuation is considered, patients taking the treatments within the time horizon,
#     except when they moved into the last health states    


#############################################################################
## Input parameters #########################################################
#############################################################################
source("./generate_base_case_values.R")
source("./generate_net_benefit.R")
source("./generate_input_parameters.R")
source("./generate_transition_matrix.R")
source("./EVPPI_MC.R")
source("./calculate_EVPI.R")
source("./EVPPI_l_p.R")

# Define global parameters that are accessed by all functions

n_samples <- 1000
n_states <- 5
n_tx <- 3
n_cycles <- 30
n_columns <- n_tx * n_states * (n_states-1) + 2*n_states + n_tx -2

# generate states and treatment names
state_names <- paste0("state_", 1:(n_states-1))
tx_names <- paste0("tx_", 1:n_tx)
transition_names <- generate_column_names(n_states)

# random states Utility
state_utility_df <- generate_state_data(n_states, 1, "Utility", "descending")
#print(state_utility_df)

# random states costs
state_costs_df <- generate_state_data(n_states, 10000,"Cost", "ascending" )
#print(state_costs_df)


# random drug costs per cycle
# Initialize empty vectors to store the drug costs and standard deviations
drug_costs <- runif(n_tx, min = 0, max = 10000)
sd_drug_costs <- runif(n_tx, min = 0.1, max = 0.2) * drug_costs

treatment_costs_df <- data.frame(
  Treatment = tx_names,
  Mean_Cost = drug_costs,
  SD_Cost = sd_drug_costs
)

treatment_costs_df <- treatment_costs_df[order(treatment_costs_df$Mean_Cost), ]
treatment_costs_df$Treatment <- tx_names

# Print the final data frame
#print(treatment_costs_df)

# random state transition matrix
state_transition_matrices <- generate_tx_transition_matrices(n_tx,n_states)

#state_transition_matrices


#############################################################################
## Cohort Simulation ########################################################
#############################################################################

# Till this point, we have
#' inputs for n_samples
#' transition matrix for n_tx,n_samples

# In this script, we want to
#' simulate a cohort all start from health state 1, and the the cohort changes with each model cycle
#' calculate the total costs and total QALYs for each cycle
#'
#'

#n_samples <- 500000

input_parameters <- generate_input_parameters(n_samples, treatment_costs_df, state_costs_df, state_utility_df,
                                              state_transition_matrices,hold_constant = c(), uncertain_level = 100)
results <- generate_net_benefit(input_parameters, n_tx, n_samples, n_cycles, n_states, threshold = 25000)

evpi_estimate <- calculate_EVPI(results)
evpi_estimate

#############################################################################
## voi Analysis of results GP #################################################
#############################################################################
# method_used <- if (n_states > 4) "gp" else "gam"
input_parameters2<-createInputs(input_parameters, print_is_linear_comb = TRUE)

#utility
evppi_utility <- voi::evppi(
  t(results$net_benefit),
  input_parameters2$mat,
  pars = paste0("utility_", state_names),
  se = TRUE,
  method="gp"
)
evppi_utility

#cost
evppi_cost <- voi::evppi(
  t(results$net_benefit),
  input_parameters2$mat,
  pars = paste0("cost_", state_names),
  se = TRUE,
  method="gp"
)
evppi_cost

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
evppi_prob

#############################################################################
## voi Analysis of results MARS #############################################
#############################################################################
# # method_used <- if (n_states > 4) "gp" else "gam"
# input_parameters2<-createInputs(input_parameters, print_is_linear_comb = TRUE)

#utility
evppi_utility <- voi::evppi(
  t(results$net_benefit),
  input_parameters2$mat,
  pars = paste0("utility_", state_names),
  se = TRUE,
  method="earth"
)
evppi_utility

#cost
evppi_cost <- voi::evppi(
  t(results$net_benefit),
  input_parameters2$mat,
  pars = paste0("cost_", state_names),
  se = TRUE,
  method="earth"
)
evppi_cost

# probs
# pattern <- paste0("state_(1|2)_to_state_[1-", n_states - 1, "]")
# 
# subset_list <- grep(pattern, input_parameters2$parameters, value = TRUE)

evppi_prob <- voi::evppi(
  t(results$net_benefit),
  input_parameters2$mat,
  pars = subset_list,
  se = TRUE,
  method= "earth"
)
evppi_prob

#############################################################################
## voi Analysis of results BART #############################################
#############################################################################
# method_used <- if (n_states > 4) "gp" else "gam"
# input_parameters2<-createInputs(input_parameters, print_is_linear_comb = TRUE)

#utility
evppi_utility <- voi::evppi(
  t(results$net_benefit),
  input_parameters2$mat,
  pars = paste0("utility_", state_names),
  se = TRUE,
  method="bart"
)
evppi_utility

#cost
evppi_cost <- voi::evppi(
  t(results$net_benefit),
  input_parameters2$mat,
  pars = paste0("cost_", state_names),
  se = TRUE,
  method="bart"
)
evppi_cost

# probs
# pattern <- paste0("state_(1|2)_to_state_[1-", n_states - 1, "]")
# 
# subset_list <- grep(pattern, input_parameters2$parameters, value = TRUE)

evppi_prob <- voi::evppi(
  t(results$net_benefit),
  input_parameters2$mat,
  pars = subset_list,
  se = TRUE,
  method= "bart"
)
evppi_prob

#############################################################################
## Analysis of results using MC #############################################
#############################################################################

## utility
evppi_utility <- EVPPI_MC(100, 1000, treatment_costs_df, state_costs_df, state_utility_df, state_transition_matrices, 
                          uncertain_level = 100,  n_tx, n_cycles, n_states,hold_constant = "state_utility", threshold = 25000)

evppi_utility

## costs
evppi_cost <- EVPPI_MC(100, 1000, treatment_costs_df, state_costs_df, state_utility_df, state_transition_matrices, 
                       uncertain_level = 100,  n_tx, n_cycles, n_states,hold_constant = "state_cost", threshold = 25000)
evppi_cost


## transition probs
evppi_prob <- EVPPI_MC(100, 1000, treatment_costs_df, state_costs_df, state_utility_df, state_transition_matrices, 
                       uncertain_level = 100,  n_tx, n_cycles, n_states, hold_constant = "transition_probability", threshold = 25000)
evppi_prob


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
                         eps.v=c(60, 30, 15,7,3,1),
                         Lmin=2, Lmax=10)


tst_cost <- mlmc.test(EVPPI_cost_l_p, N=1000,
                      L=4, N0=1000,
                      eps.v=c(60, 30, 15,7,3,1),
                      Lmin=2, Lmax=10)

tst_prob <- mlmc.test(EVPPI_prob_l_p, N=1000,
                      L=4, N0=1000,
                      eps.v=c(60, 30, 15,7,3,1),
                      Lmin=2, Lmax=10)
