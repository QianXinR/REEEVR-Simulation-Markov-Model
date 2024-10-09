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
#results$net_benefit
#############################################################################
## Analysis of results ######################################################
#############################################################################
#results$total_costs
#results$total_qalys
#results$net_benefit
# Average costs
#average_costs <- rowMeans(results$total_costs)

# Average effects
#average_effects <- rowMeans(results$total_qalys)

# Average net benefit
#average_net_benefit <- rowMeans(results$net_benefit)

# use bcea package to calculate the model results
model_bcea <- bcea(e = t(results$total_qalys),
                   c = t(results$total_costs), ref = 1,
                   interventions = tx_names,
                   k=c(20000,25000) )

summary(model_bcea)


## Calculate EVPPI use the BCEA package
input_parameters2<-createInputs(input_parameters, print_is_linear_comb = TRUE)


# Determine the method based on the number of states
# method_used <- if (n_states > 4) "gp" else "gam"

# Calculate EVPPI
## state utilities
evppi_utilities <- evppi(
  he = model_bcea,
  param_idx = paste0("utility_", state_names),
  input = input_parameters2$mat,
  method = "gp"
)
evppi_utilities$evppi

## state costs
evppi_state_costs <- evppi(
  he = model_bcea,
  param_idx = paste0("cost_", state_names),
  input = input_parameters2$mat,
  method = "gp"
)
evppi_state_costs$evppi

## transition probabilities
#param_names <- c()  # Initialize an empty vector to store parameter names

# generate parameter names
# for (tx in 1:n_tx) {
#   for (state in 1:2) {
#     for (state2 in 1:(n_states - 1)) {
#       param_name <-
#         paste0("tx_", tx, "_state_", state, "_to_state_", state2)
#       param_names <-
#         c(param_names, param_name)  # Accumulate parameter names
#     }
#   }
# }


# Remove the transition probability to the last health state
pattern <- paste0("state_(1|2)_to_state_[1-", n_states - 1, "]")
subset_list <-
  grep(pattern, input_parameters2$parameters, value = TRUE)

# Calculate EVPPI for the combined set of parameters
evppi_transition_probs <- evppi(
  he = model_bcea,
  param_idx = input_parameters2$parameters,
  input = input_parameters2$mat,
  method = 'gp'
)

#evppi_result$evppi
evppi_transition_probs$evppi

#########MLMC###########


## Utility
EVPPI_utility_std_p <- function(M, N)
{
  # N is the number of outer samples, M=2^l is the number of inner samples
  # Total number of samples is NN
  NN <- M * N
  
  input_parameters <-
    generate_input_parameters(n_samples = NN, treatment_costs_df, state_costs_df, state_utility_df,
                              state_transition_matrices,hold_constant = "state_utility",uncertain_level = 50)
  results <-
    generate_net_benefit(
      input_parameters = input_parameters,
      n_tx = n_tx,
      n_samples = NN,
      n_states,
      n_cycles = n_cycles
    )
  NetB <- t(results$net_benefit)
  
  return(NetB)
}

## state costs
EVPPI_cost_std_p <- function(M, N)
{
  # N is the number of outer samples, M=2^l is the number of inner samples
  # Total number of samples is NN
  NN <- M * N
  
  input_parameters <-
    generate_input_parameters(n_samples = NN, treatment_costs_df, state_costs_df, state_utility_df,
                              state_transition_matrices, hold_constant = "state_cost",uncertain_level = 50)
  results <-
    generate_net_benefit(
      input_parameters = input_parameters,
      n_tx = n_tx,
      n_samples = NN,
      n_states,
      n_cycles = n_cycles
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
      n_samples = NN,treatment_costs_df, state_costs_df, state_utility_df,
      state_transition_matrices,
      hold_constant = "transition_probability",
      uncertain_level = 50
    )
  results <-
    generate_net_benefit(
      input_parameters = input_parameters,
      n_tx = n_tx,
      n_samples = NN,
      n_states,
      n_cycles = n_cycles
    )
  NetB <- t(results$net_benefit)
  
  return(NetB)
}

# Wrapper function to generate the net benefit holding the parameters of interest constant
# Repeats parameters of interest M times and generates random for remainder
# Calculates net benefit based on these

#a<-EVPPI_utility_std_p(3,2)
#a
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
  
  return(list(sums = sum1, cost = M * N))
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

#EVPPI_utility_l_p(2, 5)


#set.seed(33)
tst_utility <- mlmc.test(EVPPI_utility_l_p, N=10000,
                         L=4, N0=1000,
                         eps.v=c(60, 30, 15),
                         Lmin=2, Lmax=10)

mlmc_utility <- mlmc(
  Lmin = 2,
  Lmax = 10,
  N0 = 1000,
  eps =1,
  EVPPI_utility_l_p,
  alpha = NA,
  beta = NA,
  gamma = NA,
  parallel = NA
)

tst_cost <- mlmc.test(EVPPI_cost_l_p, M=2, N=1024,
                      L=2, N0=128,
                      eps.v=c(60, 30, 15),
                      Lmin=2, Lmax=10)

tst_prob <- mlmc.test(EVPPI_prob_l_p, M=2, N=1024,
                      L=2, N0=128,
                      eps.v=c(60, 30, 15),
                      Lmin=2, Lmax=10)
