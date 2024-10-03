
calculate_EVPI <- function(results) {
  
  # Calculate the max across treatments for each sample
  Max_NetB <- apply(results$net_benefit, 2, max)
  
  
  # Calculate the mean for each treatment across all samples
  Mean_NetB <- apply(results$net_benefit, 1, mean)
  
  # Find the max of the mean across treatments
  Max_Mean_NetB <- max(Mean_NetB)
  
  # difference for each cycle
  Diff <- Max_NetB - Max_Mean_NetB
  
  # EVPI 
  EVPI <- mean(Diff)
  
  # MCSE
  MCSE <- sd(Diff) / sqrt(length(Diff))
  
  return(data.frame(EVPI, MCSE))
}

