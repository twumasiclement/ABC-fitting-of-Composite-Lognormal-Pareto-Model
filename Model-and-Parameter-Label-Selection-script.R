
# Function to simulate claim severity data from a truncated Lognormal-Pareto distribution
simulate_truncated_lognormal_pareto <- function(N, mu, sigma, alpha, y_m, truncation_threshold = 500000) {
  # Here, y_m=theta=truncation point
  claims <- numeric(N)  # Initialize vector to store claims
  claim_type<- numeric(N)
  
  # Compute survival probabilities of choosing a particular type of claim between lognormal and Pareto
  P_L <- 1 - plnorm(truncation_threshold, meanlog = mu, sdlog = sigma)  # Lognormal survival
  P_P <- (y_m / truncation_threshold)^alpha  # Pareto survival function
  
  # Adjust probability for Pareto selection
  p_corrected <- (0.3 * P_L) / (0.7 * P_P + 0.3 * P_L)
  
  # Ensure p_corrected is reasonable (avoids all claims being lognormal)
  p_corrected <- max(min(p_corrected, 0.5), 0.1)  # Keep within reasonable bounds
  
  
  for (i in 1:N) {
    repeat {
      U <- runif(1)  # Uniform random variable for deciding between lognormal & Pareto
      
      # Generate a Lognormal claim
      lognormal_claim <- rlnorm(1, meanlog = mu, sdlog = sigma)
      
      
      # Generate a Pareto claim for large severity cases
      pareto_claim <-  EnvStats::rpareto(1, location=y_m, shape=alpha)
      
      # Choose between Lognormal and Pareto based on U (Over 70% Lognormal, at most 30% Pareto)
      # If U> 1 - p_corrected, at most 30% of the time, Pareto claim is selected.
      claim <- ifelse(U > (1 - p_corrected), pareto_claim, lognormal_claim)
      claim_category<-  ifelse(U > (1 - p_corrected), "Large Claim", "Small/Moderate Claim")
      
      # Checking the claim is above the truncation threshold
      if (claim >= truncation_threshold  & lognormal_claim<y_m) {
        claims[i] <- claim/10^6 # severity data in millions
        claim_type[i]<- claim_category
        break
      }
    }
  }
  # No Copula was used here
  return(data.frame(Claim_Severity = claims, Copula = "Lognormal-Pareto", Claim_type=claim_type))
}







# A function to select the Claim Severity simulator of interest
simulation_model_selection<- function(N,copula_type, parameters){
  
  # Selecting the simulator of interest
  if (copula_type == "Lognormal-Pareto"){
    # Independent Lognormal-Pareto
    mu=exp(parameters[1]);sigma=exp(parameters[2]);alpha=exp(parameters[3]); y_m=exp(parameters[4])
    sim.output<- simulate_truncated_lognormal_pareto(N, mu, sigma, alpha, y_m)
  }
  if (copula_type == "Gaussian") {
    # Gaussian Copula-Based Lognormal-Pareto Simulation
    mu=exp(parameters[1]);sigma=exp(parameters[2]);alpha=exp(parameters[3]); 
    y_m=exp(parameters[4]);rho=exp(parameters[5])
    sim.output <- simulate_truncated_copula_LognormalPareto(N, mu, sigma, alpha, y_m,
                                                            "Gaussian", rho,truncation_threshold = 500000)
  } else if (copula_type == "Clayton") {
    # Clayton Copula-Based Lognormal-Pareto Simulation
    mu=exp(parameters[1]);sigma=exp(parameters[2]);alpha=exp(parameters[3]); 
    y_m=exp(parameters[4]);theta=exp(parameters[5])
    sim.output <- simulate_truncated_copula_LognormalPareto(N, mu, sigma, alpha,
                                                            y_m, "Clayton", theta,truncation_threshold = 500000)
  } else if (copula_type == "Gumbel") {
    # Gumbel Copula-Based Lognormal-Pareto Simulation
    mu=exp(parameters[1]);sigma=exp(parameters[2]);alpha=exp(parameters[3]); 
    y_m=exp(parameters[4]);theta=exp(parameters[5])
    sim.output <- simulate_truncated_copula_LognormalPareto(N, mu, sigma, alpha,
                                                            y_m, "Gumbel", theta,truncation_threshold = 500000)
  } else if (copula_type == "Frank") {
    # Frank Copula-Based Lognormal-Pareto Simulation
    mu=exp(parameters[1]);sigma=exp(parameters[2]);alpha=exp(parameters[3]); 
    y_m=exp(parameters[4]);theta=exp(parameters[5])
    sim.output <- simulate_truncated_copula_LognormalPareto(N, mu, sigma, alpha,
                                                            y_m, "Frank", theta,truncation_threshold = 500000)
  } else if (copula_type == "t") {
    # t-Copula-Based Lognormal-Pareto Simulation
    mu=exp(parameters[1]);sigma=exp(parameters[2]);alpha=exp(parameters[3]); 
    y_m=exp(parameters[4]);rho=exp(parameters[5]);df_t=exp(parameters[6])
    sim.output <- simulate_truncated_copula_LognormalPareto(N, mu, sigma, alpha,
                                                            y_m, "t", rho, df_t,truncation_threshold = 500000)
  } 
  
  
  return(sim.output)
}




# A function to select the parameter levels of the simulator of interest
parameter_label_selection<- function(copula_type){
  
  # Selecting the simulator of interest
  if (copula_type == "Lognormal-Pareto"){
    # Independent Lognormal-Pareto
    parameter.labels<- c("mu", "sigma", "alpha", "y_m")
  }
  if (copula_type == "Gaussian") {
    # Gaussian Copula-Based Lognormal-Pareto Simulation
    parameter.labels<- c("mu", "sigma", "alpha", "y_m", "rho")
  } else if (copula_type == "Clayton") {
    # Clayton Copula-Based Lognormal-Pareto Simulation
    parameter.labels<- c("mu", "sigma", "alpha", "y_m", "theta")
  } else if (copula_type == "Gumbel") {
    # Gumbel Copula-Based Lognormal-Pareto Simulation
    parameter.labels<- c("mu", "sigma", "alpha", "y_m", "theta")
  } else if (copula_type == "Frank") {
    # Frank Copula-Based Lognormal-Pareto Simulation
    parameter.labels<- c("mu", "sigma", "alpha", "y_m", "theta")
  } else if (copula_type == "t") {
    # t-Copula-Based Lognormal-Pareto Simulation
    parameter.labels<- c("mu", "sigma", "alpha", "y_m", "rho","df_t")
  } 
  
  
  return(parameter.labels)
}