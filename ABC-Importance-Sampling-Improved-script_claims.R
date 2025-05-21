
#A function to determine the prior based on the type of model
# Initial Prior and posterior distribution function
#prior distribution of model parameters (log scale)
prior<- function() {
  # Prior distributions on log scales
  sigma_parameter<- runif(1, min=-2 , max=0.45)
  alpha_parameter<- runif(1, min=-1 , max=1)
  y_m_parameter<- runif(1, min=13.8 , max=14.9)
  #mu <- log(y_m)-(alpha*sigma^2)
  mu_parameter_antilog<- abs(log(exp(y_m_parameter))-(exp(alpha_parameter)*(exp(sigma_parameter)^2)))
  #On log scale for the lognormal meanlog parameter
  mu_parameter<- log(mu_parameter_antilog)
  
  #rho_parameter<- runif(1, min=-2 , max=-0.145)
  #if(rho_parameter>0) rho_parameter=0 #correlation cannot exceed 1
 # theta_parameter<- runif(1, min=0 , max=1.09)
  #df_t_parameter<- runif(1, min=-1 , max=1.72)
  # return a realisation of a prior sample
  
  if(copula_type=="Lognormal-Pareto"){
    return(c(mu_parameter,sigma_parameter,alpha_parameter,y_m_parameter))
  }
  else if(copula_type=="Gaussian"){
    return(c(mu_parameter,sigma_parameter,alpha_parameter,y_m_parameter, rho_parameter))
  }
  else if (copula_type == "Clayton"|copula_type == "Gumbel"|copula_type == "Frank"){
    return(c(mu_parameter,sigma_parameter,alpha_parameter,y_m_parameter,theta_parameter))
    
  }
  else if(copula_type=="t"){
    return(c(mu_parameter,sigma_parameter,alpha_parameter,y_m_parameter,
             rho_parameter, df_t_parameter))
  } 
}





#A function to determine the number of parameters to be estimated based on the type of model
number_of_parameters_func<- function(copula_type){
  # Selecting the simulator of interest
  if (copula_type == "Lognormal-Pareto"){
    number_of_parameters<- 4 #number of parameters to be estimated
  }
  else if (copula_type == "Gaussian") {
    # Gaussian Copula-Based Lognormal-Pareto Simulation
    number_of_parameters<- 5 #number of parameters to be estimated
    
  } else if (copula_type == "Clayton"|copula_type == "Gumbel"|copula_type == "Frank") {
    # Clayton Copula-Based Lognormal-Pareto Simulation
    number_of_parameters<- 5 #number of parameters to be estimated
  } else if (copula_type == "t") {
    # t-Copula-Based Lognormal-Pareto Simulation
    number_of_parameters<- 6 #number of parameters to be estimated
  } 
  
  return(number_of_parameters)
}



#Weighted-iterative ABC via
#Sequential Monte Carlo with Importance sampling
#Function for ABC calibration
ABC <- function(fork, pftn , n) {
  # pftn is prior function
  # n is number of samples
  #dim_summaries= A total of 15 summary statistics
  
  dimS<- dim_summaries #dimension or number of ABC summary statistics
  
  
  number_of_parameters<- number_of_parameters_func(copula_type) #number of parameters to be estimated
  
  
  theta  <- matrix(nrow = n, ncol = number_of_parameters)# matrix of prior distributions
  theta_main<- list()
  #storing the summary stats across all simulation realisations
  S_i <- NULL #saving saving stats of each simulation realisaton     
  #S is a matrix(nrow = n, ncol = dimS) 
  d <- list()# weighted K-L distance
  SummaryStats_sim <- NULL
  output_without_errors<- list()

  
  for (i in 1:n) {
    print(paste("Simulation iteration=",i))
    theta[i, ] <- pftn()   
    output<- try(simulation_model_selection(N=N,copula_type, parameters=theta[i, ] ) 
      , silent = TRUE)
    
    #extract only simulations without errors
    if (!inherits(output, "try-error")) {
      k<- length(output_without_errors) + 1 #starting index (k=1) for main outputs
      theta_main[[k]]<- theta[k, ] 
      output_without_errors[[k]] <- output
    
    
    output_main<- output_without_errors[[k]]
    
    
    #Computing the summary stats for each simulation realisation
    SummaryStats_sim[[k]] <- summary_compiler(data=output_main)
    
    
    #Storing the summary stats and weighted distances (between summaries of pseudo-observed and simulated data)
    S_i[[k]] <- SummaryStats_sim[[k]]
    
    #print(paste("Simulated summary stats at","","i=", k))
    #print( S_i[[k]])
    
    #summaries_obs= summary statistics of the observed/pseudo-observed data
    d[[k]]<- KL_distance_compiler(S_i[[k]],summaries_obs)#KL distance

         }
    if (i %% 100 == 0) cat("fork", fork, "sample", i, "/", n, "\n")
  }
  
  
  S<-do.call("rbind",as.list(S_i))# summary stats matrix(nrow = n_complete, ncol = dimS) 
 
  return(list(theta= do.call("rbind",as.list(theta_main)), S=S, d=unlist(d)))
}

#Posterior function (draw from a new proposal distribution)
###Based on a multivariate normal kernel density


#Multivariate Normal kernel function given optimal bandwidth
#For peturbation
MultivNorm_rkernel<- function(Num,bandwidth_matrix){
  dim_k<- dim(bandwidth_matrix)[2]
  mean_vector<- rep(0,dim_k)
  return(tmvtnorm::rtmvnorm(n=1, mean=mean_vector, sigma=bandwidth_matrix, 
                lower=rep(-.2,dim_k),upper=rep(.2,dim_k), algorithm=c("gibbs")))
}





#Posterior function (draw from a new proposal distribution)
#Posterior function (draw from a new proposal distribution)
post <- function(samp=tha_post,importance_weight=weight,
                 optimal_bw_matrix=Sigma_optimal_t) {
  # new proposal from samp (previous prior samples)
  n <- dim(samp)[1]#dimension of proposal samples at t-1
  #randomly sampling from 1:n different sets of parameter values (importance sampling)
  #from previously accepted particles
  sample.particle<-sample(n, 1,prob=importance_weight)
  
  # Perturbation kernel: perturbing the particles for a new proposal  
  KDE_sampler<- samp[sample.particle, ]+MultivNorm_rkernel(Num=1,
                                          bandwidth_matrix=optimal_bw_matrix) 
  
  new_proposal<- KDE_sampler
  x<- new_proposal
  #Imposing this contraints in the peturbed particle: 
  #mu_parameter_antilog<- log(exp(y_m_parameter))-(exp(alpha_parameter)*(exp(sigma_parameter)^2)).
  #On log scale for the lognormal meanlog parameter 
  x[1] <- log(  abs(log(exp(x[4])) - (exp(x[3]) * (exp(x[2])^2)) ))

  
  
  return(x)
}





na.inf.zero<- function(x){
  x[is.na(x)|is.finite(x)==FALSE]<- 0
  return(x)
}


na.zero<- function(x){
  x[is.na(x)]<-0
  return(x)
}
