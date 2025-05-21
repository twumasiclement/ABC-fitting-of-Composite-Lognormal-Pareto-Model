require("kedd")

#Gaussian kernel with bandwidth delta
gauss_kernel <- function(dist, delta) {
  kern <- (1 / (sqrt(2 * pi) * delta)) * exp(- (dist^2) / (2 * delta^2))
  return(kern)
}



#R codes for Ridge regression-adjusted ABC for resulting posterior correction
Post_Ridge_reg_adj<- function(post_distn,summary_obs,copula_type=copula_type,N) {
  #k #biasing parameter or ridge regression penalty parameter
  # post_dtn is the posterior sample
  # w are weights across summary statistics for computing weighted distance
  #storing the summary stats across simulation realizations
  S_i <- NULL #saving saving stats for each simulation realizaton 
  no_of_parameters<- number_of_parameters_func(copula_type) #number of parameters to be estimated
  posterior_mean_adj<- rep(NA,no_of_parameters) #storing adjusted posterior means
  
  m<- dim(post_distn)[1]# m=number of posterior samples
  d<- rep(0, m)# weighted distances given observed data
  p<- length(summary_obs) #dimension of summary statistics
  Unadj_dist<- post_distn
  SummaryStats_sim <- NULL
  X_Design_matrix<- matrix(NA, ncol=p,nrow=m) #design matrix
  # Weights based on Gaussian kernel for local-linear regression adjustment
  W<- matrix(0, ncol=m,nrow=m)
  X_bar=numeric(length=p) #saving weighted column means of design matrix 
  beta_lasso<-list()#store regression coefficients for each dependent variables
  
  for (i in 1:m) {
    theta<- as.vector(unlist(post_distn[i,]))
    
    output_sim<- simulation_model_selection(N=N,copula_type, parameters=theta)
    
    #Computing the summary stats for each simulation realisation
    options(warn=-1)
    SummaryStats_sim[[i]] <- summary_compiler(data=output_sim)
    
    #Storing the summary stats and weighted distances (between summaries of observed and simulated data)
    SummaryStats_sim_combined<- SummaryStats_sim[[i]]
    
    
    diff<- SummaryStats_sim_combined-summary_obs
    X_Design_matrix[i, ]<- diff  #storing each row of design matrix X
    # Computing weights based on
    #Storing the average summary stats and distance (between summaries of observed and simulated data)
    S_i[[i]] <- SummaryStats_sim_combined
    # w <-apply(S_i[[i]], 2, var, na.rm = TRUE)
    #w<- w/sum(w) #normalised weight for computing distance
    #d[i] <- w_distance(S1=S_i[[i]], S2=summary_obs, weight=w)
    
    #summaries_obs= summary statistics of the observed/pseudo-observed data
    d[i]<- KL_distance_compiler(S_i[[i]],summary_obs)#KL distance
  }
  
  distances<-na.inf.zero(d)
  #Adaptively choosing the bandwidth of the Gaussian kernel based on the distances
  #bandwidth<- KernSmooth::dpik(x=distances, "normal",scalest = "minim")
  bandwidth<- kedd::h.amise(x=distances, deriv.order =0,kernel = c("gaussian"))$h
  diag(W)<- gauss_kernel(dist= distances,delta=bandwidth)
  theta_post<- as.matrix(post_distn)  
  weights<- diag(W)/sum(diag(W)) # (normalising) main diagonal of Weighting matrix
  
  
  #Transforming X and Y (posterior distribution and summary statistics)
  for(j in seq_along(posterior_mean_adj)){#For each jth model parameter, j=1,2,...23
    X<- X_Design_matrix
    Y<- theta_post[ ,j]
    
    #Step 1 (Centering X and Y)
    for (k in 1:p) X_bar[k]<- sum(weights*X[,k])
    
    for (k in 1:p) {
      X[, k]<- X[, k]-X_bar[k]
    }
    #finding the weighted mean of Y and centring
    Y_bar <- sum(weights*Y)
    Y<- Y- Y_bar
   
    
    #Step 2: scaling (centred X and Y) by weights
    
    for(k in 1:p) X[, k]<- sqrt(weights)*X[, k]
    Y<- sqrt(weights)*Y
    
    #Choose optimal value of k (the penalty paramters)
    # Using cross validation glmnet
    # Setting the range of lambda values
    options(warn = -1)
    #lambda_seq <- 10^seq(2, -2, by = -.1)
    ridge_cv <- cv.glmnet(X, Y, alpha = 0)#lambda =lambda_seq)
    # Best lambda value
    best_lambda <- ridge_cv$lambda.min
    k<-best_lambda
    #print(k)
    
    
    
    # calculate beta estimates corresponding to summary statistics X (standardised coefficients)
    beta_ridge <- solve(t(X) %*%W%*%X+ k*diag(p)) %*% t(X)%*%W%*%Y
    
    
    # calculate intercept (adjusted posterior mean estimates)
    posterior_mean_adj[j] = exp(Y_bar - X_bar%*%beta_ridge)
    
    #Adjusting the posterior distribution
    Unadj_dist[,j]<- post_distn[, j]- X_Design_matrix%*%beta_ridge
  }
  
  
  
  posterior_mean_unadj<- apply(exp(post_distn),2,mean)
  Posterior_mean_output<- data.frame(Adj_posterior_mean=posterior_mean_adj,
                                     Uadj_posterior_mean=posterior_mean_unadj)
  
  
  rownames(Posterior_mean_output)<- parameter_label_selection(copula_type)
  
  
  return(list(X_Design_matrix=X,Posterior_mean_output=Posterior_mean_output,Adjusted_posterior_dist=Unadj_dist))
}