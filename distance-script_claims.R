# Function for computing weighted distance between simulated and observed/pseudo-observed summary statistics
#Kullback–Leibler (KL) divergence
#In statistics, the Kullback–Leibler (KL) divergence is a distance metric that 
#quantifies the difference between two probability distributions.
KL_distance <- function(S1, S2)  {
  #weight is a vector/one-dimensional array (summary stats weights)
  normalised_func<- function(x) na.omit(x)/sum(x,na.rm=TRUE)
  
  # normalise the summaries to sum to unity
  P1=na.omit(normalised_func(S1))
  Q1=na.omit(normalised_func(S2))
  
  # print(paste("Length of S1:"))
  # print(length(P1))
  
  #print(paste("Length of S2:"))
  #print(length(Q1))
  
  #rbind distributions into one matrix
  x1<-rbind(P1,Q1)
  #calculate KL divergence
  KL_distance<- as.numeric(suppressMessages(philentropy::KL(x1, unit='log')))
  
  return(KL_distance)#return Kullback–Leibler (KL) divergence distance 
}


KL_distance_compiler=cmpfun(KL_distance)










