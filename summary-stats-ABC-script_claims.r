# Load required libraries
library(compiler)
library(moments)
library(ineq)
library("Dowd") # For Pickands estimator
library(DescTools)  # For Theil Index

# Function to compute the Moment Estimator for the largest k claims
moment_estimator <- function(sorted_claims, k) {
  if (length(sorted_claims) <= k) return(NA)  # Avoid issues if k is too large
  
  x_k1 <- sorted_claims[k + 1]  # (k+1)-th order statistic
  log_ratios <- log(sorted_claims[1:k] / x_k1)  # Logarithm of order statistic ratios
  
  M1 <- mean(log_ratios)  # First moment
  M2 <- mean(log_ratios^2)  # Second moment
  
  gamma_hat <- M1 + 1 - (0.5 * (1 - (M2 / M1^2))^(-1))  # Moment Estimator formula
  
  return(gamma_hat)
}

# Function to compute Expected Shortfall/Conditional Tail Expectation (CTE) at 95%
expected_shortfall <- function(x, q = 0.95) {
  threshold <- quantile(x, probs = q, na.rm = TRUE)
  excess_losses <- x[x > threshold]
  if (length(excess_losses) == 0) return(NA)
  return(mean(excess_losses, na.rm = TRUE))
}

# Function to compute Mean Excess Function Slope
mean_excess_slope <- function(x) {
  sorted_x <- sort(x, decreasing = TRUE)
  u <- quantile(sorted_x, probs = 0.90, na.rm = TRUE)  # 90th percentile as threshold
  excess_values <- sorted_x[sorted_x > u] - u
  if (length(excess_values) < 2) return(NA)
  return(mean(excess_values, na.rm = TRUE))
}

# Function to compute Pickands estimator (Robust Tail Index)
pickands_estimator <- function(x) {
  sorted_x <- sort(x, decreasing = TRUE)
  n <- length(sorted_x)
  if (n < 4) return(NA)  # Pickands requires at least 4 observations
  return(evd::pickands(sorted_x))
}

# Function to compute Theil Entropy Index
theil_index <- function(x) {
  return(DescTools::Theil(x, na.rm = TRUE))
}

# Updated summary statistics function (15 summaries)
summary_func <- function(data) {
  x <- data[, "Claim_Severity"]
  
  # Basic Distributional Statistics
  mean_value <- mean(x, na.rm = TRUE)
  sd_value <- sd(x, na.rm = TRUE)
  median_value <- median(x, na.rm = TRUE)
  Q1_value <- as.numeric(quantile(x, probs = 0.25, na.rm = TRUE))
  Q3_value <- as.numeric(quantile(x, probs = 0.75, na.rm = TRUE))
  IQR_value <- Q3_value - Q1_value
  skewness_value <- moments::skewness(x, na.rm = TRUE)
  kurtosis_value <- moments::kurtosis(x, na.rm = TRUE)
  
  # Inequality Measure
  gini_index <- ineq::Gini(x, na.rm = TRUE)
  
  # Extreme Quantiles
  extreme_quantile_99 <- as.numeric(quantile(x, probs = 0.99, na.rm = TRUE))
  
  # Expected Shortfall (CTE at 95%)
  CTE_95 <- expected_shortfall(x, q = 0.95)
  
  # Mean Excess Function Slope
  MEF_slope <- mean_excess_slope(x)
  
  # Pickands Estimator for Tail Index
  #Estimates the Value of Pickands Estimator for a specified data set and chosen tail size. 
  #Notes: (1) We estimate the Pickands Estimator by looking at the upper tail. 
  #(2) The tail size must be less than one quarter of the total sample size. (3) The tail size must be a scalar.
  PickandsEstimator_func<- function(x){
    #Estimates the Value of Pickands Estimator for a specified data set and chosen tail size. 
    #Notes: (1) We estimate the Pickands Estimator by looking at the upper tail. 
    #(2) The tail size must be less than one quarter of the total sample size. (3) The tail size must be a scalar.
    n <- length(x)
    k_pickands <- max(10, round(sqrt(n) / 2))  # Choose a reasonable k
    PickandsEstimator_est<- Dowd::PickandsEstimator(Ra=x,tail.size=k_pickands)
    return(PickandsEstimator_est)
  }
  
  tail_index_pickands <- PickandsEstimator_func(x)
  
  
  
  
  # Moment Estimator (adaptive k)
  sorted_claims <- sort(x, decreasing = TRUE)
  n <- length(sorted_claims)
  k_chosen <- max(2, round(sqrt(n)))  # Adaptive choice of k
  moment_gamma <- moment_estimator(sorted_claims, k_chosen)
  
  # Entropy Measure (Theil Index)
  # Theil Index using Entropy function
  theil_index <- function(x) {
    return(DescTools::Entropy(x, base = exp(1), na.rm = TRUE))
  }
  theil_entropy <- theil_index(x)
  

  
  # Return all statistics (a total of 15 summaries)
  return(c(mean_value, sd_value, median_value, Q1_value, Q3_value, IQR_value,
           skewness_value, kurtosis_value, gini_index, extreme_quantile_99,
           CTE_95, MEF_slope, tail_index_pickands, moment_gamma, theil_entropy))
}


summary_compiler=compiler::cmpfun(summary_func)