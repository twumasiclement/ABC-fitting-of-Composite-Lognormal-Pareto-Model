# Paper Title: Likelihood-Free Bayesian Estimation of Non-Life Insurance Claim Severity: A Sequential Monte Carlo ABC of the Truncated Composite Lognormal-Pareto Model (LNP)







# Below are the labels of the main R codes/scripts for implementing the Sequential Monte Carlo ABC (ABC-SMC) and the two penalised regression-adjusted ABC methods (based on Ridge and Lasso regularisation procedures):

1. `Model-and-Parameter-Label-Selection-script.R`: R codes of the composite LNP model (designed to also apply copula-based claim severity models)

2. `summary-stats-ABC-script_claims.r`: R Codes for computing the ABC summary statistics given a claim severity data

3. `ABC-Importance-Sampling-Improved-script_claims.R`: R codes for implementing ABC-SMC sequential importance sampling (designed to also apply copula-based claim severity models)

4. `distance-script_claims.R`: Codes for computing the Kullback-Leibler divergence-based distance metric to quantify the discrepancy between observed and simulated data 

5. `Weighted-iterative-ABC-script_claims.R`:  Main R codes for implementing the Modified ABC-SMC algorithm for fitting claim severity models (designed to also apply copula-based claim severity models)

6. `Post-Ridge-reg-adj-L2-script_claims.R`: Codes for the ridge-adjusted ABC posterior correction (L2 regularisation)

7. `Post-Lasso-reg-adj-L1-script_claims.R`: Codes for the lasso-adjusted ABC posterior correction (L1 regularisation)
