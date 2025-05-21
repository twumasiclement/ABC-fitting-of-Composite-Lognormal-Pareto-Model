# Paper Title: Likelihood-Free Bayesian Estimation of Non-Life Insurance Claim Severity: A Sequential Monte Carlo ABC of the Truncated Composite Lognormal-Pareto Model (LNP)

Robust calibration of insurance claim severity models remains a longstanding challenge due to the complexity or intractability of likelihood functions in classical maximum likelihood estimation procedures. These likelihood-based methods often require computationally intensive optimisations, are sensitive to initial parameter values and prone to convergence issues, especially in high-dimensional or poorly conditioned models. To address these limitations, this study introduces a likelihood-free calibration framework using approximate Bayesian computation (ABC). Specifically, we propose a modified sequential Monte Carlo ABC (ABC-SMC) algorithm—originally developed for a high-dimensional biological system—and adapt it for the first time to non-life insurance claim severity modelling. Using a truncated four-parameter Lognormal–Pareto (LNP) composite model as a case study, we demonstrate how the proposed framework enables efficient, likelihood-free Bayesian inference for heavy-tailed insurance loss data. Comparative simulation experiments demonstrated that the calibration performance of the ABC methods is contingent upon the model and data utilised. The ABC-SMC estimator was statistically unbiased across all parameters and varying Monte Carlo sample sizes, confirming its fidelity. Although Ridge adjustment enhanced posterior variance and mean squared error in these experiments, Lasso adjustment provided superior accuracy and empirical data coverage when fitting the LNP model to the extensively analysed Norwegian fire insurance severity dataset over time. This study contributes a flexible and scalable inferential tool for actuarial applications, enabling the calibration of sophisticated composite and mixture models without explicit likelihoods. It lays the foundation for broader adoption of robust, likelihood-free Bayesian methods in insurance modelling and paves the way for future extensions.






# Below are the labels of the main R codes/scripts for implementing the Sequential Monte Carlo ABC (ABC-SMC) and the two penalised regression-adjusted ABC methods (based on Ridge and Lasso regularisation procedures):

1. `Model-and-Parameter-Label-Selection-script.R`: R codes of the composite LNP model (designed to also apply copula-based claim severity models)

2. `summary-stats-ABC-script_claims.r`: R Codes for computing the ABC summary statistics given a claim severity data

3. `ABC-Importance-Sampling-Improved-script_claims.R`: R codes for implementing ABC-SMC sequential importance sampling (designed to also apply copula-based claim severity models)

4. `distance-script_claims.R`: Codes for computing the Kullback-Leibler divergence-based distance metric to quantify the discrepancy between observed and simulated data 

5. `Weighted-iterative-ABC-script_claims.R`:  Main R codes for implementing the Modified ABC-SMC algorithm for fitting claim severity models (designed to also apply copula-based claim severity models)

6. `Post-Ridge-reg-adj-L2-script_claims.R`: Codes for the ridge-adjusted ABC posterior correction (L2 regularisation)

7. `Post-Lasso-reg-adj-L1-script_claims.R`: Codes for the lasso-adjusted ABC posterior correction (L1 regularisation)
