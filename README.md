# **GLeMa: Maximum Likelihood Estimation for Generalized Linear Models**

**GLeMa** is a two-part analytical and computational project focused on deriving and implementing **Maximum Likelihood Estimators (MLEs)** for **Generalized Linear Models (GLMs)**.
The project combines mathematical derivations with practical computation in **R**, and applies the concepts to a real-world dataset from Spotify and YouTube.

## **Analytical Component - MLE for a GLM**

The main goal is to understand and derive the MLE for parameters of a selected **Generalized Linear Model**.
The working example is a **Poisson Regression Model**, a GLM used for modeling count data (e.g., YouTube video views, likes, comments).

Key tasks:

* Introduce the GLM framework:

  * **Random component:** probability distribution of the response (Poisson for counts)
  * **Systematic component:** linear predictor
  * **Link function:** typically log link for Poisson GLMs

* Derive:

  * Poisson log-likelihood function
  * Score function and Hessian
  * Analytical form of the MLE for, if available, or the iterative procedure (e.g., Newton-Raphson)

* Evaluate whether the estimator is **biased or unbiased**.



## **Computational Component - Simulation and Estimation in R**

This part validates the analytical results using numerical experiments.

Steps:

* Simulate a dataset from a Poisson GLM.
* Implement the likelihood and optimization procedure in R.
* Compare:

  * Analytical MLE vs. simulated MLE
  * Estimated parameters vs. true parameters
