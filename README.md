# **GLeMa**

**Maximum Likelihood Estimators (MLEs)** for **Generalized Linear Models (GLMs)**. The project combines mathematical derivations with practical computation in **R**.

## **Analytical Component - MLE for a GLM**

The main goal is to understand and derive the MLE for parameters of a selected **Generalized Linear Model**.
The working example is a **Binomial Regression Model**, a GLM used for modeling count data (e.g., YouTube video views, likes, comments).

Key tasks:

* Introduce the GLM framework:

  * **Random component:** probability distribution of the response (Poisson for counts)
  * **Systematic component:** linear predictor
  * **Link function:** logit link for Binomial GLMs

* Derive:

  * Binomial likelihood function
  * Score function

* Analytical form of the MLE for, if available, or the iterative procedure (spoiler: I derived iterative procedure)

Evaluate whether the estimator is **biased or unbiased**.



## Computational Component - Simulation and Estimation in R

This section details the implementation of an **Iteratively Reweighted Least Squares (IRLS)** algorithm for estimating coefficients in a **Binomial Generalized Linear Model (GLM)** with a logit link. 

### What's it about?

- **Reproducibility**: The random seed is set for consistent results.
- **IRLS Function**: The core function estimates coefficients iteratively until convergence, using a specified maximum number of iterations.
- **Data Generation**: Synthetic binomial data is created based on true parameters, along with predictor variables.
- **Model Fitting**: The IRLS method is applied to the simulated data, with results compared to R's built-in **`glm()`** output.
- **Performance Evaluation**: Differences in estimates and root mean squared error (RMSE) are computed to assess accuracy.
- **Visualization**: A plot shows the relationship between true and fitted probabilities, further validating the model's performance.
