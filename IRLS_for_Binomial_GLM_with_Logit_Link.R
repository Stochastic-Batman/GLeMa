# for reproducibility
set.seed(95)

# IRLS for Binomial GLM with Logit Link 
# Check out the derivation for this algorithm in "GLeMa documentation.pdf"
IRLS_Binomial <- function(X, y, n, eps = 1e-8, M = 100) {
  # X: design matrix (N x (p + 1))
  # y: response vector (counts)
  # n: number of trials for each observation
  # eps: convergence tolerance
  # M: maximum number of iterations
  
  N <- nrow(X)
  p <- ncol(X) - 1
  
  beta <- rep(0, p + 1)
  
  for (r in 1:M) {
    # Step 1: Compute linear predictors
    # I had not used %*% before, but it turned out to be R's built-in matrix multiplication
    # The pseudocode I derived uses for loop, but it can be simplified to the matrix multiplication
    eta <- X %*% beta 
    
    # Step 2: Compute fitted probabilities
    p_r <- exp(eta) / (1 + exp(eta))
    mu <- n * p_r
    
    # Step 3: Compute working responses
    z <- eta + (y - mu) / (n * p_r * (1 - p_r))
    
    # Step 4: Compute iterative weights
    w <- n * p_r * (1 - p_r)
    W <- diag(as.vector(w))
    
    # Step 5: Weighted least squares update
    # t(X) is transpose: X^T
    # solve() finds the inverse by solving a %*% x = b; if b is missing, as it is here, it is assuemd to be identity matrix
    beta_new <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% z
    
    # Step 6: Check convergence
    delta <- sqrt(sum((beta_new - beta)^2))
    beta <- as.vector(beta_new)
    
    if (delta < eps) {
      cat("Converged in", r, "iteiterations\n")
      break
    }
  }
  
  if (r == M) {
    cat("Warning: Maximum iterations reached\n")
  }
  
  return(list(beta = beta, iterations = r))
}


# Generating synthetic data with random seed 95 for reproducibility
N <- 500  # number of observations
p <- 2    # number of predictors (excluding intercept)

# True parameter values (known)
beta_true <- c(1.5, -2.7, 3.14)  # in this order: intercept, beta1, beta2

# Generate predictor variables
X1 <- rnorm(N, mean = 0, sd = 1)
X2 <- rnorm(N, mean = 0, sd = 1)

# Create design matrix (add intercept column)
X <- cbind(1, X1, X2)

# Generate number of trials (randomly between 1 and 20)
n_trials <- sample(1:20, N, replace = TRUE)

# Compute true linear predictor
eta_true <- X %*% beta_true

# Compute true probabilities (inverse logit)
p_true <- exp(eta_true) / (1 + exp(eta_true))

# Generate binomial responses
y <- rbinom(N, size = n_trials, prob = p_true)

# Fit model using my IRLS binomial specific algorithm
cat("Fitting with IRLS\n")
result <- IRLS_Binomial(X, y, n_trials)

# Compare with R's built-in glm function
cat("\nFitting with R's glm\n")
glm_fit <- glm(cbind(y, n_trials - y) ~ X1 + X2, family = binomial(link = "logit"))

cat("\nLet's compare!\n")
cat("\nTrue Parameters:\n")
print(sprintf("%.4f", beta_true))

cat("\nIRLS Estimates:\n")
print(sprintf("%.4f", result$beta))

cat("\nR's glm Estimates:\n")
print(sprintf("%.4f", coef(glm_fit)))

cat("\nDifference (IRLS - True):\n")
print(sprintf("%.4f", result$beta - beta_true))

cat("\nDifference (glm - True):\n")
print(sprintf("%.4f", coef(glm_fit) - beta_true))

cat("\nDifference (IRLS - glm):\n")
print(sprintf("%.4f", result$beta - coef(glm_fit)))

cat("\nSummary\n")
cat("Root Mean Squared Error (IRLS):", sprintf("%.4f", sqrt(mean((result$beta - beta_true)^2))), "\n")
cat("Root Mean Squared Error (glm):", sprintf("%.4f", sqrt(mean((coef(glm_fit) - beta_true)^2))), "\n")
cat("Max Absolute Difference (IRLS vs glm):", sprintf("%.4f", max(abs(result$beta - coef(glm_fit)))), "\n")

fitted_probs <- exp(X %*% result$beta) / (1 + exp(X %*% result$beta))

png("true_vs_fitted_probabilities.png", width = 800, height = 600)
plot(p_true, fitted_probs,xlab = "True Probabilities", ylab = "Fitted Probabilities (IRLS)", main = "True vs Fitted Probabilities", pch = 16, col = rgb(0, 0, 1, 0.3))
abline(0, 1, col = "red", lwd = 2)
dev.off()