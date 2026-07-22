# ============================================================================
# FILE 3: code/01_simulation_NK2.R 
# ============================================================================
# Monte Carlo simulation for Table 4 (NK2: mu=2.0, sigma_u=0.5)
# R = 1000 for n = 100, 200, and 1000 (as in manuscript)
# ============================================================================
# ----------------------------------------------------------------------------
# Function to run simulation for a given sample size
# ----------------------------------------------------------------------------
run_simulation_NK2 <- function(n, R, 
                               true_beta0 = 1.0, true_beta1 = 0.6, true_beta2 = 0.4,
                               true_sigma_v = 0.3, true_mu = 2.0, true_sigma_u = 0.5,
                               true_TE = 0.634) {
  
  results_beta0 <- c()
  results_beta1 <- c()
  results_beta2 <- c()
  results_sigma_v <- c()
  results_sigma_u <- c()
  results_mu <- c()
  results_mean_TE <- c()
  
  for (r in 1:R) {
    set.seed(12345 + r)
    
    # Generate data
    x1 <- runif(n, 1, 10)
    x2 <- rnorm(n, 5, 2)
    v <- rnorm(n, 0, true_sigma_v)
    u <- sqrt(qgamma(runif(n), shape = true_mu, rate = true_mu/true_sigma_u^2))
    y <- true_beta0 + true_beta1*x1 + true_beta2*x2 + v - u
    X_mat <- cbind(1, x1, x2)
    
    # Initial parameters
    init_params <- c(mean(y), 0.6, 0.4, log(0.3), log(0.5), log(1.8))
    
    # Optimization
    opt_result <- tryCatch({
      optim(init_params, neg_loglik_nakagami, y = y, X = X_mat,
            method = "Nelder-Mead", 
            control = list(maxit = 1000, trace = 0, reltol = 1e-8))
    }, error = function(e) NULL)
    
    if (!is.null(opt_result) && opt_result$convergence == 0) {
      params_opt <- opt_result$par
      K <- ncol(X_mat)
      
      results_beta0 <- c(results_beta0, params_opt[1])
      results_beta1 <- c(results_beta1, params_opt[2])
      results_beta2 <- c(results_beta2, params_opt[3])
      results_sigma_v <- c(results_sigma_v, exp(params_opt[K+1]))
      results_sigma_u <- c(results_sigma_u, exp(params_opt[K+2]))
      results_mu <- c(results_mu, exp(params_opt[K+3]))
      
      # Compute mean technical efficiency
      TE_i <- numeric(n)
      for (i in 1:n) {
        eps_i <- y[i] - sum(X_mat[i, ] * params_opt[1:3])
        TE_i[i] <- compute_TE(eps_i, exp(params_opt[K+1]), exp(params_opt[K+3]), 
                              exp(params_opt[K+2]), params_opt[1:3], X_mat[i, ])
      }
      results_mean_TE <- c(results_mean_TE, mean(TE_i, na.rm = TRUE))
    }
    
    if (r %% 200 == 0) cat("Replication", r, "of", R, "completed. Successes:", length(results_beta0), "\n")
  }
  
  # Compute summary statistics
  n_success <- length(results_beta0)
  if (n_success == 0) {
    cat("ERROR: No successful replications!\n")
    return(NULL)
  }
  
  bias_beta0 <- mean(results_beta0, na.rm = TRUE) - true_beta0
  bias_beta1 <- mean(results_beta1, na.rm = TRUE) - true_beta1
  bias_beta2 <- mean(results_beta2, na.rm = TRUE) - true_beta2
  bias_sigma_v <- mean(results_sigma_v, na.rm = TRUE) - true_sigma_v
  bias_sigma_u <- mean(results_sigma_u, na.rm = TRUE) - true_sigma_u
  bias_mu <- mean(results_mu, na.rm = TRUE) - true_mu
  
  rmse_beta0 <- sqrt(mean((results_beta0 - true_beta0)^2, na.rm = TRUE))
  rmse_beta1 <- sqrt(mean((results_beta1 - true_beta1)^2, na.rm = TRUE))
  rmse_beta2 <- sqrt(mean((results_beta2 - true_beta2)^2, na.rm = TRUE))
  rmse_sigma_v <- sqrt(mean((results_sigma_v - true_sigma_v)^2, na.rm = TRUE))
  rmse_sigma_u <- sqrt(mean((results_sigma_u - true_sigma_u)^2, na.rm = TRUE))
  rmse_mu <- sqrt(mean((results_mu - true_mu)^2, na.rm = TRUE))
  
  mae_beta0 <- mean(abs(results_beta0 - true_beta0), na.rm = TRUE)
  mae_beta1 <- mean(abs(results_beta1 - true_beta1), na.rm = TRUE)
  mae_beta2 <- mean(abs(results_beta2 - true_beta2), na.rm = TRUE)
  
  mean_TE_est <- mean(results_mean_TE, na.rm = TRUE)
  rmse_TE <- sqrt(mean((results_mean_TE - true_TE)^2, na.rm = TRUE))
  
  # Output results
  cat("\n========== Simulation Results n =", n, "==========\n")
  cat("Successful replications:", n_success, "of", R, "\n\n")
  cat("--- Bias ---\n")
  cat(sprintf("beta0: %.4f\n", bias_beta0))
  cat(sprintf("beta1: %.4f\n", bias_beta1))
  cat(sprintf("beta2: %.4f\n", bias_beta2))
  cat(sprintf("sigma_v: %.4f\n", bias_sigma_v))
  cat(sprintf("sigma_u: %.4f\n", bias_sigma_u))
  cat(sprintf("mu: %.4f\n", bias_mu))
  
  cat("\n--- RMSE ---\n")
  cat(sprintf("beta0: %.4f\n", rmse_beta0))
  cat(sprintf("beta1: %.4f\n", rmse_beta1))
  cat(sprintf("beta2: %.4f\n", rmse_beta2))
  cat(sprintf("sigma_v: %.4f\n", rmse_sigma_v))
  cat(sprintf("sigma_u: %.4f\n", rmse_sigma_u))
  cat(sprintf("mu: %.4f\n", rmse_mu))
  
  cat("\n--- MAE (for betas) ---\n")
  cat(sprintf("beta0: %.4f\n", mae_beta0))
  cat(sprintf("beta1: %.4f\n", mae_beta1))
  cat(sprintf("beta2: %.4f\n", mae_beta2))
  
  cat("\n--- Technical Efficiency ---\n")
  cat(sprintf("Estimated Mean TE: %.3f\n", mean_TE_est))
  cat(sprintf("RMSE(TE): %.4f\n", rmse_TE))
  
  # Return results as data frame
  results_df <- data.frame(
    Parameter = c("beta0", "beta1", "beta2", "sigma_v", "sigma_u", "mu", "Mean_TE"),
    Bias = c(bias_beta0, bias_beta1, bias_beta2, bias_sigma_v, bias_sigma_u, bias_mu, NA),
    RMSE = c(rmse_beta0, rmse_beta1, rmse_beta2, rmse_sigma_v, rmse_sigma_u, rmse_mu, rmse_TE),
    MAE = c(mae_beta0, mae_beta1, mae_beta2, NA, NA, NA, NA)
  )
  
  return(results_df)
}

# ----------------------------------------------------------------------------
# Run simulations for n = 100, 200, 1000 (R = 1000 for all as in manuscript)
# ----------------------------------------------------------------------------

cat("\n========================================\n")
cat("Table 3: NK2 Simulation (mu=2.0, sigma_u=0.5)\n")
cat("========================================\n")

# Create outputs directory if it doesn't exist
if (!dir.exists("outputs")) dir.create("outputs")

cat("\nRunning simulation for n = 100 with R = 1000...\n")
results_n100 <- run_simulation_NK2(n = 100, R = 1000)
if (!is.null(results_n100)) write.csv(results_n100, "outputs/table3_n100.csv", row.names = FALSE)

cat("\nRunning simulation for n = 200 with R = 1000...\n")
results_n200 <- run_simulation_NK2(n = 200, R = 1000)
if (!is.null(results_n200)) write.csv(results_n200, "outputs/table3_n200.csv", row.names = FALSE)

cat("\nRunning simulation for n = 1000 with R = 1000...\n")
results_n1000 <- run_simulation_NK2(n = 1000, R = 1000)
if (!is.null(results_n1000)) write.csv(results_n1000, "outputs/table3_n1000.csv", row.names = FALSE)

cat("\nTable 3 simulation completed successfully.\n")
