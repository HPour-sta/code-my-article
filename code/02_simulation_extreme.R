% ============================================================================
% FILE 4: code/02_simulation_extreme.R
% ============================================================================
# Monte Carlo simulation for Table 5 (Extreme Nakagami cases)
# NK-Low: mu=0.6, sigma_u=0.5
# NK-High: mu=5.0, sigma_u=1.0
# R = 1000 as in manuscript
# ----------------------------------------------------------------------------
# Function to run simulation for extreme cases
# ----------------------------------------------------------------------------
run_simulation_extreme <- function(n = 200, R = 1000, 
                                   true_beta0 = 1.0, true_beta1 = 0.6, true_beta2 = 0.4,
                                   true_sigma_v = 0.3, true_mu, true_sigma_u
                                  ) {
  
  results_beta0 <- c()
  results_beta1 <- c()
  results_beta2 <- c()
  results_mean_TE <- c()
  
set.seed(999)
  u_sim <- sqrt(qgamma(runif(50000), shape = true_mu, rate = true_mu/true_sigma_u^2))
  true_TE <- mean(exp(-u_sim))

  for (r in 1:R) {
    set.seed(23456 + r)
    
    # Generate data
    x1 <- runif(n, 1, 10)
    x2 <- rnorm(n, 5, 2)
    v <- rnorm(n, 0, true_sigma_v)
    u <- sqrt(qgamma(runif(n), shape = true_mu, rate = true_mu/true_sigma_u^2))
    y <- true_beta0 + true_beta1*x1 + true_beta2*x2 + v - u
    X_mat <- cbind(1, x1, x2)
    
    # Initial parameters
    init_params <- c(mean(y), 0.6, 0.4, log(0.3), log(true_sigma_u), log(true_mu))
    
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
  
  n_success <- length(results_beta0)
  if (n_success == 0) {
    cat("ERROR: No successful replications!\n")
    return(NULL)
  }
  
  # Compute summary statistics
  bias_beta0 <- mean(results_beta0, na.rm = TRUE) - true_beta0
  bias_beta1 <- mean(results_beta1, na.rm = TRUE) - true_beta1
  bias_beta2 <- mean(results_beta2, na.rm = TRUE) - true_beta2
  
  rmse_beta0 <- sqrt(mean((results_beta0 - true_beta0)^2, na.rm = TRUE))
  rmse_beta1 <- sqrt(mean((results_beta1 - true_beta1)^2, na.rm = TRUE))
  rmse_beta2 <- sqrt(mean((results_beta2 - true_beta2)^2, na.rm = TRUE))
  
  mae_beta0 <- mean(abs(results_beta0 - true_beta0), na.rm = TRUE)
  mae_beta1 <- mean(abs(results_beta1 - true_beta1), na.rm = TRUE)
  mae_beta2 <- mean(abs(results_beta2 - true_beta2), na.rm = TRUE)
  
  bias_TE <- mean(results_mean_TE, na.rm = TRUE) - true_TE
  rmse_TE <- sqrt(mean((results_mean_TE - true_TE)^2, na.rm = TRUE))
  mean_TE_est <- mean(results_mean_TE, na.rm = TRUE)
  
  # Output results
  cat("\n========== Results for mu =", true_mu, "==========\n")
  cat("Successful replications:", n_success, "of", R, "\n\n")
  cat("--- Bias ---\n")
  cat(sprintf("beta0: %.4f\n", bias_beta0))
  cat(sprintf("beta1: %.4f\n", bias_beta1))
  cat(sprintf("beta2: %.4f\n", bias_beta2))
  
  cat("\n--- RMSE ---\n")
  cat(sprintf("beta0: %.4f\n", rmse_beta0))
  cat(sprintf("beta1: %.4f\n", rmse_beta1))
  cat(sprintf("beta2: %.4f\n", rmse_beta2))
  
  cat("\n--- MAE ---\n")
  cat(sprintf("beta0: %.4f\n", mae_beta0))
  cat(sprintf("beta1: %.4f\n", mae_beta1))
  cat(sprintf("beta2: %.4f\n", mae_beta2))
  
  cat("\n--- Technical Efficiency ---\n")
  cat(sprintf("Estimated Mean TE: %.3f\n", mean_TE_est))
  cat(sprintf("Bias TE: %.4f\n", bias_TE))
  cat(sprintf("RMSE(TE): %.4f\n", rmse_TE))
  
  # Return results
  results_df <- data.frame(
    DGP = c(paste0("NK-Low (mu=", true_mu, ")"), "", "", "", "", "", ""),
    Parameter = c("beta0", "beta1", "beta2", "sigma_v", "sigma_u", "mu", "Mean_TE"),
    Bias = c(bias_beta0, bias_beta1, bias_beta2, NA, NA, NA, bias_TE),
    RMSE = c(rmse_beta0, rmse_beta1, rmse_beta2, NA, NA, NA, rmse_TE),
    MAE = c(mae_beta0, mae_beta1, mae_beta2, NA, NA, NA, NA)
  )
  
  return(results_df)
}

# ----------------------------------------------------------------------------
# Run simulations for extreme cases (R = 1000 as in manuscript)
# ----------------------------------------------------------------------------

cat("\n========================================\n")
cat("Table 5: Extreme Nakagami Cases\n")
cat("========================================\n")

if (!dir.exists("outputs")) dir.create("outputs")

cat("\nRunning simulation for NK-Low (mu=0.6, sigma_u=0.5) with R = 1000...\n")
results_low <- run_simulation_extreme(n = 200, R = 1000, 
                                      true_mu = 0.6, true_sigma_u = 0.5, 
                                      true_TE = 0.689)
if (!is.null(results_low)) write.csv(results_low, "outputs/table5_low.csv", row.names = FALSE)

cat("\nRunning simulation for NK-High (mu=5.0, sigma_u=1.0) with R = 1000...\n")
results_high <- run_simulation_extreme(n = 200, R = 1000, 
                                       true_mu = 5.0, true_sigma_u = 1.0, 
                                       true_TE = 0.387)
if (!is.null(results_high)) write.csv(results_high, "outputs/table5_high.csv", row.names = FALSE)

cat("\nTable 5 simulation completed successfully.\n")


