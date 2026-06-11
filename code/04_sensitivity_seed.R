% ============================================================================
% FILE 6: code/04_sensitivity_seed.R
% ============================================================================
# Seed sensitivity analysis for simulation results
# Compares results across four different random seeds
# R = 500 for each seed (sufficient for sensitivity comparison)
# ----------------------------------------------------------------------------
# Function to run simulation for a given seed
# ----------------------------------------------------------------------------
run_seed_sensitivity <- function(n = 200, R = 500, seed_value, 
                                 true_beta0 = 1.0, true_beta1 = 0.6, true_beta2 = 0.4,
                                 true_sigma_v = 0.3, true_mu = 2.0, true_sigma_u = 0.5) {
  
  set.seed(seed_value)
  
  results_beta0 <- c()
  results_beta1 <- c()
  results_beta2 <- c()
  results_mu <- c()
  
  for (r in 1:R) {
    # Use a secondary seed within the loop for reproducibility across seeds
    set.seed(seed_value * 10000 + r)
    
    x1 <- runif(n, 1, 10)
    x2 <- rnorm(n, 5, 2)
    v <- rnorm(n, 0, true_sigma_v)
    u <- sqrt(qgamma(runif(n), shape = true_mu, rate = true_mu/true_sigma_u^2))
    y <- true_beta0 + true_beta1*x1 + true_beta2*x2 + v - u
    X_mat <- cbind(1, x1, x2)
    
    init_params <- c(mean(y), 0.6, 0.4, log(0.3), log(0.5), log(1.8))
    
    opt_result <- tryCatch({
      optim(init_params, neg_loglik_nakagami, y = y, X = X_mat,
            method = "Nelder-Mead", 
            control = list(maxit = 500, trace = 0, reltol = 1e-8))
    }, error = function(e) NULL)
    
    if (!is.null(opt_result) && opt_result$convergence == 0) {
      params_opt <- opt_result$par
      K <- ncol(X_mat)
      results_beta0 <- c(results_beta0, params_opt[1])
      results_beta1 <- c(results_beta1, params_opt[2])
      results_beta2 <- c(results_beta2, params_opt[3])
      results_mu <- c(results_mu, exp(params_opt[K+3]))
    }
    
    if (r %% 100 == 0) cat("  Seed", seed_value, "- replication", r, "of", R, "completed. Successes:", length(results_beta0), "\n")
  }
  
  if (length(results_beta0) > 0) {
    return(list(
      seed = seed_value,
      mean_beta0 = mean(results_beta0, na.rm = TRUE),
      mean_beta1 = mean(results_beta1, na.rm = TRUE),
      mean_beta2 = mean(results_beta2, na.rm = TRUE),
      mean_mu = mean(results_mu, na.rm = TRUE),
      bias_beta0 = mean(results_beta0, na.rm = TRUE) - true_beta0,
      bias_beta1 = mean(results_beta1, na.rm = TRUE) - true_beta1,
      bias_beta2 = mean(results_beta2, na.rm = TRUE) - true_beta2,
      bias_mu = mean(results_mu, na.rm = TRUE) - true_mu,
      n_success = length(results_beta0)
    ))
  } else {
    return(list(seed = seed_value, mean_beta0 = NA, bias_beta0 = NA, n_success = 0))
  }
}

# ----------------------------------------------------------------------------
# Run sensitivity analysis with four different seeds (R = 500 each)
# ----------------------------------------------------------------------------

cat("\n========================================\n")
cat("Seed Sensitivity Analysis\n")
cat("Each simulation: n=200, R=500\n")
cat("========================================\n")

seeds <- c(12345, 54321, 99999, 77777)
results_list <- list()

for (s in seeds) {
  cat("\nRunning simulation with seed =", s, "...\n")
  res <- run_seed_sensitivity(n = 200, R = 500, seed_value = s)
  results_list <- c(results_list, list(res))
  if (!is.na(res$mean_beta0)) {
    cat(sprintf("  Mean beta0: %.4f\n", res$mean_beta0))
    cat(sprintf("  Bias beta0: %.4f\n", res$bias_beta0))
    cat(sprintf("  Mean mu: %.4f\n", res$mean_mu))
    cat(sprintf("  Successful replications: %d of 500\n", res$n_success))
  } else {
    cat("  WARNING: No successful replications for this seed!\n")
  }
}

# Create summary table
summary_df <- data.frame(
  Seed = sapply(results_list, function(x) x$seed),
  Mean_Beta0 = sapply(results_list, function(x) round(x$mean_beta0, 4)),
  Bias_Beta0 = sapply(results_list, function(x) round(x$bias_beta0, 4)),
  Mean_Beta1 = sapply(results_list, function(x) round(x$mean_beta1, 4)),
  Mean_Beta2 = sapply(results_list, function(x) round(x$mean_beta2, 4)),
  Mean_Mu = sapply(results_list, function(x) round(x$mean_mu, 4)),
  Success_Rate = sapply(results_list, function(x) x$n_success)
)

cat("\n========== Seed Sensitivity Summary ==========\n")
print(summary_df)

# Save results to CSV
if (!dir.exists("outputs")) dir.create("outputs")
write.csv(summary_df, "outputs/seed_sensitivity.csv", row.names = FALSE)

cat("\nSeed sensitivity analysis completed successfully.\n")
cat("Results saved to outputs/seed_sensitivity.csv\n")

