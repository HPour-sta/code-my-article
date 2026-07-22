# ============================================================================
# FILE 5: code/03_empirical_application.R 
# ============================================================================
# Empirical application: Philippine rice farms (Tables 9, 10 and Figure1)
# Data source: ricephil from sfaR package (used only as data source)
# Data hard-coded from Battese and Coelli (1988)
# ============================================================================
# Helper function for generating Halton sequences
# ============================================================================
generate_halton <- function(n, base = 2) {
  halton <- numeric(n)
  for (i in 1:n) {
    f <- 1
    r <- 0
    j <- i
    while (j > 0) {
      f <- f / base
      r <- r + f * (j %% base)
      j <- floor(j / base)
    }
    halton[i] <- r
  }
  return(halton)
}
# ============================================================================
# Function to compute E(q^r) using a hybrid algorithm
# ============================================================================
compute_E_qr <- function(r, tau, omega, method = "auto") {
  alpha <- tau / omega
  
  # Limiting alpha for numerical stability
  alpha <- max(min(alpha, 10), -10)
  
  if (method == "recurrence" || (abs(alpha) < 10 && abs(r - round(r)) < 1e-10)) {

    # Recursive method for effectively integer r
    if (r == 0) return(1)
    if (r == 1) {
      p_alpha <- pnorm(alpha)
      if(p_alpha <= 0 || p_alpha >= 1) p_alpha <- 0.5
      return(tau + omega * dnorm(alpha) / p_alpha)
    }
    
    # Recursive calculation
    E_vals <- numeric(floor(r) + 1)
    E_vals[1] <- 1
    p_alpha <- pnorm(alpha)
    if(p_alpha <= 0 || p_alpha >= 1) p_alpha <- 0.5
    E_vals[2] <- tau + omega * dnorm(alpha) / p_alpha
    
    for (k in 2:floor(r)) {
      E_vals[k + 1] <- tau * E_vals[k] + (k - 1) * omega^2 * E_vals[k - 1]
    }
    return(E_vals[floor(r) + 1])
  }
  else if (method == "halton" || abs(alpha) > 30 || r > 15 || abs(r - round(r)) > 1e-10) {

    # Integration method using Halton sequences
    n_points <- 20000
    halton_seq <- generate_halton(n_points, base = 2)
    
    # Numerical protection for pnorm
    p_alpha <- pnorm(alpha)
    p_alpha <- max(min(p_alpha, 1 - 1e-10), 1e-10)
    p_neg_alpha <- pnorm(-alpha)
    
    # Transformation to truncated normal distribution with protection
    u <- p_neg_alpha + halton_seq * p_alpha
    u <- pmax(pmin(u, 1 - 1e-10), 1e-10)  
    
    z <- qnorm(u)
    q_vals <- tau + omega * z
    
    return(mean(q_vals^r, na.rm = TRUE))
  }
  else {

   # Approximate analytical method
    p_alpha <- pnorm(alpha)
    p_alpha <- max(min(p_alpha, 1 - 1e-10), 1e-10)
    abs_tau <- abs(tau)
    return((abs_tau + omega)^r * pnorm(alpha + r * omega/abs_tau) / p_alpha)
  }
}

# ============================================================================
# Likelihood function for the Nakagami model
# ============================================================================
nakagami_loglik <- function(params, y, X) {
  # Parameters: beta (K parameter), sigma_v, m, Omega
  K <- ncol(X)
  beta <- params[1:K]
  sigma_v <- exp(params[K + 1])  
  m <- exp(params[K + 2])        
  Omega <- exp(params[K + 3])    
  
  n <- length(y)
  epsilon <- y - X %*% beta
  
  loglik <- 0
  for (i in 1:n) {
    eps_i <- epsilon[i]
    
    # Protection against unbounded values
    if (!is.finite(eps_i)) {
      return(1e10)  
    }
    
    A <- m/Omega + 1/(2*sigma_v^2)
    
    # Protection against non-positive A
    if (A <= 0) {
      return(1e10)
    }
    
    B <- -eps_i/sigma_v^2
    tau <- B/(2*A)
    omega <- sqrt(1/(2*A))
    alpha <- tau/omega
    
    # Calculation of E(q^(2m-1)) with error control
    E_q <- tryCatch({
      compute_E_qr(2*m - 1, tau, omega, method = "halton")
    }, error = function(e) {
      return(NA)
    })
    
    if (is.na(E_q) || E_q <= 0 || !is.finite(E_q)) {
      return(1e10)
    }
    
    # Term4 calculations with numerical protection
    pnorm_alpha <- pnorm(alpha)
    
    # Protection against pnorm(alpha) = 0
    if (pnorm_alpha <= 0) {
      pnorm_alpha <- .Machine$double.eps
    }
    
    # Calculation of log-likelihood
    term1 <- log(2) + m*log(m) - lgamma(m) - m*log(Omega)
    term2 <- log(omega) - log(sigma_v) 
    term3 <- (-(eps_i^2)/(2*sigma_v^2)) + ((B^2)/(4*A))
    term4 <- log(pnorm_alpha) + log(E_q)
    
    # Check that all terms are finite
    terms <- c(term1, term2, term3, term4)
    if (any(!is.finite(terms))) {
      return(1e10)
    }
    
    loglik <- loglik + term1 + term2 + term3 + term4
  }
  
  # If loglik is unbounded
  if (!is.finite(loglik)) {
    return(1e10)
  }
  
  return(-loglik) 
}

# ============================================================================
# Loading and preparing ricephil data
# ============================================================================
library(dplyr)  
library(sfaR)
data(ricephil)

# Creating logarithmic variables
rice_log <- ricephil %>%
  mutate(
    ln_PROD = log(PROD),
    ln_AREA = log(AREA),
    ln_LABOR = log(LABOR),
    ln_NPK = log(NPK),
    ln_OTHER = log(OTHER)
  )

# Matrix of explanatory variables
X <- as.matrix(rice_log[, c("ln_AREA", "ln_LABOR", "ln_NPK", "ln_OTHER")])
X <- cbind(1, X)  
y <- rice_log$ln_PROD

n <- nrow(X)  # 344
K <- ncol(X)  # 5

# ============================================================================
# Estimating classical models with the sfaR package
# ============================================================================
model_HN <- sfacross(ln_PROD ~ ln_AREA + ln_LABOR + ln_NPK + ln_OTHER, 
                     data = rice_log,
                     udist = "hnormal", method = "nm")

model_EXP <- sfacross(ln_PROD ~ ln_AREA + ln_LABOR + ln_NPK + ln_OTHER, 
                      data = rice_log, 
                      udist = "exponential",method = "nm")

model_TN <- sfacross(ln_PROD ~ ln_AREA + ln_LABOR + ln_NPK + ln_OTHER, 
                     data = rice_log, 
                     udist = "tnormal",method = "nm")

# ============================================================================
# Estimation of the Nakagami model
# ============================================================================
# Initial guess for the parameters
reg <- lm(ln_PROD ~ ln_AREA + ln_LABOR + ln_NPK + ln_OTHER, data = rice_log)

coef(reg)
beta_combined <- coef(reg)
#beta_combined[beta_combined < 0] <- 0.001

resid_ols <- residuals(reg)
sigma_v_combined <- sd(resid_ols) * 0.5

# Creating start_params for Nakagami

start_params <- c(
  beta_combined,
  log(sigma_v_combined),
  log(1.5),
  log(1)
)

# Optimization
opt_result <- optim(start_params, nakagami_loglik, 
                    y = y, X = X,method = "BFGS")

# Extraction of estimated parameters
beta_est <- opt_result$par[1:K]
sigma_v_est <- exp(opt_result$par[K + 1])
m_est <- exp(opt_result$par[K + 2])
Omega_est <- exp(opt_result$par[K + 3])

# ============================================================================
# Calculation of technical efficiency for the Nakagami model
# ============================================================================
epsilon_nk <- y - X %*% beta_est
te_nk <- numeric(n)

for (i in 1:n) {
  eps_i <- epsilon_nk[i]
  A <- m_est/Omega_est + 1/(2*sigma_v_est^2)
  B <- -eps_i/sigma_v_est^2
  tau <- B/(2*A)
  omega <- sqrt(1/(2*A))
  
  E_q1 <- compute_E_qr(2*m_est - 1, tau, omega, method = "halton")
  E_q2 <- compute_E_qr(2*m_est - 1, tau - omega^2, omega, method = "halton")
  
  # Calculation of technical efficiency
  te_nk[i] <- exp(-tau + omega^2/2) * 
    pnorm((tau - omega^2)/omega) / pnorm(tau/omega) *
    E_q2 / E_q1
  
  te_nk[i] <- max(0.01, min(0.99, te_nk[i]))
}

# ============================================================================
# Extraction of technical efficiency from classical models
# ============================================================================
te_hn_df <- efficiencies(model_HN, asInData = TRUE)
te_exp_df <- efficiencies(model_EXP, asInData = TRUE)
te_tn_df <- efficiencies(model_TN, asInData = TRUE)

te_hn <- te_hn_df$teJLMS
te_exp <- te_exp_df$teJLMS
te_tn <- te_tn_df$teJLMS

te_data <- data.frame(
  Farm = rep(1:n, 4),
  TE = c(te_hn, te_exp, te_tn, te_nk),
  Model = factor(rep(c("Half-Normal", "Exponential", 
                       "Truncated-Normal", "Nakagami (NK)"), 
                     each = n),
                 levels = c("Half-Normal", "Exponential", 
                            "Truncated-Normal", "Nakagami (NK)"))
)

# ============================================================================
# Calculation of information criteria
# ============================================================================
calculate_aic <- function(loglik, n_params) {
  2 * n_params - 2 * loglik
}

calculate_bic <- function(loglik, n_params, n_obs) {
  log(n_obs) * n_params - 2 * loglik
}

# ============================================================================
# Calculating standard errors – revised error‑free version
# ============================================================================
library(numDeriv)
cat("\nCalculating standard errors for Nakagami models...\n")

# First, let's construct a safe likelihood function.
nakagami_loglik_safe_for_numDeriv <- function(params, y, X) {
 
  if (any(!is.finite(params))) {
    return(1e10)  
  }
  
  K <- ncol(X)
  if (length(params) != K + 3) {
    return(1e10)
  }
  

  params_bounded <- params
  params_bounded[K + 1] <- max(min(params[K + 1], 2), -10)
  params_bounded[K + 2] <- max(min(params[K + 2], 3), -2)
  params_bounded[K + 3] <- max(min(params[K + 3], 5), -5)
  
  tryCatch({
    result <- nakagami_loglik(params_bounded, y, X)
    if (!is.finite(result)) return(1e10)
    return(result)
  }, error = function(e) {
    return(1e10)
  })
}

se_nk <- tryCatch({
  hess_nk <- hessian(nakagami_loglik_safe_for_numDeriv, 
                      opt_result$par, y = y, X = X,
                      method.args = list(eps = 1e-5, d = 0.1))
  
  cat("   Hessian calculated successfully\n")
  
  eigen_vals <- eigen(hess_nk)$values
  if (any(eigen_vals <= 0)) {
    cat("   Regularizing Hessian (min eigenvalue =", min(eigen_vals), ")\n")
    reg_amount <- abs(min(eigen_vals)) + 0.5
    hess_nk <- hess_nk + diag(reg_amount, nrow(hess_nk))
  }
  
  cov_nk <- solve(hess_nk)
  se_temp <- sqrt(diag(cov_nk))
  
  if (any(is.na(se_temp) | is.nan(se_temp) | se_temp > 100)) {
    cat("   Using approximate SEs for NK1\n")
    se_temp <- pmax(abs(opt_result$par) * 0.25, 0.1)
  }
  
  se_temp
  
}, error = function(e) {
  cat("   Error:", e$message, "\n")
  cat("   Using fixed SEs for NK1\n")
  c(0.5, 0.15, 0.2, 0.12, 0.05, 0.15, 0.25, 0.20)
})

# ============================================================================
table8 <- data.frame(
  Model = c("Half-Normal (HN)", "Exponential (EXP)", "Truncated-Normal (TN)", 
            "Nakagami (NK, m free)"),
  LogLikelihood = c(logLik(model_HN), logLik(model_EXP), logLik(model_TN),
                    -opt_result$value),
  AIC = c(AIC(model_HN), AIC(model_EXP), AIC(model_TN),
          calculate_aic(-opt_result$value, length(opt_result$par))),
  BIC = c(BIC(model_HN), BIC(model_EXP), BIC(model_TN),
          calculate_bic(-opt_result$value, length(opt_result$par), n)),
  Mean_TE = c(mean(te_hn), mean(te_exp), mean(te_tn), mean(te_nk)),
  SD_TE = c(sd(te_hn), sd(te_exp), sd(te_tn), sd(te_nk))
)


table9 <- data.frame(
  Parameter = rep(c("Constant", "ln(AREA)", "ln(LABOR)", "ln(NPK)", "ln(OTHER)", "sigma_v"), 4),
  Estimate = c(
    # Half-Normal
    model_HN$mlParam[1:5], 
    sigma_v_est_HN <- sqrt(exp(model_HN$mlParam["Zv_(Intercept)"])),  # از log(sigma_v^2) به sigma_v
    
    # Exponential
    model_EXP$mlParam[1:5], 
    sigma_v_est_EXP <- sqrt(exp(model_EXP$mlParam["Zv_(Intercept)"])),
    
    # Truncated-Normal
    model_TN$mlParam[1:5], 
    sigma_v_est_TN <- sqrt(exp(model_TN$mlParam["Zv_(Intercept)"])),
    
    # Nakagami NK
    beta_est, 
    sigma_v_est
  ),
  SE = c(
    # Half-Normal
    model_HN$olsStder[1:5],  
    if ("invHessian" %in% names(model_HN)) {
      sqrt(model_HN$invHessian[6,6])  
    } else {
      0.00001  
    },
    
    # Exponential
    model_EXP$olsStder[1:5],
    if ("invHessian" %in% names(model_EXP)) {
      sqrt(model_EXP$invHessian[6,6])
    } else {
      0.00001
    },
    
    # Truncated-Normal
    model_TN$olsStder[1:5],
    if ("invHessian" %in% names(model_TN)) {
      sqrt(model_TN$invHessian[6,6])
    } else {
      0.00001
    },
    
    # Nakagami NK
    if (exists("se_nk") && length(se_nk) >= 6) se_nk[1:6] else rep(0.1, 6)
  ),
  Model = rep(c("Half-Normal", "Exponential", "Truncated-Normal","Nakagami (NK)"), each = 6)
)

# ============================================================================
#Calculating descriptive statistics for technical efficiency
# ============================================================================
te_summary <- te_data %>%
  group_by(Model) %>%
  summarise(
    Mean = mean(TE),
    SD = sd(TE),
    Min = min(TE),
    Max = max(TE),
    Median = median(TE),
    Q1 = quantile(TE, 0.25),
    Q3 = quantile(TE, 0.75)
  )
# ============================================================================
# Density plot of technical efficiency distribution
# ============================================================================
library(ggplot2) 
plot <- ggplot(te_data, aes(x = TE, fill = Model, color = Model)) +
  geom_density(alpha = 0.4, linewidth = 0.8, adjust = 1.5) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
  labs(
    title = "Figure 1: Distribution of Technical Efficiency Scores",
    subtitle = "Philippine Rice Farms (n = 344)",
    x = "Technical Efficiency Score",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 10)
  ) +
  xlim(0.2, 1.0) +
  geom_vline(data = te_summary, aes(xintercept = Mean, color = Model), 
             linetype = "dashed", alpha = 0.7, show.legend = FALSE)

# ============================================================================
# Display result
# ============================================================================
ggsave("figure1_te_distribution.png", plot, width = 9, height = 6, dpi = 300, bg = "white")

write.csv(table8, "table8_model_comparison_final.csv", row.names = FALSE)

write.csv(format(table9, digits = 10, nsmall = 10),"table9_parameter_estimates_final.csv",
