% ============================================================================
% FILE 2: code/00_hybrid_algorithm.R
% ============================================================================
library(randtoolbox)
library(pracma)

phi_std <- function(x) dnorm(x)
Phi_std <- function(x) pnorm(x)

# ----------------------------------------------------------------------------
# Calculation of torque using numerical integration
# ----------------------------------------------------------------------------
compute_moment_numerical <- function(tau, omega, r, alpha) {
  if (alpha > 10) {
    # For large alpha, use the asymptotic approximation
    return(compute_moment_asymptotic(tau, omega, r, alpha))
  }
  
  # Simple numerical integration
  z_upper <- max(20, abs(alpha) + 10)
  z_points <- seq(0, z_upper, length.out = 5000)
  dz <- z_points[2] - z_points[1]
  
  integrand <- (z_points^r) * dnorm(z_points - alpha)
  integral <- sum(integrand * dz)
  
  numerator <- (omega^r) * integral
  denominator <- Phi_std(alpha)
  
  if (denominator > 1e-10 && is.finite(numerator) && numerator > 0) {
    return(numerator / denominator)
  }
  return(NA)
}

# ----------------------------------------------------------------------------
# Asymptotic approximation for large alpha
# ----------------------------------------------------------------------------
compute_moment_asymptotic <- function(tau, omega, r, alpha) {
  abs_alpha <- abs(alpha)
  sign_alpha <- sign(alpha)
  
  if (alpha > 0) {
    # Large positive alpha: Q rarely becomes zero
    # The torque is approximately equal to tau^r
    return(tau^r)
  } else {
    # alpha منفی بزرگ: توزیع heavily truncated
    # Use the Halton sequence
    N_halton <- 20000
    halton_points <- halton(N_halton, dim = 1, init = TRUE)
    u_min <- Phi_std(-alpha)
    u_max <- 0.999999
    u_trans <- u_min + halton_points * (u_max - u_min)
    z_samples <- qnorm(u_trans)
    q_samples <- tau + omega * z_samples
    q_samples <- q_samples[q_samples >= 0 & is.finite(q_samples)]
    
    if (length(q_samples) > 100) {
      return(mean(q_samples^r, na.rm = TRUE))
    }
    return(NA)
  }
}

# ----------------------------------------------------------------------------
# Main function for hybrid torque calculation
# ----------------------------------------------------------------------------
compute_moment_hybrid <- function(tau, omega, r, alpha) {
  
  # Check inputs
  if (is.na(tau) || is.na(omega) || omega <= 0 || is.na(r) || r <= 0) {
    return(NA)
  }
  
  #Path 1: r is integer and small - use the recurrence relation
  if (r >= 1 && abs(r - round(r)) < 1e-6 && r <= 5) {
    m1 <- tau + omega * phi_std(alpha) / max(Phi_std(alpha), 1e-10)
    if (abs(r - 1) < 1e-6) return(m1)
    
    moments <- numeric(floor(r) + 1)
    moments[1] <- 1
    moments[2] <- m1
    for (k in 1:(floor(r)-1)) {
      moments[k+2] <- tau * moments[k+1] + k * (omega^2) * moments[k]
    }
    return(moments[floor(r)+1])
  }
  
  #Path 2: Use numerical integration
  result <- compute_moment_numerical(tau, omega, r, alpha)
  if (!is.na(result) && is.finite(result) && result > 0) {
    return(result)
  }
  
  # Path 3: Halton quasi-Monte Carlo as a last resort
  N_halton <- 30000
  halton_points <- halton(N_halton, dim = 1, init = TRUE)
  u_min <- Phi_std(-alpha)
  u_max <- 0.999999
  u_trans <- u_min + halton_points * (u_max - u_min)
  z_samples <- qnorm(pmax(pmin(u_trans, 0.999999), 0.000001))
  q_samples <- tau + omega * z_samples
  q_samples <- q_samples[q_samples >= 0 & is.finite(q_samples) & !is.na(q_samples)]
  
  if (length(q_samples) > 100) {
    result <- mean(q_samples^r, na.rm = TRUE)
    if (is.finite(result) && result > 0) {
      return(result)
    }
  }
  
  return(NA)
}

# ----------------------------------------------------------------------------
# Compute E(Q^(2mu-1))
# ----------------------------------------------------------------------------
compute_E_Q_power <- function(epsilon, sigma_v, mu, sigma_u) {
  if (sigma_v <= 0 || sigma_u <= 0 || mu <= 0) return(NA)
  
  A_val <- mu/(sigma_u^2) + 1/(2 * sigma_v^2)
  B_val <- -epsilon / (sigma_v^2)
  tau <- B_val / (2 * A_val)
  omega <- sqrt(1/(2 * A_val))
  alpha <- tau / omega
  r <- 2 * mu - 1
  
  compute_moment_hybrid(tau, omega, r, alpha)
}

# ----------------------------------------------------------------------------
# Negative log-likelihood function
# ----------------------------------------------------------------------------
neg_loglik_nakagami <- function(params, y, X) {
  n <- length(y)
  K <- ncol(X)
  beta <- params[1:K]
  sigma_v <- exp(params[K+1])
  sigma_u <- exp(params[K+2])
  mu_val <- exp(params[K+3])
  
  # Basic constraints
  if (sigma_v < 0.01 || sigma_v > 10) return(1e10)
  if (sigma_u < 0.01 || sigma_u > 10) return(1e10)
  if (mu_val < 0.2 || mu_val > 10) return(1e10)
  
  epsilon <- y - X %*% beta
  
  ll_i <- numeric(n)
  for (i in 1:n) {
    E_Q <- compute_E_Q_power(epsilon[i], sigma_v, mu_val, sigma_u)
    if (is.na(E_Q) || E_Q <= 1e-10 || !is.finite(E_Q)) return(1e10)
    
    A_val <- mu_val/(sigma_u^2) + 1/(2 * sigma_v^2)
    tau_i <- -epsilon[i]/(sigma_v^2) / (2 * A_val)
    omega_i <- sqrt(1/(2 * A_val))
    alpha_i <- tau_i / omega_i
    
    term1 <- log(2) + mu_val * log(mu_val) - lgamma(mu_val) - 
             mu_val * log(sigma_u^2) - 0.5 * log(sigma_v^2) - mu_val * log(2 * A_val)
    term2 <- -epsilon[i]^2/(2 * sigma_v^2) + epsilon[i]^2/(4 * A_val * sigma_v^4)
    term3 <- log(max(Phi_std(alpha_i), 1e-10))
    term4 <- log(max(E_Q, 1e-10))
    
    ll_i[i] <- term1 + term2 + term3 + term4
    
    if (!is.finite(ll_i[i]) || is.na(ll_i[i])) return(1e10)
  }
  
  -sum(ll_i)
}

# ----------------------------------------------------------------------------
# Calculation of technical efficiency
# ----------------------------------------------------------------------------
compute_TE <- function(epsilon, sigma_v, mu, sigma_u, beta, X_row) {
  A_val <- mu/(sigma_u^2) + 1/(2 * sigma_v^2)
  tau <- -epsilon/(sigma_v^2) / (2 * A_val)
  omega <- sqrt(1/(2 * A_val))
  alpha <- tau / omega
  
  r <- 2 * mu - 1
  E_Q_power <- compute_moment_hybrid(tau, omega, r, alpha)
  
  tau_star <- tau - omega^2
  alpha_star <- tau_star / omega
  E_Qstar_power <- compute_moment_hybrid(tau_star, omega, r, alpha_star)
  
  if (is.na(E_Q_power) || is.na(E_Qstar_power) || E_Q_power <= 0) {
    return(0.5)
  }
  
  TE <- exp(-tau + omega^2/2) * (Phi_std(alpha_star) / max(Phi_std(alpha), 1e-10)) * 
        (E_Qstar_power / E_Q_power)
  
  return(min(max(TE, 0.01), 0.99))
}

cat("Core hybrid algorithm functions loaded successfully (stable version).\n")

