# Likelihood and Distribution Functions

estimate_alpha_beta <- function(sds) {
  var_sqs <- sds^2
  
  if (length(var_sqs) < 2 || any(is.na(var_sqs))) {
    mu_s <- mean(var_sqs, na.rm = TRUE)
    alpha <- 10
    beta <- (alpha - 1) * mu_s
  } else {
    mu_s <- mean(var_sqs)
    sigma_s2 <- var(var_sqs)
    
    if (is.na(mu_s) || is.na(sigma_s2) || mu_s <= 0 || sigma_s2 <= 0) {
      alpha <- 10
      beta <- (alpha - 1) * mu_s
    } else {
      alpha <- 2 + (mu_s^2) / sigma_s2
      #alpha<-3
      beta <- (alpha - 1) * mu_s
    }
  }
  
  return(list(alpha = alpha, beta = beta))
}


dnig <- function(x, sigma_sq, mu, lambda, alpha, beta, log = FALSE) {
  log_density <- 0.5 * log(lambda) - 0.5 * log(2 * pi * sigma_sq) + 
    alpha * log(beta) - lgamma(alpha) - 
    (alpha + 1) * log(sigma_sq) - 
    (2 * beta + lambda * (x - mu)^2) / (2 * sigma_sq)
  
  if (log) return(log_density) else return(exp(log_density))
}


nig_log_likelihood <- function(cluster_assignments, cluster_means, std_errors, dates, n, lambda = 10, age_ranges = NULL) {
  K <- length(cluster_means)
  
  
  if (!is.null(age_ranges) && K > 1) {
    for (i in 1:n) {
      cluster <- cluster_assignments[i]
      joint_date <- cluster_means[cluster]
      
      # Check if joint admixture date is more recent than lower bound of sampling age- for time travelers and give them minus infinity likelihood
      if (!is.na(joint_date) && joint_date < age_ranges[i, 1]) {
        return(list(
          total = -Inf,
          per_individual = rep(-Inf, n)
        ))
      }
    }
  }
  
  cluster_sds_list <- split(std_errors, cluster_assignments)
  cluster_params <- lapply(cluster_sds_list, estimate_alpha_beta)
  per_ind_ll <- numeric(n)
  
  for (i in 1:n) {
    cluster <- cluster_assignments[i]
    mu <- cluster_means[cluster]
    sigma_sq <- std_errors[i]^2
    
    alpha <- cluster_params[[as.character(cluster)]]$alpha
    beta  <- cluster_params[[as.character(cluster)]]$beta
    
    per_ind_ll[i] <- dnig(dates[i], sigma_sq, mu, lambda, alpha, beta, log = TRUE)
  }
  
  total_ll <- sum(per_ind_ll)
  return(list(
    total = total_ll,
    per_individual = per_ind_ll
  ))
}


current_log_likelihood <- function(cluster_assignments, cluster_means,std_errors,dates, n) {
  log_likelihood <- 0
  for (i in 1:n) {
    log_likelihood <- log_likelihood + dnorm(dates[i], mean = cluster_means[cluster_assignments[i]], 
                                             sd = std_errors[i], log = TRUE)
  }
  return(log_likelihood)
}