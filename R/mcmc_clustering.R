#' MCMC clustering for admixture dates
#' 
#' @param dates Vector of admixture dates
#' @param std_errors Standard errors for dates
#' @param init_cluster_means Initial cluster means
#' @param cluster_assignments Initial cluster assignments
#' @param num_iter Number of MCMC iterations
#' @param burn_in Number of burn-in iterations
#' @param k Number of clusters
#' @param thin Thinning interval
#' @param sample_ages Whether to sample ages
#' @param age_ranges Age range matrix
#' @param sampling_ages Initial sampling ages
#' @param lambda Lambda parameter
#' @return List with MCMC results
#' @importFrom stats quantile runif sd
#' @importFrom utils flush.console
#' @keywords internal
mcmc_clustering <- function(dates, std_errors, init_cluster_means, cluster_assignments,
                            num_iter = 200000, burn_in = 100000, k = 2, thin = 100,
                            sample_ages = FALSE, age_ranges = NULL, sampling_ages = NULL,
                            lambda = 10) {
  
  n <- length(dates)
  lambda_age <- 1
  lambda_cluster <- 10
  
  # Helper function for confidence intervals
  get_ci <- function(x) {
    q <- stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
    list(mean = mean(x, na.rm = TRUE), lower = q[1], upper = q[2])
  }
  
  # Identify fixed and free age individuals
  fixed_age_inds <- which(age_ranges[, 1] == age_ranges[, 2])
  free_age_inds <- which(age_ranges[, 1] != age_ranges[, 2])
  
  # Handle k=1 case
  if (k == 1 && (!sample_ages || length(free_age_inds) == 0)) {
    return(handle_k1_case(dates, std_errors, sampling_ages, n, lambda))
  }
  
  # Initialize MCMC variables
  acceptance_rate <- 0
  time_traveler_log <- list()
  
  # Lambda-specific likelihood tracking (KEPT AS REQUESTED)
  log_likelihoods_lambda1 <- c()
  log_likelihoods_lambda10 <- c()
  
  current_cluster_means <- init_cluster_means
  adjusted_dates <- dates + sampling_ages
  
  # Calculate initial likelihood
  log_likelihood <- nig_log_likelihood(cluster_assignments, current_cluster_means,
                                       std_errors, adjusted_dates, n, lambda)$total
  
  # Storage for post burn-in samples
  log_likelihoods_lambda1_post_burn_in <- numeric(0)
  log_likelihoods_lambda10_post_burn_in <- numeric(0)
  sampling_age_samples <- matrix(NA, nrow = (num_iter - burn_in)/thin, ncol = n)
  cluster_mean_samples <- matrix(NA, nrow = (num_iter - burn_in)/thin, ncol = k)
  co_clustering_matrix <- matrix(0, nrow = n, ncol = n)
  
  # Best solution tracking
  proposed_sampling_ages <- sampling_ages
  proposed_means <- init_cluster_means
  best_log_likelihood <- log_likelihood
  best_cluster_assignments <- cluster_assignments
  best_sampling_ages <- sampling_ages
  best_cluster_means <- current_cluster_means
  
  # MCMC Loop
  for (iter in 1:num_iter) {
    if (iter %% 1000 == 0 || iter == 1) {
      cat(sprintf("\rMCMC Progress: %d/%d (%.1f%%)     ",
                  iter, num_iter, 100 * iter / num_iter))
      utils::flush.console()
    }
    
    # Propose move
    proposed_assignments <- cluster_assignments
    proposed_sampling_ages <- sampling_ages
    move_type <- NULL
    
    # Choose move type and set appropriate lambda
    if (sample_ages && length(free_age_inds) && (k == 1 || stats::runif(1) < 0.5)) {
      move_type <- "age"
      i <- sample(free_age_inds, 1)
      proposed_age <- stats::runif(1, age_ranges[i, 1], age_ranges[i, 2])
      proposed_sampling_ages[i] <- proposed_age
      lambda_used <- lambda_age
    } else {
      move_type <- "cluster"
      i <- sample(1:n, 1)
      proposed_cluster <- sample(setdiff(1:k, cluster_assignments[i]), 1)
      proposed_assignments[i] <- proposed_cluster
      lambda_used <- lambda_cluster
    }
    
    proposed_adjusted_dates <- dates + proposed_sampling_ages
    
    # Calculate current likelihood
    current_log_likelihood <- nig_log_likelihood(cluster_assignments, current_cluster_means,
                                                 std_errors, adjusted_dates, n, lambda_used,
                                                 age_ranges = age_ranges)$total
    
    # Update proposed cluster means and check for time travelers
    for (j in 1:k) {
      if (sum(proposed_assignments == j) > 0) {
        cluster_dates <- proposed_adjusted_dates[proposed_assignments == j]
        cluster_errors <- std_errors[proposed_assignments == j]
        weights <- 1 / cluster_errors^2
        proposed_means[j] <- sum(cluster_dates * weights) / sum(weights)
        
        # Check for time travelers
        if (!is.null(age_ranges) && !is.null(sampling_ages)) {
          cluster_inds <- which(proposed_assignments == j)
          time_travelers <- which(proposed_means[j] < age_ranges[cluster_inds, 1])
          if (length(time_travelers) > 0) {
            time_traveler_inds <- cluster_inds[time_travelers]
            time_traveler_log[[length(time_traveler_log) + 1]] <- list(
              iter = iter,
              cluster = j,
              individuals = time_traveler_inds
            )
          }
        }
      }
    }
    
    # Relabel clusters by mean
    valid <- !is.na(proposed_means)
    ordered_clusters <- order(proposed_means[valid])
    relabel_map <- rep(NA, k)
    relabel_map[which(valid)] <- match(which(valid), ordered_clusters)
    proposed_assignments <- relabel_map[proposed_assignments]
    proposed_means <- proposed_means[order(proposed_means, na.last = TRUE)]
    
    # Calculate new likelihood
    new_log_likelihood <- nig_log_likelihood(proposed_assignments, proposed_means,
                                             std_errors, proposed_adjusted_dates,
                                             n, lambda_used, age_ranges = age_ranges)$total
    
    # Accept/reject step
    log_accept_ratio <- new_log_likelihood - current_log_likelihood
    if (log_accept_ratio > 0 || log(stats::runif(1)) < log_accept_ratio) {
      acceptance_rate <- acceptance_rate + 1
      cluster_assignments <- proposed_assignments
      if (sample_ages) sampling_ages <- proposed_sampling_ages
      adjusted_dates <- proposed_adjusted_dates
      current_cluster_means <- proposed_means
      log_likelihood <- new_log_likelihood
      
      # Update best solution
      if (log_likelihood > best_log_likelihood) {
        best_log_likelihood <- log_likelihood
        best_cluster_assignments <- cluster_assignments
        if (sample_ages) best_sampling_ages <- sampling_ages
        best_cluster_means <- current_cluster_means
      }
    }
    
    # Store likelihoods based on move type (KEPT AS REQUESTED)
    if (move_type == "age") {
      log_likelihoods_lambda1 <- c(log_likelihoods_lambda1, new_log_likelihood)
    } else {
      log_likelihoods_lambda10 <- c(log_likelihoods_lambda10, new_log_likelihood)
    }
    
    # Store samples after burn-in
    if (iter > burn_in && (iter - burn_in) %% thin == 0) {
      # Co-clustering matrix
      for (i_pair in 1:(n - 1)) {
        for (j_pair in (i_pair + 1):n) {
          if (cluster_assignments[i_pair] == cluster_assignments[j_pair]) {
            co_clustering_matrix[i_pair, j_pair] <- co_clustering_matrix[i_pair, j_pair] + 1
            co_clustering_matrix[j_pair, i_pair] <- co_clustering_matrix[j_pair, i_pair] + 1
          }
        }
      }
      
      # Store lambda-specific likelihoods post burn-in
      if (move_type == "age") {
        log_likelihoods_lambda1_post_burn_in <- c(log_likelihoods_lambda1_post_burn_in, log_likelihood)
      } else {
        log_likelihoods_lambda10_post_burn_in <- c(log_likelihoods_lambda10_post_burn_in, log_likelihood)
      }
      
      # Store samples
      if (sample_ages) sampling_age_samples[(iter - burn_in)/thin, ] <- sampling_ages
      cluster_mean_samples[(iter - burn_in)/thin, ] <- current_cluster_means
    }
  }
  
  # Calculate mean likelihoods by lambda (AS IN ORIGINAL)
  mean_log_likelihood_lambda1 <- mean(log_likelihoods_lambda1_post_burn_in)
  mean_log_likelihood_lambda10 <- mean(log_likelihoods_lambda10_post_burn_in)
  
  # Post-processing
  return(process_mcmc_results(
    n, k, dates, std_errors, best_cluster_assignments, best_sampling_ages, 
    best_cluster_means, sampling_age_samples, cluster_mean_samples, 
    co_clustering_matrix, acceptance_rate, num_iter, burn_in, thin, 
    sample_ages, age_ranges, time_traveler_log, mean_log_likelihood_lambda1, 
    mean_log_likelihood_lambda10
  ))
}