# Clustering Algorithms
#' @importFrom stats dnorm var runif quantile sd kmeans
#' @importFrom utils flush.console head


# Original EM clustering function (without time traveler checks)
em_clustering_original <- function(dates, std_errors, k = 2, cluster_means=sample(dates, k), max_iter = 10000, epsilon = 1e-6, n_init = 10) {
  n <- length(dates)
  best_log_likelihood <- -Inf
  best_cluster_means <- NULL
  best_responsibilities <- NULL
  for (init in 1:n_init) {
    # Random initialization of cluster means

    log_likelihood <- -Inf
    diff <- Inf
    iter <- 0

    while (!is.na(diff) && diff > epsilon && iter < max_iter) {
      iter <- iter + 1

      # E-step: Calculate responsibilities
      responsibilities <- matrix(0, n, k)
      for (i in 1:n) {
        for (j in 1:k) {
          responsibilities[i, j] <- dnorm(dates[i], mean = cluster_means[j], sd = std_errors[i])
        }
        total_responsibilities <- sum(responsibilities[i, ])
        if (total_responsibilities == 0) total_responsibilities <- 1e-10
        responsibilities[i, ] <- responsibilities[i, ] / total_responsibilities
      }

      # M-step: Update cluster means
      new_cluster_means <- numeric(k)
      for (j in 1:k) {
        total_responsibilities_j <- sum(responsibilities[, j])
        if (total_responsibilities_j == 0) total_responsibilities_j <- 1e-10  # Avoid division by zero
        new_cluster_means[j] <- sum(responsibilities[, j] * dates) / total_responsibilities_j
      }

      # Check for convergence with log-likelihood calculation
      new_log_likelihood <- sum(log(rowSums(responsibilities + 1e-10)))
      diff <- abs(new_log_likelihood - log_likelihood)
      if (is.nan(diff)) break
      log_likelihood <- new_log_likelihood
      cluster_means <- new_cluster_means
    }

    if (log_likelihood > best_log_likelihood) {
      best_log_likelihood <- log_likelihood
      best_cluster_means <- cluster_means
      best_responsibilities <- responsibilities
    }
  }
  list(
    cluster_means = best_cluster_means,
    responsibilities = best_responsibilities
  )
}

# Two-stage EM clustering with CORRECT time traveler logic
em_clustering <- function(dates, std_errors, k = 2, cluster_means=sample(dates, k),
                                    max_iter = 10000, epsilon = 1e-6, n_init = 10,
                                    age_ranges = NULL, sampling_ages = NULL) {

  cat("=== Stage 1: Initial EM clustering without constraints ===\n")


  stage1_result <- em_clustering_original(dates, std_errors, k, cluster_means, max_iter, epsilon, n_init)

  # Check for time travelers in Stage 1 results
  if (!is.null(age_ranges) && !is.null(sampling_ages)) {
    cat("\n=== Checking for time travelers in initial clustering ===\n")

    cluster_assignments <- apply(stage1_result$responsibilities, 1, which.max)
    time_traveler_found <- FALSE
    all_time_travelers <- c()

    for (j in 1:k) {
      cluster_inds <- which(cluster_assignments == j)
      if (length(cluster_inds) > 0) {
        cluster_mean <- stage1_result$cluster_means[j]
        # CORRECT: Time travelers are when joint date > sampling age upper bound
        # (joint date further back than when individual lived)
        time_travelers <- which(cluster_mean < age_ranges[cluster_inds, 2])

        if (length(time_travelers) > 0) {
          time_traveler_inds <- cluster_inds[time_travelers]
          all_time_travelers <- c(all_time_travelers, time_traveler_inds)
         # cat("Cluster", j, "mean:", round(cluster_mean, 2), "years ago")
        #  cat(" - Time travelers:", paste(time_traveler_inds, collapse = ", "))
         # cat(" (sampling ages:", paste(round(age_ranges[time_traveler_inds, 2], 2), collapse = ", "), ")\n")
          time_traveler_found <- TRUE
        }
      }
    }

    if (time_traveler_found) {
      unique_time_travelers <- unique(all_time_travelers)
      cat("\nWARNING: Time travelers detected! Individual(s)", paste(unique_time_travelers, collapse = ", "),
          "would have joint admixture date more recent than their sampling age.\n")
    #  cat("We suggest removing individual(s)", paste(unique_time_travelers, collapse = ", "),
      #   "OR increasing k to", k+1, "\n")

      cat("\n=== Stage 2: Re-running EM clustering with time traveler constraints ===\n")

      # Stage 2: Run EM clustering with time traveler constraints
      stage2_result <- em_clustering_safe(dates, std_errors, k, cluster_means, max_iter, epsilon, n_init, age_ranges, sampling_ages)

      cat("EM clustering completed with time traveler constraints.\n")
      return(list(
        cluster_means = stage2_result$cluster_means,
        responsibilities = stage2_result$responsibilities,
        time_travelers_detected = unique_time_travelers,
        used_constraints = TRUE
      ))

    } else {
      cat("No time travelers detected. Using initial clustering results.\n")
      return(list(
        cluster_means = stage1_result$cluster_means,
        responsibilities = stage1_result$responsibilities,
        time_travelers_detected = NULL,
        used_constraints = FALSE
      ))
    }
  } else {
    # No age constraints provided, return Stage 1 results
    return(list(
      cluster_means = stage1_result$cluster_means,
      responsibilities = stage1_result$responsibilities,
      time_travelers_detected = NULL,
      used_constraints = FALSE
    ))
  }
}

# Also fix the safe EM clustering function
em_clustering_safe <- function(dates, std_errors, k = 2, cluster_means=sample(dates, k), max_iter = 10000, epsilon = 1e-6, n_init = 10, age_ranges = NULL, sampling_ages = NULL) {
  n <- length(dates)
  best_log_likelihood <- -Inf
  best_cluster_means <- NULL
  best_responsibilities <- NULL

  for (init in 1:n_init) {
    log_likelihood <- -Inf
    diff <- Inf
    iter <- 0

    while (!is.na(diff) && diff > epsilon && iter < max_iter) {
      iter <- iter + 1

      # E-step: Calculate responsibilities
      responsibilities <- matrix(0, n, k)
      for (i in 1:n) {
        for (j in 1:k) {
          responsibilities[i, j] <- dnorm(dates[i], mean = cluster_means[j], sd = std_errors[i])
          # ADDED: Set responsibility to 0 if this would be a time traveler assignment
          if (!is.null(age_ranges) && !is.null(sampling_ages)) {
            if (age_ranges[i, 2] > 0 && cluster_means[j] < age_ranges[i, 2]) {
              responsibilities[i, j] <- 0
            }
          }

        }
        total_responsibilities <- sum(responsibilities[i, ])
        if (total_responsibilities == 0) total_responsibilities <- 1e-10
        responsibilities[i, ] <- responsibilities[i, ] / total_responsibilities
      }

      # M-step: Update cluster means with time traveler check
      new_cluster_means <- numeric(k)
      for (j in 1:k) {
        total_responsibilities_j <- sum(responsibilities[, j])
        if (total_responsibilities_j == 0) total_responsibilities_j <- 1e-10


        proposed_mean <- sum(responsibilities[, j] * dates) / total_responsibilities_j

        # Check for time travelers if age_ranges provided
        if (!is.null(age_ranges) && !is.null(sampling_ages)) {
          current_assignments <- apply(responsibilities, 1, which.max)
          cluster_inds <- which(current_assignments == j)

          if (length(cluster_inds) > 0) {
            # CORRECT: Time travelers are when joint date > sampling age upper bound
            time_travelers <- which(proposed_mean < age_ranges[cluster_inds, 2])

            if (length(time_travelers) > 0) {
              time_traveler_inds <- cluster_inds[time_travelers]
              # Set responsibilities to zero for time travelers in this cluster
              responsibilities[time_traveler_inds, j] <- 0
              # Renormalize responsibilities for affected individuals
              for (tt in time_traveler_inds) {
                total_resp <- sum(responsibilities[tt, ])
                if (total_resp > 0) {
                  responsibilities[tt, ] <- responsibilities[tt, ] / total_resp
                } else {
                  responsibilities[tt, sample(1:k, 1)] <- 1
                }
              }
              # Recompute mean without time travelers
              total_responsibilities_j <- sum(responsibilities[, j])
              if (total_responsibilities_j > 1e-10) {
                new_cluster_means[j] <- sum(responsibilities[, j] * dates) / total_responsibilities_j
              } else {
                new_cluster_means[j] <- mean(dates)
              }
            } else {
              new_cluster_means[j] <- proposed_mean
            }
          } else {
            new_cluster_means[j] <- proposed_mean
          }
        } else {
          new_cluster_means[j] <- proposed_mean
        }
      }

      # Check for convergence
      new_log_likelihood <- sum(log(rowSums(responsibilities + 1e-10)))
      diff <- abs(new_log_likelihood - log_likelihood)
      if (is.nan(diff)) break
      log_likelihood <- new_log_likelihood
      cluster_means <- new_cluster_means
    }

    if (log_likelihood > best_log_likelihood) {
      best_log_likelihood <- log_likelihood
      best_cluster_means <- cluster_means
      best_responsibilities <- responsibilities
    }
  }

  list(
    cluster_means = best_cluster_means,
    responsibilities = best_responsibilities
  )
}



mcmc_clustering <- function(dates, std_errors, init_cluster_means, cluster_assignments,
                            num_iter = 200000, burn_in = 100000, k = 2, thin=100,
                            sample_ages = FALSE, age_ranges = NULL, sampling_ages = NULL,
                            lambda = 10) {
  n <- length(dates)
  lambda_age <- 1
  lambda_cluster <- 10
  log_likelihoods_lambda1 <- c()
  log_likelihoods_lambda10 <- c()
  get_ci <- function(x) {
    q <- quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
    list(mean = mean(x, na.rm = TRUE), lower = q[1], upper = q[2])
  }
  fixed_age_inds <- which(age_ranges[, 1] == age_ranges[, 2])
  free_age_inds <- which(age_ranges[, 1] != age_ranges[, 2])
  if (k == 1 && (!sample_ages || length(free_age_inds) == 0)) {
    dates <- dates + sampling_ages
    weights <- 1 / std_errors^2
    weighted_mean <- sum(dates * weights) / sum(weights)
    se_weighted <- sqrt(1 / sum(weights))
    ci <- weighted_mean + c(-1.96, 1.96) * se_weighted

    final_ll <- nig_log_likelihood(cluster_assignments = rep(1, n),
                                   cluster_means = weighted_mean,
                                   std_errors = std_errors,
                                   dates = dates,
                                   n = n, lambda)$total
    final_ll_perind <- nig_log_likelihood(cluster_assignments = rep(1, n),
                                          cluster_means = weighted_mean,
                                          std_errors = std_errors,
                                          dates = dates,
                                          n = n, lambda)$per_individual

    result_df <- data.frame(
      individual = 1:n,
      cluster = rep(1, n),

      sampling_age_mean = sampling_ages,
      sampling_age_best = sampling_ages,
      sampling_age_mean_lower = rep(NA, n),
      sampling_age_mean_upper = rep(NA, n),

      # Joint dates based on mean sampling ages
      joint_date_mean = rep(weighted_mean, n),
      joint_date_lower_mean = rep(ci[1], n),
      joint_date_upper_mean = rep(ci[2], n),

      # Joint dates based on best sampling ages
      joint_date_best = rep(weighted_mean, n),
      joint_date_lower_best = rep(ci[1], n),
      joint_date_upper_best = rep(ci[2], n),

      # Uncertainties
      weighted_joint_date_lower = rep(ci[1], n),
      weighted_joint_date_upper = rep(ci[2], n),
      combined_joint_date_lower_mean = rep(ci[1], n),
      combined_joint_date_upper_mean = rep(ci[2], n),
      combined_joint_date_lower_best = rep(ci[1], n),
      combined_joint_date_upper_best = rep(ci[2], n),

      joint_date_sd_weighted = rep(se_weighted, n),
      sample_ages_sd_mcmc = rep(se_weighted, n),
      joint_date_sd_combined = rep(se_weighted, n),

      # Likelihoods
      mean_log_likelihood = rep(final_ll, n),
      final_log_likelihood_mean_sampling_ages = rep(final_ll, n),
      final_log_likelihood_mean_sampling_ages_perind = final_ll_perind,
      final_log_likelihood_best_sampling_ages = rep(final_ll, n),
      final_log_likelihood_best_sampling_ages_perind = final_ll_perind,

      # MCMC
      acceptance_rate = rep(NA, n)
    )

    co_clustering_df <- NULL

    return(list(
      result = result_df,
      co_clustering_matrix = co_clustering_df,
      cluster_mean_samples = NULL,
      sampling_age_samples = NULL
    ))
  }
  acceptance_rate <- 0
  time_traveler_log <- list()

#  cluster_assignments <- apply(init_responsibilities, 1, which.max)
  proposed_assignments <- cluster_assignments
  current_cluster_means <- init_cluster_means
  adjusted_dates <- dates + sampling_ages
  log_likelihood <- nig_log_likelihood(cluster_assignments, current_cluster_means,
                                       std_errors, adjusted_dates, n, lambda)$total
  log_likelihoods_lambda1_post_burn_in <- numeric(0)
  log_likelihoods_lambda10_post_burn_in <- numeric(0)
  log_likelihoods_post_burn_in <- numeric(0)
  sampling_age_samples <- matrix(NA, nrow = (num_iter - burn_in)/thin, ncol = n)
  cluster_mean_samples <- matrix(NA, nrow = (num_iter - burn_in)/thin, ncol = k)
  co_clustering_matrix <- matrix(0, nrow = n, ncol = n)
  proposed_sampling_ages <- sampling_ages
  proposed_means <- init_cluster_means
  best_log_likelihood <- log_likelihood
  best_cluster_assignments <- cluster_assignments
  best_sampling_ages <- sampling_ages
  best_cluster_means <- current_cluster_means
  for (iter in 1:num_iter) {
    if (iter %% 1000 == 0 || iter == 1) {
      cat(sprintf("\rMCMC Progress: %d/%d (%.1f%%)     ",
                  iter, num_iter, 100 * iter / num_iter))
      flush.console() }

    proposed_assignments <- cluster_assignments
    proposed_sampling_ages <- sampling_ages
    move_type <- NULL

    if (sample_ages && length(free_age_inds) && (k == 1 || runif(1) < 0.5)) {
      move_type <- "age"
      i <- sample(free_age_inds, 1)
      proposed_age <- runif(1, age_ranges[i, 1], age_ranges[i, 2])
      proposed_sampling_ages[i] <- proposed_age
      proposed_adjusted_dates <- dates + proposed_sampling_ages
      lambda_used <- lambda_age
    } else {
      move_type <- "cluster"
      i <- sample(1:n, 1)
      proposed_cluster <- sample(setdiff(1:k, cluster_assignments[i]), 1)
      proposed_assignments[i] <- proposed_cluster
      proposed_adjusted_dates <- dates + proposed_sampling_ages
      lambda_used <- lambda_cluster
    }

    current_log_likelihood <- nig_log_likelihood(cluster_assignments, current_cluster_means,
                                                 std_errors, adjusted_dates, n, lambda_used, age_ranges = age.matrix)$total

    for (j in 1:k) {
      if (sum(proposed_assignments == j) > 0) {
        cluster_dates <- proposed_adjusted_dates[proposed_assignments == j]
        cluster_errors <- std_errors[proposed_assignments == j]
        weights <- 1 / cluster_errors^2
        proposed_means[j] <- sum(cluster_dates * weights) / sum(weights)

        if (!is.null(age_ranges) && !is.null(sampling_ages)) {
          cluster_inds <- which(proposed_assignments == j)
          time_travelers <- which(proposed_means[j] < age_ranges[cluster_inds, 1])
          if (length(time_travelers) > 0) {
            time_traveler_inds <- cluster_inds[time_travelers]

            # Store instead of warning
            time_traveler_log[[length(time_traveler_log) + 1]] <- list(
              iter = iter,
              cluster = j,
              individuals = time_traveler_inds
            )
          }
        }

        #proposed_means[j] <- mean(cluster_dates)
      }
    }
    valid <- !is.na(proposed_means)
    ordered_clusters <- order(proposed_means[valid])
    relabel_map <- rep(NA, k)
    relabel_map[which(valid)] <- match(which(valid), ordered_clusters)

    proposed_assignments <- relabel_map[proposed_assignments]
    proposed_means <- proposed_means[order(proposed_means, na.last = TRUE)]

    new_log_likelihood <- nig_log_likelihood(proposed_assignments, proposed_means,
                                             std_errors, proposed_adjusted_dates,
                                             n, lambda_used, age_ranges = age.matrix)$total
    #print(new_log_likelihood)
    log_accept_ratio <- new_log_likelihood - current_log_likelihood
    #print(log_accept_ratio)
    if (log_accept_ratio > 0 || log(runif(1)) < log_accept_ratio) {
      acceptance_rate <- acceptance_rate + 1
      cluster_assignments <- proposed_assignments
      if (sample_ages) sampling_ages <- proposed_sampling_ages
      adjusted_dates <- proposed_adjusted_dates
      current_cluster_means <- proposed_means
      log_likelihood <- new_log_likelihood
      if (log_likelihood > best_log_likelihood) {
        best_log_likelihood <- log_likelihood
        best_cluster_assignments <- cluster_assignments
        if (sample_ages) best_sampling_ages <- sampling_ages
        best_cluster_means <- current_cluster_means
      }
    }
    # Store likelihoods based on move type
    if (move_type == "age") {
      log_likelihoods_lambda1 <- c(log_likelihoods_lambda1, new_log_likelihood)
    } else {
      log_likelihoods_lambda10 <- c(log_likelihoods_lambda10, new_log_likelihood)
    }

    if (iter > burn_in && (iter - burn_in) %% thin == 0) {
      for (i_pair in 1:(n - 1)) {
        for (j_pair in (i_pair + 1):n) {
          if (cluster_assignments[i_pair] == cluster_assignments[j_pair]) {
            co_clustering_matrix[i_pair, j_pair] <- co_clustering_matrix[i_pair, j_pair] + 1
            co_clustering_matrix[j_pair, i_pair] <- co_clustering_matrix[j_pair, i_pair] + 1
          }
        }
      }
      if (move_type == "age") {
        log_likelihoods_lambda1_post_burn_in <- c(log_likelihoods_lambda1_post_burn_in, log_likelihood)
      } else {
        log_likelihoods_lambda10_post_burn_in <- c(log_likelihoods_lambda10_post_burn_in, log_likelihood)
      }
      if (sample_ages) sampling_age_samples[(iter - burn_in)/thin, ] <- sampling_ages
      cluster_mean_samples[(iter - burn_in)/thin, ] <- current_cluster_means
    }
  }
  # --- After the MCMC loop finishes ---
  mean_log_likelihood_lambda1 <- mean(log_likelihoods_lambda1_post_burn_in)
  mean_log_likelihood_lambda10 <- mean(log_likelihoods_lambda10_post_burn_in)
  acceptance_rate <- acceptance_rate / num_iter
  co_clustering_matrix <- co_clustering_matrix / ((num_iter - burn_in)/thin)
  # Final majority cluster assignments = best
  majority_cluster_assignments <- best_cluster_assignments
  # Summarize sampling ages
  age_summaries <- if (sample_ages) {
    apply(sampling_age_samples, 2, get_ci)
  } else {
    lapply(sampling_ages, function(x) list(mean = x, lower = x, upper = x))
  }
  # --- 1. Compute FINAL cluster means based on mean sampling ages across MCMC ---
  final_cluster_means <- numeric(k)
  final_cluster_sds <- numeric(k)
  for (j in 1:k) {
    inds_in_cluster <- which(majority_cluster_assignments == j)
    if (length(inds_in_cluster) > 0) {
      cluster_dates <- dates[inds_in_cluster] + sapply(age_summaries[inds_in_cluster], `[[`, "mean")
      cluster_errors <- std_errors[inds_in_cluster]
      weights <- 1 / cluster_errors^2
      final_cluster_means[j] <- sum(cluster_dates * weights) / sum(weights)
      final_cluster_sds[j] <- sqrt(1 / sum(weights))

      # Check for time travelers in final results
      if (!is.null(age_ranges) && !is.null(sampling_ages)) {
        time_travelers <- which(final_cluster_means[j] < age_ranges[inds_in_cluster, 1])
        if (length(time_travelers) > 0) {
          time_traveler_inds <- inds_in_cluster[time_travelers]
          warning(paste("Time travelers detected in final results! Individual(s)", paste(time_traveler_inds, collapse = ", "),
                        "have joint admixture date more recent than sampling age in cluster", j,
                        ". We suggest removing individual(s)", paste(time_traveler_inds, collapse = ", "),
                        "OR increasing k to", k+1))
        }
      }
    } else {
      final_cluster_means[j] <- NA
      final_cluster_sds[j] <- NA
    }
  }
  # --- 2. Compute cluster means using BEST sampling ages ---
  best_cluster_means_all <- numeric(k)
  for (j in 1:k) {
    inds_in_cluster <- which(majority_cluster_assignments == j)
    if (length(inds_in_cluster) > 0) {
      best_cluster_dates <- dates[inds_in_cluster] + best_sampling_ages[inds_in_cluster]
      cluster_errors <- std_errors[inds_in_cluster]
      weights <- 1 / cluster_errors^2
      best_cluster_means_all[j] <- sum(best_cluster_dates * weights) / sum(weights)
    } else {
      best_cluster_means_all[j] <- NA
    }
  }
  # --- 3. Compute joint date stats for each individual ---
  # (A) Based on mean sampling ages
  joint_date_stats <- sapply(1:n, function(i) {
    cl <- majority_cluster_assignments[i]
    mean_joint <- final_cluster_means[cl]
    se_joint <- final_cluster_sds[cl]
    lower_joint <- mean_joint - 1.96 * se_joint
    upper_joint <- mean_joint + 1.96 * se_joint
    c(mean_joint, lower_joint, upper_joint)
  })
  # (B) Based on BEST sampling ages
  joint_date_stats_best <- sapply(1:n, function(i) {
    cl <- majority_cluster_assignments[i]
    mean_joint <- best_cluster_means_all[cl]
    se_joint <- final_cluster_sds[cl]
    lower_joint <- mean_joint - 1.96 * se_joint
    upper_joint <- mean_joint + 1.96 * se_joint
    c(mean_joint, lower_joint, upper_joint)
  })
  # --- 4. Recalculate log likelihoods based on final means ---
  ll_with_mean_samples <- nig_log_likelihood(
    cluster_assignments = majority_cluster_assignments,
    cluster_means = final_cluster_means,
    std_errors = std_errors,
    dates = dates+sapply(age_summaries, `[[`, "mean"), n=n
  )$total
  ll_with_mean_samples_perind <- nig_log_likelihood(
    cluster_assignments = majority_cluster_assignments,
    cluster_means = final_cluster_means,
    std_errors = std_errors,
    dates = dates+sapply(age_summaries, `[[`, "mean"), n=n
  )$per_individual
  ll_with_best_samples <- nig_log_likelihood(
    cluster_assignments = majority_cluster_assignments,
    cluster_means = best_cluster_means_all,
    std_errors = std_errors,
    dates = dates+best_sampling_ages , n = n
  )$total
  ll_with_best_samples_perind <- nig_log_likelihood(
    cluster_assignments = majority_cluster_assignments,
    cluster_means = best_cluster_means_all,
    std_errors = std_errors,
    dates = dates+best_sampling_ages , n = n
  )$per_individual
  # --- 5. Calculate uncertainty (SDs) ---
  if (sample_ages) {
    individual_age_sd <- apply(sampling_age_samples, 2, sd)
  } else {
    individual_age_sd <- rep(0, n)
  }
  tiny_constant <- 1e-6
  joint_date_sd_mcmc <- sapply(1:k, function(cl) {
    mean_samples_for_cluster <- cluster_mean_samples[, cl]
    sd(mean_samples_for_cluster, na.rm = TRUE)
  })
  joint_date_sd_mcmc_individual <- sapply(1:n, function(i) {
    cl <- majority_cluster_assignments[i]
    joint_date_sd_mcmc[cl]
  })
  joint_date_sd_weighted <- sapply(1:n, function(i) {
    cl <- majority_cluster_assignments[i]
    final_cluster_sds[cl]
  })
  sample_ages_sd_by_cluster <- sapply(1:k, function(cl) {
    inds_in_cluster <- which(majority_cluster_assignments == j)
    if (length(inds_in_cluster) > 0) {
      cluster_individual_age_sd <- individual_age_sd[inds_in_cluster]
      cluster_individual_age_sd[cluster_individual_age_sd == 0] <- tiny_constant
      sqrt(1 / sum(1 / cluster_individual_age_sd^2))
    } else {
      NA
    }
  })
  joint_date_sd_combined <- sapply(1:n, function(i) {
    cl <- majority_cluster_assignments[i]
    sqrt(joint_date_sd_weighted[i]^2 + sample_ages_sd_by_cluster[cl]^2)
  })
  sample_ages_sd_mcmc <- sapply(1:n, function(i) {
    cl <- majority_cluster_assignments[i]
    sample_ages_sd_by_cluster[cl]
  })
  # --- 6. Combined stats for uncertainty ---
  combined_stats <- sapply(1:n, function(i) {
    mean_joint <- joint_date_stats[1, i]
    lower_combined <- mean_joint - 1.96 * joint_date_sd_combined[i]
    upper_combined <- mean_joint + 1.96 * joint_date_sd_combined[i]
    c(lower_combined, upper_combined)
  })
  combined_stats_best <- sapply(1:n, function(i) {
    mean_joint <- joint_date_stats_best[1, i]
    lower_combined <- mean_joint - 1.96 * joint_date_sd_combined[i]
    upper_combined <- mean_joint + 1.96 * joint_date_sd_combined[i]
    c(lower_combined, upper_combined)
  })
  # --- 7. Create Final Output  ---
  # At the end of your MCMC function, before returning:
  if (length(time_traveler_log) > 0) {
    # Count how often each individual was a time traveler
    all_time_travelers <- unlist(lapply(time_traveler_log, function(x) x$individuals))
    tt_counts <- sort(table(all_time_travelers), decreasing = TRUE)

    cat("\n=== Time Traveler Summary ===\n")
    cat("Time travelers detected in", length(time_traveler_log), "MCMC iterations\n")
    cat("Most problematic individuals:\n")
    print(head(tt_counts, 5))

    # Suggest action
    most_problematic <- as.numeric(names(tt_counts)[1:min(3, length(tt_counts))])
    cat("SUGGESTION: Remove individual(s)", paste(most_problematic, collapse = ", "),
        "OR increase k to", k+1, "\n")

  }
  result_df <- data.frame(
    individual = 1:n,
    cluster = majority_cluster_assignments,

    # Sampling ages
    sampling_age_mean = sapply(age_summaries, `[[`, "mean"),
    sampling_age_best = best_sampling_ages,
    sampling_age_mean_lower = sapply(age_summaries, `[[`, "lower"),
    sampling_age_mean_upper = sapply(age_summaries, `[[`, "upper"),

    # Joint dates based on mean sampling ages
    joint_date_mean = joint_date_stats[1, ],
    joint_date_lower_mean = joint_date_stats[2, ],
    joint_date_upper_mean = joint_date_stats[3, ],

    # Joint dates based on best sampling ages
    joint_date_best = joint_date_stats_best[1, ],
    joint_date_lower_best = joint_date_stats_best[2, ],
    joint_date_upper_best = joint_date_stats_best[3, ],

    # Uncertainties
    weighted_joint_date_lower = joint_date_stats[2, ],
    weighted_joint_date_upper = joint_date_stats[3, ],
    combined_joint_date_lower_mean = combined_stats[1, ],
    combined_joint_date_upper_mean = combined_stats[2, ],
    combined_joint_date_lower_best = combined_stats_best[1, ],
    combined_joint_date_upper_best = combined_stats_best[2, ],

    joint_date_sd_weighted = joint_date_sd_weighted,
    sample_ages_sd_mcmc = sample_ages_sd_mcmc,
    joint_date_sd_combined = joint_date_sd_combined,

    # Likelihoods
    mean_log_likelihood_lambda1 = mean_log_likelihood_lambda1,
    mean_log_likelihood_lambda10 = mean_log_likelihood_lambda10,
    final_log_likelihood_mean_sampling_ages = ll_with_mean_samples,
    final_log_likelihood_mean_sampling_ages_perind = ll_with_mean_samples_perind,
    final_log_likelihood_best_sampling_ages = ll_with_best_samples,
    final_log_likelihood_best_sampling_ages_perind = ll_with_best_samples_perind,

    # MCMC
    acceptance_rate = acceptance_rate
  )
  co_clustering_df <- as.data.frame(co_clustering_matrix)
  rownames(co_clustering_df) <- colnames(co_clustering_df) <- paste0("ind", 1:n)
  # --- 8. Return everything ---
  return(list(
    result = result_df,
    co_clustering_matrix = co_clustering_df,
    cluster_mean_samples = cluster_mean_samples,
    sampling_age_samples = if (sample_ages) sampling_age_samples else NULL
  ))
}


#' Process single k value clustering
#'
#' @param k_current Current k value to process
#' @param prepared_data List from prepare_clustering_data
#' @param chains Number of MCMC chains
#' @param sample_age_est Whether to estimate sample ages
#' @param plot Whether to create plots
#' @param outfile Output file prefix
#' @return MCMC results for this k value
#' @keywords internal
process_single_k <- function(k_current, prepared_data, chains, sample_age_est, plot, outfile) {
  cat("\n>>> Processing k =", k_current, "<<<\n")

  mcmc_results <- list()

  # Check if there are any ancients (non-zero age ranges)
  has_ancients <- sum(prepared_data$age.matrix[,2]) > 0

  if (!has_ancients) {
    # No ancients - go straight to EM initialization
    cat("No ancients detected - using EM initialization\n")

    em_result <- em_clustering_original(
      prepared_data$dates_original + prepared_data$midpoint_ages,
      std_errors = prepared_data$std_errors,
      n_init = 10000,
      k = k_current
    )
    initial_means <- em_result$cluster_means
    initial_responsibilities <- em_result$responsibilities
    cluster_assignments <- apply(em_result$responsibilities, 1, which.max)

    # Run MCMC with EM initialization
    for (chain in 1:chains) {
      print(chain)
      mcmc_results[[chain]] <- mcmc_clustering(
        dates = prepared_data$dates_original,
        std_errors = prepared_data$std_errors,
        init_cluster_means = initial_means,
        cluster_assignments = cluster_assignments,
        k = k_current,
        sampling_ages = runif(length(prepared_data$dates_original),
                              prepared_data$age.matrix[, 1],
                              prepared_data$age.matrix[, 2]),
        sample_ages = sample_age_est,
        age_ranges = prepared_data$age.matrix
      )
    }

  } else {
    # Has ancients - need to check feasibility first
    cat("Ancients detected - checking feasibility\n")

    mc_res <- mc_feasibility_with_solutions(
      prepared_data$dates_original,
      prepared_data$std_errors,
      k = k_current,
      age_ranges = prepared_data$age.matrix,
      n_draws = 200,
      verbose = FALSE

    )

    cat(sprintf(
      "Feasible in %d / %d draws (%.1f%%)\n",
      mc_res$n_feasible, mc_res$n_draws, 100 * mc_res$prop_feasible
    ))

    if (mc_res$n_feasible == 0) {
      cat("No feasible solutions found for k =", k_current, ". Skipping to next k.\n")
      return(NULL)
    }

    if (mc_res$prop_feasible == 1.0) {
      # 100% feasible - use EM initialization
      cat("100% feasible - using EM initialization\n")

      em_result <- em_clustering_original(
        prepared_data$dates_original + prepared_data$midpoint_ages,
        std_errors = prepared_data$std_errors,
        n_init = 10000,
        k = k_current
      )
      initial_means <- em_result$cluster_means
      initial_responsibilities <- em_result$responsibilities
      cluster_assignments <- apply(em_result$responsibilities, 1, which.max)

      # Run MCMC with EM initialization
      for (chain in 1:chains) {
        print(chain)
        mcmc_results[[chain]] <- mcmc_clustering(
          dates = prepared_data$dates_original,
          std_errors = prepared_data$std_errors,
          init_cluster_means = initial_means,
          cluster_assignments = cluster_assignments,
          k = k_current,
          sampling_ages = runif(length(prepared_data$dates_original),
                                prepared_data$age.matrix[, 1],
                                prepared_data$age.matrix[, 2]),
          sample_ages = sample_age_est,
          age_ranges = prepared_data$age.matrix
        )
      }

    } else {
      # Not 100% feasible - use random feasible solutions
      cat("Not 100% feasible - using random feasible solutions\n")

      for (chain in 1:chains) {
        # Pick a random feasible solution
        random_idx <- sample(1:mc_res$n_feasible, 1)
        solution <- mc_res$solutions[[random_idx]]

        # Convert assignments to initial means and responsibilities
        initial_means <- compute_cluster_means(
          solution$assignments,
          prepared_data$dates_original + solution$sampling_ages,
          prepared_data$std_errors,
          k_current
        )
        mcmc_results[[chain]] <- mcmc_clustering(
          dates = prepared_data$dates_original,
          std_errors = prepared_data$std_errors,
          init_cluster_means = initial_means,
          cluster_assignments = solution$assignments,
          k = k_current,
          sampling_ages = solution$sampling_ages,
          sample_ages = sample_age_est,
          age_ranges = prepared_data$age.matrix
        )
      }
    }
  }

  # Pick best result (common to all paths)
  ll_vals <- sapply(mcmc_results, function(x) x$result$final_log_likelihood_best_sampling_ages[1])
  best_idx <- which.max(ll_vals)
  mcmc_result_best <- mcmc_results[[best_idx]]

  if (plot) {
    plot_mcmc_results(mcmc_result_best, k = k_current, outfile_prefix = outfile, plot == TRUE)
  }

  return(mcmc_result_best)
}


