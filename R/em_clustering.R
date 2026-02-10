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
    cat("No ancients detected - using EM initialisation\n")

    em_result <- em_clustering_original(
      prepared_data$dates_original + prepared_data$midpoint_ages,
      std_errors = prepared_data$std_errors,
      n_init = 10000,
      k = k_current
    )
    initial_means <- em_result$cluster_means
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
        sampling_ages = runif(
          length(prepared_data$dates_original),
          prepared_data$age.matrix[, 1],
          prepared_data$age.matrix[, 2]
        ),
        sample_ages = sample_age_est,
        age_ranges = prepared_data$age.matrix
      )
    }

  } else {
    # Has ancients - need to check feasibility first
    cat("Ancients detected - checking feasibility\n")

    mc_res <- mc_feasibility_with_solutions(
      dates      = prepared_data$dates_original,
      std_errors = prepared_data$std_errors,
      k          = k_current,
      age_ranges = prepared_data$age.matrix,
      n_draws    = 200,
      verbose    = FALSE
    )

    cat(sprintf(
      "Feasible in %d / %d draws (%.1f%%)\n",
      mc_res$n_feasible, mc_res$n_draws, 100 * mc_res$prop_feasible
    ))

    cat(sprintf(
      "Draws with at least one infeasible assignment: %d / %d (%.1f%%)\n",
      mc_res$n_draws_with_infeasible_assignments,
      mc_res$n_draws,
      100 * mc_res$prop_draws_with_infeasible
    ))

    if (mc_res$n_feasible == 0) {
      cat("No feasible solutions found for k =", k_current, ". Skipping to next k.\n")
      return(NULL)
    }

    # Strong “all solutions feasible” criterion:
    # every draw has some feasible assignment, and we never saw an infeasible one
    all_solutions_look_feasible <- (
      mc_res$prop_feasible == 1.0 &&
        mc_res$prop_draws_with_infeasible == 0
    )

    if (all_solutions_look_feasible) {
      cat("All tested assignments were feasible across all draws - using EM initialisation\n")

      em_result <- em_clustering_original(
        prepared_data$dates_original + prepared_data$midpoint_ages,
        std_errors = prepared_data$std_errors,
        n_init = 10000,
        k = k_current
      )
      initial_means <- em_result$cluster_means
      cluster_assignments <- apply(em_result$responsibilities, 1, which.max)

      for (chain in 1:chains) {
        print(chain)
        mcmc_results[[chain]] <- mcmc_clustering(
          dates = prepared_data$dates_original,
          std_errors = prepared_data$std_errors,
          init_cluster_means = initial_means,
          cluster_assignments = cluster_assignments,
          k = k_current,
          sampling_ages = runif(
            length(prepared_data$dates_original),
            prepared_data$age.matrix[, 1],
            prepared_data$age.matrix[, 2]
          ),
          sample_ages = sample_age_est,
          age_ranges = prepared_data$age.matrix
        )
      }

    } else {
      cat("Feasible solutions exist, but infeasible assignments were observed - using random feasible solutions for initialisation\n")

      for (chain in 1:chains) {
        # Pick a random feasible solution
        random_idx <- sample(1:mc_res$n_feasible, 1)
        solution   <- mc_res$solutions[[random_idx]]

        # Compute initial means from that feasible solution
        initial_means <- compute_cluster_means(
          solution$assignments,
          prepared_data$dates_original + solution$sampling_ages,
          prepared_data$std_errors,
          k_current
        )

        mcmc_results[[chain]] <- mcmc_clustering(
          dates               = prepared_data$dates_original,
          std_errors          = prepared_data$std_errors,
          init_cluster_means  = initial_means,
          cluster_assignments = solution$assignments,
          k                   = k_current,
          sampling_ages       = solution$sampling_ages,
          sample_ages         = sample_age_est,
          age_ranges          = prepared_data$age.matrix
        )
      }
    }
  }

  # Pick best result (common to all paths)
  ll_vals <- sapply(
    mcmc_results,
    function(x) x$result$final_log_likelihood_best_sampling_ages[1]
  )
  best_idx <- which.max(ll_vals)
  mcmc_result_best <- mcmc_results[[best_idx]]

  if (plot) {
    plot_mcmc_results(
      mcmc_result_best,
      k = k_current,
      outfile_prefix = outfile,
      plot == TRUE
    )
  }

  return(mcmc_result_best)
}

