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
process_single_k <- function(k_current, prepared_data, chains,
                             sample_age_est, plot, outfile) {
  cat("\n>>> Processing k =", k_current, "<<<\n")
  n_inds <- length(prepared_data$dates_original)

  if (n_inds < k_current) {
    cat(sprintf(
      "Skipping k=%d: only %d individual(s) available (need n >= k).\n",
      k_current, n_inds
    ))
    return(NULL)
  }

  # 1. Always run EM first.
  em_result <- em_clustering_original(
    prepared_data$dates_original + prepared_data$midpoint_ages,
    std_errors = prepared_data$std_errors,
    n_init = 10000,
    k = k_current
  )

  em_assignments <- apply(
    em_result$responsibilities, 1, which.max
  )

  # Filled only if a chain's EM/greedy start fails validation.
  feasible_starts <- NULL
  mcmc_results <- vector("list", chains)

  for (chain in seq_len(chains)) {
    cluster_assignments <- em_assignments

    sampling_ages <- stats::runif(
      n_inds,
      prepared_data$age.matrix[, 1],
      prepared_data$age.matrix[, 2]
    )

    # 2. Apply the original greedy repair if any cluster is empty.
    if (k_current > 1) {
      init_counts <- tabulate(
        cluster_assignments, nbins = k_current
      )

      if (any(init_counts == 0)) {
        cat(sprintf(
          "Chain %d: empty EM cluster(s); reassigning greedily.\n",
          chain
        ))

        ord <- order(prepared_data$dates_original)
        cluster_assignments <- rep(NA, n_inds)
        cluster_assignments[ord[1:k_current]] <- 1:k_current

        remaining <- setdiff(
          seq_len(n_inds), ord[1:k_current]
        )

        if (length(remaining) > 0) {
          cluster_assignments[remaining] <- sample(
            1:k_current,
            length(remaining),
            replace = TRUE
          )
        }
      }
    }

    # Recompute means after any repair, using this chain's actual ages.
    initial_means <- compute_cluster_means(
      cluster_assignments,
      prepared_data$dates_original + sampling_ages,
      prepared_data$std_errors,
      k_current
    )

    # 3. Validate the EM-derived or greedy-repaired starting state.
    if (!valid_mcmc_start(
      cluster_assignments,
      initial_means,
      sampling_ages,
      prepared_data$age.matrix,
      k_current
    )) {
      cat(sprintf(
        "Chain %d: EM/greedy start is invalid; using feasibility search.\n",
        chain
      ))

      # Run feasibility at most once per k, only when needed.
      if (is.null(feasible_starts)) {
        mc_res <- mc_feasibility_with_solutions(
          dates = prepared_data$dates_original,
          std_errors = prepared_data$std_errors,
          k = k_current,
          age_ranges = prepared_data$age.matrix,
          n_draws = 200,
          verbose = FALSE
        )

        # Validate the stored states, not just their feasibility flags.
        feasible_starts <- list()

        for (solution in mc_res$solutions) {
          solution_means <- compute_cluster_means(
            solution$assignments,
            prepared_data$dates_original + solution$sampling_ages,
            prepared_data$std_errors,
            k_current
          )

          if (valid_mcmc_start(
            solution$assignments,
            solution_means,
            solution$sampling_ages,
            prepared_data$age.matrix,
            k_current
          )) {
            solution$cluster_means <- solution_means

            feasible_starts[[length(feasible_starts) + 1L]] <- solution
          }
        }

        if (length(feasible_starts) == 0L) {
          cat(sprintf(
            "No valid feasible starts found for k=%d; skipping this k.\n",
            k_current
          ))
          return(NULL)
        }
      }

      # Select a checked solution, keeping its corresponding ages.
      solution <- feasible_starts[[
        sample.int(length(feasible_starts), 1L)
      ]]

      cluster_assignments <- solution$assignments
      sampling_ages <- solution$sampling_ages
      initial_means <- solution$cluster_means
    }

    # Final check of the exact state being passed to MCMC.
    if (!valid_mcmc_start(
      cluster_assignments,
      initial_means,
      sampling_ages,
      prepared_data$age.matrix,
      k_current
    )) {
      stop(
        "Invalid MCMC initial state after EM/greedy/feasibility selection."
      )
    }

    mcmc_results[[chain]] <- mcmc_clustering(
      dates = prepared_data$dates_original,
      std_errors = prepared_data$std_errors,
      init_cluster_means = initial_means,
      cluster_assignments = cluster_assignments,
      k = k_current,
      sampling_ages = sampling_ages,
      sample_ages = sample_age_est,
      age_ranges = prepared_data$age.matrix
    )
  }

  # Keep the existing best-chain selection and plotting behaviour.
  ll_vals <- vapply(
    mcmc_results,
    function(x) {
      x$result$final_log_likelihood_best_sampling_ages[1]
    },
    numeric(1)
  )

  mcmc_result_best <- mcmc_results[[which.max(ll_vals)]]

  if (plot) {
    plot_mcmc_results(
      mcmc_result_best,
      k = k_current,
      outfile_prefix = outfile,
      plot = TRUE
    )
  }

  return(mcmc_result_best)
}
