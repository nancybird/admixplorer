# Clustering Algorithms
#' @importFrom stats dnorm var runif quantile sd kmeans
#' @importFrom utils flush.console head


# Original EM clustering function (without time traveler checks)
em_clustering_original <- function(dates, std_errors, k = 2, cluster_means=NULL, max_iter = 10000, epsilon = 1e-6, n_init = 10) {
  # Validate inputs.
  positive_integer <- function(x) {
    is.numeric(x) &&
      length(x) == 1L &&
      is.finite(x) &&
      x >= 1 &&
      x == floor(x)
  }

  if (
    !is.numeric(dates) ||
    length(dates) == 0L ||
    any(!is.finite(dates))
  ) {
    stop("dates must be a nonempty, finite numeric vector.")
  }

  n <- length(dates)

  if (
    !is.numeric(std_errors) ||
    length(std_errors) != n ||
    any(!is.finite(std_errors)) ||
    any(std_errors <= 0)
  ) {
    stop("std_errors must contain one finite, positive value per date.")
  }

  if (!positive_integer(k) || k > n) {
    stop("k must be a positive integer no greater than length(dates).")
  }

  if (!positive_integer(max_iter) || !positive_integer(n_init)) {
    stop("max_iter and n_init must be positive integers.")
  }

  if (
    !is.numeric(epsilon) ||
    length(epsilon) != 1L ||
    !is.finite(epsilon) ||
    epsilon <= 0
  ) {
    stop("epsilon must be finite and positive.")
  }

  if (
    !is.null(cluster_means) &&
    (
      !is.numeric(cluster_means) ||
      length(cluster_means) != k ||
      any(!is.finite(cluster_means))
    )
  ) {
    stop("cluster_means must be NULL or a finite numeric vector of length k.")
  }

  # E-step and observed-data log-likelihood.
  expectation <- function(means) {
    log_density <- matrix(NA_real_, nrow = n, ncol = k)

    for (j in seq_len(k)) {
      log_density[, j] <- stats::dnorm(
        dates,
        mean = means[j],
        sd = std_errors,
        log = TRUE
      )
    }

    # Row-wise log-sum-exp.
    row_max <- apply(log_density, 1L, max)

    if (any(!is.finite(row_max))) {
      stop(
        "Numerical overflow in EM: check the scale of dates and std_errors."
      )
    }

    shifted <- sweep(log_density, 1L, row_max, "-")
    log_normalizer <- row_max + log(rowSums(exp(shifted)))

    log_responsibilities <- sweep(
      log_density, 1L, log_normalizer, "-"
    )

    # Equal mixing probabilities contribute -log(k) per individual.
    log_likelihood <- sum(log_normalizer - log(k))

    if (!is.finite(log_likelihood)) {
      stop("Nonfinite EM log-likelihood: check the input scale.")
    }

    list(
      log_responsibilities = log_responsibilities,
      log_likelihood = log_likelihood
    )
  }

  log_precision <- -2 * log(std_errors)

  best_fit <- NULL
  start_log_likelihoods <- rep(NA_real_, n_init)
  start_converged <- rep(FALSE, n_init)

  for (init in seq_len(n_init)) {
    # Explicit means, if supplied, are used for the first start only.
    # Every additional start independently samples observed date indices.
    if (init == 1L && !is.null(cluster_means)) {
      means <- as.numeric(cluster_means)
    } else {
      means <- dates[sample.int(n, size = k, replace = FALSE)]
    }

    current <- expectation(means)
    trace <- numeric(max_iter + 1L)
    trace[1L] <- current$log_likelihood
    converged <- FALSE

    for (iter in seq_len(max_iter)) {
      new_means <- means

      # Precision- and responsibility-weighted M-step.
      for (j in seq_len(k)) {
        log_weights <-
          current$log_responsibilities[, j] + log_precision

        max_log_weight <- max(log_weights)

        # If a component has numerically zero support, retain its mean.
        # Any empty hard-assignment cluster is repaired by the caller.
        if (!is.finite(max_log_weight)) {
          next
        }

        weights <- exp(log_weights - max_log_weight)
        weights <- weights / sum(weights)

        new_means[j] <- sum(weights * dates)
      }

      # Evaluate responsibilities and likelihood at the UPDATED means.
      updated <- expectation(new_means)

      improvement <-
        updated$log_likelihood - current$log_likelihood

      # Allow only small decreases attributable to floating-point rounding.
      rounding_tolerance <-
        100 * .Machine$double.eps *
        (1 + abs(current$log_likelihood))

      if (improvement < -rounding_tolerance) {
        stop(
          "EM log-likelihood decreased beyond numerical tolerance."
        )
      }

      means <- new_means
      current <- updated
      trace[iter + 1L] <- current$log_likelihood

      # Absolute change in the actual observed-data log-likelihood.
      if (abs(improvement) <= epsilon) {
        converged <- TRUE
        break
      }
    }

    start_log_likelihoods[init] <- current$log_likelihood
    start_converged[init] <- converged

    if (
      is.null(best_fit) ||
      current$log_likelihood > best_fit$log_likelihood
    ) {
      best_fit <- list(
        cluster_means = means,
        responsibilities = exp(current$log_responsibilities),
        log_likelihood = current$log_likelihood,
        log_likelihood_trace = trace[seq_len(iter + 1L)],
        converged = converged,
        iterations = iter,
        best_start = init
      )
    }
  }

  best_fit$start_log_likelihoods <- start_log_likelihoods
  best_fit$start_converged <- start_converged


  return(best_fit)
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

    # 2. Preserve EM assignments, moving one individual per empty cluster.
    if (k_current > 1) {
      cluster_counts <- tabulate(
        cluster_assignments,
        nbins = k_current
      )

      empty_clusters <- which(cluster_counts == 0L)

      for (empty_cluster in empty_clusters) {
        # Only move individuals whose current cluster has >1 member.
        eligible_inds <- which(
          cluster_counts[cluster_assignments] > 1L
        )

        if (length(eligible_inds) == 0L) {
          stop("Cannot fill empty clusters without emptying a donor cluster.")
        }

        # sample.int avoids sample()'s special behaviour for one numeric value.
        moved_ind <- eligible_inds[
          sample.int(length(eligible_inds), 1L)
        ]

        donor_cluster <- cluster_assignments[moved_ind]
        cluster_assignments[moved_ind] <- empty_cluster

        # Update counts before repairing the next empty cluster.
        cluster_counts[donor_cluster] <-
          cluster_counts[donor_cluster] - 1L

        cluster_counts[empty_cluster] <- 1L

        cat(sprintf(
          "Chain %d: moved individual %d from cluster %d to empty cluster %d.\n",
          chain, moved_ind, donor_cluster, empty_cluster
        ))
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
