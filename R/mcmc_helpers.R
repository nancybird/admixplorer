# All the missing helper functions for the MCMC function:

#' Handle k=1 case
#' @keywords internal
handle_k1_case <- function(dates, std_errors, sampling_ages, n, lambda) {
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
    joint_date_mean = rep(weighted_mean, n),
    joint_date_lower_mean = rep(ci[1], n),
    joint_date_upper_mean = rep(ci[2], n),
    joint_date_best = rep(weighted_mean, n),
    joint_date_lower_best = rep(ci[1], n),
    joint_date_upper_best = rep(ci[2], n),
    weighted_joint_date_lower = rep(ci[1], n),
    weighted_joint_date_upper = rep(ci[2], n),
    combined_joint_date_lower_mean = rep(ci[1], n),
    combined_joint_date_upper_mean = rep(ci[2], n),
    combined_joint_date_lower_best = rep(ci[1], n),
    combined_joint_date_upper_best = rep(ci[2], n),
    joint_date_sd_weighted = rep(se_weighted, n),
    sample_ages_sd_mcmc = rep(se_weighted, n),
    joint_date_sd_combined = rep(se_weighted, n),
    mean_log_likelihood_lambda1 = rep(final_ll, n),
    mean_log_likelihood_lambda10 = rep(final_ll, n),
    final_log_likelihood_mean_sampling_ages = rep(final_ll, n),
    final_log_likelihood_mean_sampling_ages_perind = final_ll_perind,
    final_log_likelihood_best_sampling_ages = rep(final_ll, n),
    final_log_likelihood_best_sampling_ages_perind = final_ll_perind,
    acceptance_rate = rep(NA, n)
  )

  return(list(
    result = result_df,
    co_clustering_matrix = NULL,
    cluster_mean_samples = NULL,
    sampling_age_samples = NULL
  ))
}

#' Calculate joint date statistics for individuals
#' @keywords internal
calculate_joint_date_stats <- function(n, cluster_assignments, cluster_means, cluster_sds) {
  joint_date_stats <- sapply(1:n, function(i) {
    cl <- cluster_assignments[i]
    mean_joint <- cluster_means[cl]
    se_joint <- cluster_sds[cl]
    lower_joint <- mean_joint - 1.96 * se_joint
    upper_joint <- mean_joint + 1.96 * se_joint
    c(mean_joint, lower_joint, upper_joint)
  })
  return(joint_date_stats)
}

#' Calculate final likelihoods with different sampling age approaches
#' @keywords internal
calculate_final_likelihoods <- function(cluster_assignments, final_cluster_means,
                                        best_cluster_means_all, std_errors, dates,
                                        age_summaries, best_sampling_ages, n) {

  # Likelihood with mean sampling ages
  ll_with_mean_samples <- nig_log_likelihood(
    cluster_assignments = cluster_assignments,
    cluster_means = final_cluster_means,
    std_errors = std_errors,
    dates = dates + sapply(age_summaries, `[[`, "mean"),
    n = n
  )$total

  ll_with_mean_samples_perind <- nig_log_likelihood(
    cluster_assignments = cluster_assignments,
    cluster_means = final_cluster_means,
    std_errors = std_errors,
    dates = dates + sapply(age_summaries, `[[`, "mean"),
    n = n
  )$per_individual

  # Likelihood with best sampling ages
  ll_with_best_samples <- nig_log_likelihood(
    cluster_assignments = cluster_assignments,
    cluster_means = best_cluster_means_all,
    std_errors = std_errors,
    dates = dates + best_sampling_ages,
    n = n
  )$total

  ll_with_best_samples_perind <- nig_log_likelihood(
    cluster_assignments = cluster_assignments,
    cluster_means = best_cluster_means_all,
    std_errors = std_errors,
    dates = dates + best_sampling_ages,
    n = n
  )$per_individual

  return(list(
    ll_mean = ll_with_mean_samples,
    ll_mean_perind = ll_with_mean_samples_perind,
    ll_best = ll_with_best_samples,
    ll_best_perind = ll_with_best_samples_perind
  ))
}

#' Calculate uncertainties for sampling ages and joint dates
#' @importFrom stats sd
#' @keywords internal
calculate_uncertainties <- function(n, k, cluster_assignments, sampling_age_samples,
                                    cluster_mean_samples, final_cluster_sds, sample_ages) {

  # Individual age standard deviations
  if (sample_ages) {
    individual_age_sd <- apply(sampling_age_samples, 2, stats::sd)
  } else {
    individual_age_sd <- rep(0, n)
  }

  tiny_constant <- 1e-6

  # Joint date SD from MCMC
  joint_date_sd_mcmc <- sapply(1:k, function(cl) {
    mean_samples_for_cluster <- cluster_mean_samples[, cl]
    stats::sd(mean_samples_for_cluster, na.rm = TRUE)
  })

  joint_date_sd_mcmc_individual <- sapply(1:n, function(i) {
    cl <- cluster_assignments[i]
    joint_date_sd_mcmc[cl]
  })

  # Joint date SD from weighted estimation
  joint_date_sd_weighted <- sapply(1:n, function(i) {
    cl <- cluster_assignments[i]
    final_cluster_sds[cl]
  })

  # Sample ages SD by cluster
  sample_ages_sd_by_cluster <- sapply(1:k, function(cl) {
    inds_in_cluster <- which(cluster_assignments == cl)
    if (length(inds_in_cluster) > 0) {
      cluster_individual_age_sd <- individual_age_sd[inds_in_cluster]
      cluster_individual_age_sd[cluster_individual_age_sd == 0] <- tiny_constant
      sqrt(1 / sum(1 / cluster_individual_age_sd^2))
    } else {
      NA
    }
  })

  # Combined joint date SD
  joint_date_sd_combined <- sapply(1:n, function(i) {
    cl <- cluster_assignments[i]
    sqrt(joint_date_sd_weighted[i]^2 + sample_ages_sd_by_cluster[cl]^2)
  })

  # Sample ages SD per individual
  sample_ages_sd_mcmc <- sapply(1:n, function(i) {
    cl <- cluster_assignments[i]
    sample_ages_sd_by_cluster[cl]
  })

  return(list(
    joint_date_sd_weighted = joint_date_sd_weighted,
    sample_ages_sd_mcmc = sample_ages_sd_mcmc,
    joint_date_sd_combined = joint_date_sd_combined
  ))
}

#' Print time traveler summary
#' @keywords internal
print_time_traveler_summary <- function(time_traveler_log, k) {
  if (length(time_traveler_log) > 0) {
    # Count how often each individual was a time traveler
    all_time_travelers <- unlist(lapply(time_traveler_log, function(x) x$individuals))
    tt_counts <- sort(table(all_time_travelers), decreasing = TRUE)

    warning("Time travelers detected in MCMC iterations. See co_clustering_matrix for details.")
    print(head(tt_counts, 5))

    # Suggest action
    most_problematic <- as.numeric(names(tt_counts)[1:min(3, length(tt_counts))])
    cat("SUGGESTION: Remove individual(s)", paste(most_problematic, collapse = ", "),
        "OR increase k to", k+1, "\n")
  }
}

#' Create final result dataframe
#' @keywords internal
create_result_dataframe <- function(n, cluster_assignments, age_summaries, best_sampling_ages,
                                    joint_date_stats, joint_date_stats_best, uncertainty_results,
                                    ll_results, acceptance_rate, mean_log_likelihood_lambda1,
                                    mean_log_likelihood_lambda10, combined_stats, combined_stats_best) {

  result_df <- data.frame(
    individual = 1:n,
    cluster = cluster_assignments,

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

    joint_date_sd_weighted = uncertainty_results$joint_date_sd_weighted,
    sample_ages_sd_mcmc = uncertainty_results$sample_ages_sd_mcmc,
    joint_date_sd_combined = uncertainty_results$joint_date_sd_combined,

    # Likelihoods (INCLUDING LAMBDA-SPECIFIC ONES)
    mean_log_likelihood_lambda1 = rep(mean_log_likelihood_lambda1, n),
    mean_log_likelihood_lambda10 = rep(mean_log_likelihood_lambda10, n),
    final_log_likelihood_mean_sampling_ages = rep(ll_results$ll_mean, n),
    final_log_likelihood_mean_sampling_ages_perind = ll_results$ll_mean_perind,
    final_log_likelihood_best_sampling_ages = rep(ll_results$ll_best, n),
    final_log_likelihood_best_sampling_ages_perind = ll_results$ll_best_perind,

    # MCMC
    acceptance_rate = rep(acceptance_rate, n)
  )

  return(result_df)
}

#' Process MCMC results after sampling (MAIN PROCESSING FUNCTION)
#' @importFrom stats sd quantile
#' @keywords internal
process_mcmc_results <- function(n, k, dates, std_errors, best_cluster_assignments,
                                 best_sampling_ages, best_cluster_means, sampling_age_samples,
                                 cluster_mean_samples, co_clustering_matrix, acceptance_rate,
                                 num_iter, burn_in, thin, sample_ages, age_ranges,
                                 time_traveler_log, mean_log_likelihood_lambda1,
                                 mean_log_likelihood_lambda10) {

  get_ci <- function(x) {
    q <- stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
    list(mean = mean(x, na.rm = TRUE), lower = q[1], upper = q[2])
  }

  acceptance_rate <- acceptance_rate / num_iter
  co_clustering_matrix <- co_clustering_matrix / ((num_iter - burn_in)/thin)

  # Summarize sampling ages
  age_summaries <- if (sample_ages) {
    apply(sampling_age_samples, 2, get_ci)
  } else {
    lapply(best_sampling_ages, function(x) list(mean = x, lower = x, upper = x))
  }

  # Calculate final cluster means and statistics
  final_cluster_means <- numeric(k)
  final_cluster_sds <- numeric(k)
  best_cluster_means_all <- numeric(k)

  for (j in 1:k) {
    inds_in_cluster <- which(best_cluster_assignments == j)
    if (length(inds_in_cluster) > 0) {
      # Using mean sampling ages
      cluster_dates <- dates[inds_in_cluster] + sapply(age_summaries[inds_in_cluster], `[[`, "mean")
      cluster_errors <- std_errors[inds_in_cluster]
      weights <- 1 / cluster_errors^2
      final_cluster_means[j] <- sum(cluster_dates * weights) / sum(weights)
      final_cluster_sds[j] <- sqrt(1 / sum(weights))

      # Using best sampling ages
      best_cluster_dates <- dates[inds_in_cluster] + best_sampling_ages[inds_in_cluster]
      best_cluster_means_all[j] <- sum(best_cluster_dates * weights) / sum(weights)

      # Final check for time travelers
      if (!is.null(age_ranges)) {
        time_travelers <- which(final_cluster_means[j] < age_ranges[inds_in_cluster, 1])
        if (length(time_travelers) > 0) {
          time_traveler_inds <- inds_in_cluster[time_travelers]
          warning(paste("Time travelers detected in final results! Individual(s)",
                        paste(time_traveler_inds, collapse = ", "),
                        "have joint admixture date more recent than sampling age in cluster", j,
                        ". We suggest removing individual(s)", paste(time_traveler_inds, collapse = ", "),
                        "OR increasing k to", k+1))
        }
      }
    } else {
      final_cluster_means[j] <- NA
      final_cluster_sds[j] <- NA
      best_cluster_means_all[j] <- NA
    }
  }

  # Calculate joint date statistics
  joint_date_stats <- calculate_joint_date_stats(n, best_cluster_assignments, final_cluster_means, final_cluster_sds)
  joint_date_stats_best <- calculate_joint_date_stats(n, best_cluster_assignments, best_cluster_means_all, final_cluster_sds)

  # Calculate likelihoods
  ll_results <- calculate_final_likelihoods(best_cluster_assignments, final_cluster_means,
                                            best_cluster_means_all, std_errors, dates,
                                            age_summaries, best_sampling_ages, n)

  # Calculate uncertainties
  uncertainty_results <- calculate_uncertainties(n, k, best_cluster_assignments, sampling_age_samples,
                                                  cluster_mean_samples, final_cluster_sds, sample_ages)

  # Combined uncertainty statistics
  combined_stats <- sapply(1:n, function(i) {
    mean_joint <- joint_date_stats[1, i]
    lower_combined <- mean_joint - 1.96 * uncertainty_results$joint_date_sd_combined[i]
    upper_combined <- mean_joint + 1.96 * uncertainty_results$joint_date_sd_combined[i]
    c(lower_combined, upper_combined)
  })

  combined_stats_best <- sapply(1:n, function(i) {
    mean_joint <- joint_date_stats_best[1, i]
    lower_combined <- mean_joint - 1.96 * uncertainty_results$joint_date_sd_combined[i]
    upper_combined <- mean_joint + 1.96 * uncertainty_results$joint_date_sd_combined[i]
    c(lower_combined, upper_combined)
  })

  # Print time traveler summary
  print_time_traveler_summary(time_traveler_log, k)

  # Create final result dataframe
  result_df <- create_result_dataframe(n, best_cluster_assignments, age_summaries, best_sampling_ages,
                                       joint_date_stats, joint_date_stats_best, uncertainty_results,
                                       ll_results, acceptance_rate, mean_log_likelihood_lambda1,
                                       mean_log_likelihood_lambda10, combined_stats, combined_stats_best)

  co_clustering_df <- as.data.frame(co_clustering_matrix)
  rownames(co_clustering_df) <- colnames(co_clustering_df) <- paste0("ind", 1:n)

  return(list(
    result = result_df,
    co_clustering_matrix = co_clustering_df,
    cluster_mean_samples = cluster_mean_samples,
    sampling_age_samples = if (sample_ages) sampling_age_samples else NULL
  ))
}
