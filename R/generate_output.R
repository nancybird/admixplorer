#' Generate final output
#'
#' @param all_mcmc_results List of MCMC results
#' @param recommended_k Recommended k value
#' @param original_data Original filtered data
#' @param pop.vec Population vector
#' @param dates_original Original dates
#' @param sample_age_est Whether sample ages were estimated
#' @param outfile Output file prefix
#' @param method Method used
#' @param thresholds_info Threshold information from analysis
#' @return Combined output dataframe
#' @export
#' @importFrom utils write.table
generate_final_output <- function(all_mcmc_results, recommended_k, original_data, pop.vec,
                                  dates_original, sample_age_est, outfile, method,
                                  thresholds_info) {

  # Build final output for all completed k values
  output_list <- list()
  completed_ks <- names(all_mcmc_results)[!sapply(all_mcmc_results, is.null)]

  for (k in completed_ks) {
    k_num <- as.numeric(k)
    mcmc_result <- all_mcmc_results[[k]]

    # Extract results
    likelihood_mean <- mcmc_result$result$final_log_likelihood_best_sampling_ages[1]
    likelihood_per_ind <- mcmc_result$result$final_log_likelihood_best_sampling_ages_perind

    # Add recommendation flag
    is_recommended <- (k == recommended_k) || (recommended_k == "4+" && k_num >= 4)

    # Build output dataframe - SIMPLE VERSION
    output <- data.frame(
      model_k = rep(k_num, length(pop.vec)),
      pop = pop.vec,
      cluster = if (k_num == 1) rep(1, length(dates_original)) else mcmc_result$result$cluster,
      joint_date_est_best = mcmc_result$result$joint_date_best,
      joint_date_lowerci = mcmc_result$result$combined_joint_date_lower_best,
      joint_date_upperci = mcmc_result$result$combined_joint_date_upper_best,
      likelihood_mean = likelihood_mean,
      likelihood_per_ind = likelihood_per_ind,
      recommended_k = recommended_k,  # Just the recommended k
      is_recommended = is_recommended,  # Simple boolean flag
      ind_date_est = dates_original + mcmc_result$result$sampling_age_best,
      stringsAsFactors = FALSE
    )

    # Add sampling ages if estimated
    if (sample_age_est == TRUE) {
      output$sample_age_est_best <- mcmc_result$result$sampling_age_best
      output$sample_age_lowerci_mean <- mcmc_result$result$sampling_age_mean_lower
      output$sample_age_upperci_mean <- mcmc_result$result$sampling_age_mean_upper
    }

    output_list[[k]] <- output
  }

  # Combine and write output
  combined_output <- do.call(rbind, output_list)

  # Enhanced summary header
  summary_header <- paste0(
    "# LIKELIHOOD-BASED K SELECTION RESULTS (", toupper(method), ")\n",
    "# Recommended K: ", recommended_k, "\n",
    "# Method: ", method, "\n",
    "# Scaled thresholds: k1 to k2(k=2):", round(thresholds_info$thresholds_used$k1_to_k2_for_k2, 3),
    " k1 to k2(k=3):", round(thresholds_info$thresholds_used$k1_to_k2_for_k3, 3),
    " k2 to k3:", round(thresholds_info$thresholds_used$k2_to_k3, 3),
    " k3 to k4:", round(thresholds_info$thresholds_used$k3_to_k4, 3), "\n",
    "#\n"
  )

  # Write output with header
  output.outfile <- paste0(outfile, ".output.txt")
  cat(summary_header, file = output.outfile)
  utils::write.table(combined_output, output.outfile, row.names = FALSE, col.names = TRUE,
                     quote = FALSE, append = TRUE)

  cat("Analysis complete. Results written to", output.outfile, "\n")
  cat("File includes recommended k =", recommended_k, "based on likelihood thresholds\n")

  return(combined_output)
}
