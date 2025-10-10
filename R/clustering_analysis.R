#' Calculate clustering strength from co-clustering matrix
#'
#' @param coincidence_matrix Co-clustering matrix from MCMC results
#' @return Clustering strength score between 0 and 1
#' @export
clustering_strength <- function(coincidence_matrix) {
  # Get upper triangle indices (avoid double counting)
  upper_tri <- upper.tri(coincidence_matrix)
  values <- coincidence_matrix[upper_tri]
  certainty_scores <- abs(values - 0.5) * 2  # Scale to [0,1]
  return(mean(certainty_scores))
}

#' Apply threshold-based k selection with clustering strength
#'
#' @param improvements List from calculate_likelihood_improvements
#' @param method Method used (GLOBETROTTER or DATES)
#' @param cv Coefficient of variation for scaling
#' @param all_mcmc_results List of all MCMC results (needed for clustering strength)
#' @return Recommended k value
#' @export
apply_threshold_selection <- function(improvements, method, cv, all_mcmc_results) {
  cat("\n>>> FINAL ANALYSIS WITH LIKELIHOOD THRESHOLDS <<<\n")

  # Get base thresholds for the method
  if (method %in% names(thresholds_config)) {
    base_thresholds <- thresholds_config[[method]]
    cat(sprintf("\n=== THRESHOLD-BASED K SELECTION (%s) ===\n", method))
  } else {
    base_thresholds <- thresholds_config[["GLOBETROTTER"]]  # Default
    cat(sprintf("\n=== THRESHOLD-BASED K SELECTION (DEFAULT - GLOBETROTTER) ===\n"))
    cat("Warning: Unknown method '", method, "', using GLOBETROTTER thresholds\n")
  }

  # Apply CV scaling to improvements
  scaled_improvements <- lapply(improvements$improvements, function(x) x / cv)

  threshold_k1_to_k2_for_k2 <- base_thresholds$k1_to_k2_for_k2
  threshold_k1_to_k2_for_k3 <- base_thresholds$k1_to_k2_for_k3
  threshold_k2_to_k3 <- base_thresholds$k2_to_k3
  threshold_k3_to_k4 <- base_thresholds$k3_to_k4

  cat(sprintf("Method: %s\n", method))
  cat(sprintf("Likelihood thresholds: k1 to k2(k=2): %.2f, k1 to k2(k=3): %.2f, k2 to k3: %.2f, k3 to k4: %.2f\n",
              threshold_k1_to_k2_for_k2, threshold_k1_to_k2_for_k3,
              threshold_k2_to_k3, threshold_k3_to_k4))

  recommended_k <- "1"  # Start with k=1

  # NEW: Check clustering strength for k=2 first
  clustering_strength_threshold <- 0.77

  if ("2" %in% names(all_mcmc_results) && !is.null(all_mcmc_results[["2"]]$co_clustering_matrix)) {
    k2_strength <- clustering_strength(all_mcmc_results[["2"]]$co_clustering_matrix)
    cat(sprintf("\n=== CLUSTERING STRENGTH ANALYSIS ===\n"))
    cat(sprintf("K=2 clustering strength: %.3f\n", k2_strength))
    cat(sprintf("Clustering strength threshold: %.2f\n", clustering_strength_threshold))

    if (k2_strength > clustering_strength_threshold) {
      recommended_k <- "2"
      cat(sprintf("K=2 clustering strength (%.3f) > threshold (%.2f): Consider k>=2\n",
                  k2_strength, clustering_strength_threshold))

      # Now use likelihood thresholds for k>=2 decisions
      cat(sprintf("\n=== LIKELIHOOD THRESHOLD ANALYSIS (k>=2) ===\n"))

      # Check k2 to k3
      if ("k2_to_k3" %in% names(scaled_improvements)) {
        if (scaled_improvements$k2_to_k3 > threshold_k2_to_k3) {
          # For k=3, we need BOTH good k2→k3 AND the original k1→k2 was strong enough
          if ("k1_to_k2" %in% names(scaled_improvements) &&
              scaled_improvements$k1_to_k2 > threshold_k1_to_k2_for_k3) {
            recommended_k <- "3"
            cat(sprintf("k2 to k3 improvement (%.4f) > threshold (%.2f) AND k1 to k2 (%.4f) > k=3 threshold (%.2f): Consider k=3\n",
                        scaled_improvements$k2_to_k3, threshold_k2_to_k3,
                        scaled_improvements$k1_to_k2, threshold_k1_to_k2_for_k3))

            # Check k3 to k4
            if ("k3_to_k4" %in% names(scaled_improvements)) {
              if (scaled_improvements$k3_to_k4 > threshold_k3_to_k4) {
                recommended_k <- "4+"  # 4+ because no higher thresholds defined
                cat(sprintf("k3 to k4 improvement (%.4f) > threshold (%.2f): Consider k=4 or higher\n",
                            scaled_improvements$k3_to_k4, threshold_k3_to_k4))
              } else {
                cat(sprintf("k3 to k4 improvement (%.4f) <= threshold (%.2f): Stop at k=3\n",
                            scaled_improvements$k3_to_k4, threshold_k3_to_k4))
              }
            }
          } else {
            cat(sprintf("k2 to k3 improvement (%.4f) > threshold (%.2f) BUT k1 to k2 not strong enough for k=3: Stop at k=2\n",
                        scaled_improvements$k2_to_k3, threshold_k2_to_k3))
          }
        } else {
          cat(sprintf("k2 to k3 improvement (%.4f) <= threshold (%.2f): Stop at k=2\n",
                      scaled_improvements$k2_to_k3, threshold_k2_to_k3))
        }
      }
    } else {
      cat(sprintf("K=2 clustering strength (%.3f) <= threshold (%.2f): Stop at k=1\n",
                  k2_strength, clustering_strength_threshold))
    }
  } else {
    cat("\nWarning: No k=2 results found or no co-clustering matrix available\n")
    cat("Falling back to likelihood-based k=1 vs k=2 decision\n")

    # Fallback to original likelihood method if clustering strength can't be calculated
    if ("k1_to_k2" %in% names(scaled_improvements)) {
      if (scaled_improvements$k1_to_k2 > threshold_k1_to_k2_for_k2) {
        recommended_k <- "2"
        cat(sprintf("k1 to k2 improvement (%.4f) > threshold for k=2 (%.2f): Consider k=2\n",
                    scaled_improvements$k1_to_k2, threshold_k1_to_k2_for_k2))
      }
    }
  }

  cat(sprintf("\n*** RECOMMENDED K: %s ***\n", recommended_k))

  list(
    recommended_k = recommended_k,
    thresholds_used = base_thresholds,
    clustering_strength_k2 = if("2" %in% names(all_mcmc_results)) clustering_strength(all_mcmc_results[["2"]]$co_clustering_matrix) else NA
  )
}



calculate_likelihood_improvements <- function(all_mcmc_results, n_individuals) {
  completed_ks <- names(all_mcmc_results)[!sapply(all_mcmc_results, is.null)]

  likelihoods <- sapply(completed_ks, function(k) {
    all_mcmc_results[[k]]$result$final_log_likelihood_best_sampling_ages[1]
  })

  likelihoods_per_ind <- likelihoods / n_individuals

  # Print likelihoods
  for (k in completed_ks) {
    cat(sprintf("Likelihood k=%s: %.2f \n", k, likelihoods[k]))
  }

  # Calculate improvements between consecutive k values
  improvements <- list()
  completed_ks_num <- sort(as.numeric(completed_ks))

  for (i in 1:(length(completed_ks_num)-1)) {
    k_current <- as.character(completed_ks_num[i])
    k_next <- as.character(completed_ks_num[i+1])

    improvement <- likelihoods_per_ind[k_next] - likelihoods_per_ind[k_current]
    improvements[[paste0("k", k_current, "_to_k", k_next)]] <- improvement

    cat(sprintf("Improvement k%s to k%s: %.4f\n", k_current, k_next, improvement))
  }

  list(
    completed_ks = completed_ks,
    likelihoods = likelihoods,
    likelihoods_per_ind = likelihoods_per_ind,
    improvements = improvements
  )
}


#' Run clustering analysis for all k values
#'
#' @param prepared_data List from prepare_clustering_data
#' @param ks Vector of k values to try
#' @param chains Number of MCMC chains
#' @param sample_age_est Whether to estimate sample ages
#' @param plot Whether to create plots
#' @param outfile Output file prefix
#' @return List of MCMC results for each k
#' @export
run_clustering_analysis <- function(prepared_data, ks, chains, sample_age_est, plot, outfile) {
  print("CLUSTERING STAGE")

  all_mcmc_results <- list()

  for (k_current in ks) {
    result <- process_single_k(k_current, prepared_data, chains, sample_age_est, plot, outfile)
    if (!is.null(result)) {
      all_mcmc_results[[as.character(k_current)]] <- result
    }
  }

  return(all_mcmc_results)
}
