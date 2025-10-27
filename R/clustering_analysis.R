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

  # 1. Get base thresholds for the method
  if (method %in% names(thresholds_config)) {
    base_thresholds <- thresholds_config[[method]]
    cat(sprintf("\n=== THRESHOLD-BASED K SELECTION (%s) ===\n", method))
  } else {
    base_thresholds <- thresholds_config[["GLOBETROTTER"]]
    cat(sprintf("\n=== THRESHOLD-BASED K SELECTION (DEFAULT - GLOBETROTTER) ===\n"))
    cat("Warning: Unknown method '", method, "', using GLOBETROTTER thresholds\n")
  }

  # 2. Apply CV scaling to improvements
  scaled_improvements <- lapply(improvements$improvements, function(x) x / cv)
  completed_ks       <- improvements$completed_ks

  # 3. Setup cascade thresholds and CV cutoff
  threshold_k1_k2_for_k2 <- base_thresholds$k1_to_k2_for_k2
  threshold_k1_k2_for_k3 <- base_thresholds$k1_to_k2_for_k3
  threshold_k2_k3        <- base_thresholds$k2_to_k3
  threshold_k3_k4        <- base_thresholds$k3_to_k4
  clustering_strength_threshold <- base_thresholds$clustering_strength
  cv_cutoff <- base_thresholds$cv_threshold %||% 1

  # 4. Helper for checking improvements
  check_improvement <- function(name, threshold) {
    if (!is.null(scaled_improvements[[name]])) {
      scaled_improvements[[name]] > threshold
    } else {
      warning(sprintf("No %s improvement available; skipping test", name))
      FALSE
    }
  }

  # Helper to print the comparison for a given likelihood jump
  print_lik_comparison <- function(name, threshold) {
    val <- scaled_improvements[[name]]
    if (is.null(val) || is.na(val)) {
      status <- "FAIL"
      val_print <- NA_real_
    } else {
      val_print <- val
      status <- if (val > threshold) "PASS" else "FAIL"
    }
    cat(sprintf("  %s: %.3f vs threshold %.3f [%s]\n",
                name, val_print, threshold, status))
  }

  # 5. Handle missing k=1 case
  if (!"1" %in% completed_ks) {
    if (length(completed_ks) == 0) {
      stop("No feasible solutions found for any k values.")
    }
    min_k <- as.numeric(min(as.numeric(completed_ks)))
    cat(sprintf("No feasible solutions for k=1. Minimum feasible k: %g\n", min_k))

    recommended_k <- "2"
    if ("2" %in% completed_ks && !is.null(all_mcmc_results[["2"]]$co_clustering_matrix)) {
      k2_str <- clustering_strength(all_mcmc_results[["2"]]$co_clustering_matrix)
      cat("\n=== CLUSTERING STRENGTH ANALYSIS ===\n")
      cat(sprintf("K=2 strength: %.3f vs threshold %.2f [%s]\n",
                  k2_str, clustering_strength_threshold,
                  ifelse(k2_str > clustering_strength_threshold, "PASS", "FAIL")))
    } else {
      k2_str <- NA
    }

    # Continue to cascade even if k=1 missing
    use_lik <- cv >= cv_cutoff
    if (use_lik) {
      cat("\n=== LIKELIHOOD CASCADE (k=1 missing, using k>=2) ===\n")
      print_lik_comparison("k2_to_k3", threshold_k2_k3)
      print_lik_comparison("k1_to_k2", threshold_k1_k2_for_k3)

      if (check_improvement("k2_to_k3", threshold_k2_k3) &&
          check_improvement("k1_to_k2", threshold_k1_k2_for_k3)) {
        recommended_k <- "3"
        cat("=> Both criteria pass: selecting k=3\n")

        print_lik_comparison("k3_to_k4", threshold_k3_k4)
        if (check_improvement("k3_to_k4", threshold_k3_k4)) {
          recommended_k <- "4+"
          cat("=> k3->k4 passes: selecting k=4+\n")
        } else {
          cat("=> k3->k4 fails: stopping at k=3\n")
        }
      } else {
        cat("=> Criteria not met: stopping at k=2\n")
      }
    }

    cat(sprintf("\n*** RECOMMENDED K: %s ***\n", recommended_k))
    return(list(
      recommended_k          = recommended_k,
      thresholds_used        = base_thresholds,
      clustering_strength_k2 = k2_str
    ))
  }

  # 6. Normal workflow when k=1 exists
  cat(sprintf("CV: %.3f (cutoff %.3f)\n", cv, cv_cutoff))
  recommended_k <- "1"

  # 7. K=2 clustering strength check
  if ("2" %in% names(all_mcmc_results) && !is.null(all_mcmc_results[["2"]]$co_clustering_matrix)) {
    k2_str <- clustering_strength(all_mcmc_results[["2"]]$co_clustering_matrix)
    cat("\n=== CLUSTERING STRENGTH ANALYSIS ===\n")
    cat(sprintf("K=2 strength: %.3f vs threshold %.2f [%s]\n",
                k2_str, clustering_strength_threshold,
                ifelse(k2_str > clustering_strength_threshold, "PASS", "FAIL")))

    use_likelihood <- cv >= cv_cutoff

    # 7a. Decide on k=2
    if (!use_likelihood) {
      # Low CV: use clustering strength only
      if (k2_str > clustering_strength_threshold) {
        recommended_k <- "2"
        cat("=> CV < cutoff: using clustering strength only => k=2\n")
      } else {
        cat("=> CV < cutoff: clustering strength fails => stay at k=1\n")
      }
    } else {
      # High CV: use clustering strength OR likelihood
      cat("\n=== LIKELIHOOD TEST FOR K=2 ===\n")
      print_lik_comparison("k1_to_k2", threshold_k1_k2_for_k2)

      if (k2_str > clustering_strength_threshold ||
          check_improvement("k1_to_k2", threshold_k1_k2_for_k2)) {
        recommended_k <- "2"
        cat("=> Clustering strength or k1->k2 passes => k=2\n")
      } else {
        cat("=> Neither criterion met => stay at k=1\n")
      }
    }

    # 8. Cascade to k>=3 if k=2 selected
    if (recommended_k == "2" && use_likelihood) {
      cat("\n=== LIKELIHOOD CASCADE FOR K>=3 ===\n")
      print_lik_comparison("k2_to_k3", threshold_k2_k3)
      print_lik_comparison("k1_to_k2", threshold_k1_k2_for_k3)

      if (check_improvement("k2_to_k3", threshold_k2_k3) &&
          check_improvement("k1_to_k2", threshold_k1_k2_for_k3)) {
        recommended_k <- "3"
        cat("=> Both k2->k3 and k1->k2 (for k3) pass => k=3\n")

        print_lik_comparison("k3_to_k4", threshold_k3_k4)
        if (check_improvement("k3_to_k4", threshold_k3_k4)) {
          recommended_k <- "4+"
          cat("=> k3->k4 passes => k=4+\n")
        } else {
          cat("=> k3->k4 fails => stop at k=3\n")
        }
      } else {
        cat("=> k2->k3 criteria not met => stop at k=2\n")
      }
    }
  } else {
    warning("No k=2 co-clustering matrix; falling back to likelihood only")
    print_lik_comparison("k1_to_k2", threshold_k1_k2_for_k2)
    if (check_improvement("k1_to_k2", threshold_k1_k2_for_k2)) {
      recommended_k <- "2"
      cat("=> Likelihood test passes => k=2\n")
    }
  }

  cat(sprintf("\n*** RECOMMENDED K: %s ***\n", recommended_k))
  if (grepl("^[0-9]+$", recommended_k) && as.integer(recommended_k) >= 4) {
    recommended_k <- "4+"
  }

  return(list(
    recommended_k          = recommended_k,
    thresholds_used        = base_thresholds,
    clustering_strength_k2 = if (exists("k2_str")) k2_str else NA,
    decision_log           = list(
      k2_strength      = if (exists("k2_str")) k2_str else NA,
      pass_k1_to_k2    = check_improvement("k1_to_k2", threshold_k1_k2_for_k2),
      pass_k2_to_k3    = check_improvement("k2_to_k3", threshold_k2_k3),
      pass_k3_to_k4    = check_improvement("k3_to_k4", threshold_k3_k4)
    )
  ))
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
