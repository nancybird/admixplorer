#' Calculate likelihood improvements between k values
#' 
#' @param all_mcmc_results List of MCMC results from run_clustering_analysis
#' @param n_individuals Number of individuals (for per-individual likelihood)
#' @return List of improvements and related statistics
#' @export
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

#' Calculate FST-based scaling factor
#' 
#' @param input_fst FST between admixing populations
#' @param method Method used (GLOBETROTTER or DATES)
#' @return List with scaling information
#' @keywords internal
calculate_fst_scaling <- function(input_fst, method) {
  if (method == "GLOBETROTTER") {
    baseline_fst <- 0.11  # Spanish-Japanese FST (your training data)
  } else if (method == "DATES") {
    baseline_fst <- 0.11  # Assuming same training populations
  } else {
    baseline_fst <- 0.11  # Default
  }
  
  # Scale thresholds: higher FST = more differentiated = easier to cluster = lower thresholds needed
  # Lower FST = less differentiated = harder to cluster = higher thresholds needed
  scaling_factor <- baseline_fst / input_fst
  
  # Add bounds to prevent extreme scaling
  scaling_factor <- max(0.5, min(2.0, scaling_factor))  # Limit scaling to 0.5x - 2.0x
  
  return(list(
    scaling_factor = scaling_factor,
    baseline_fst = baseline_fst,
    input_fst = input_fst
  ))
}

#' Apply threshold-based k selection
#' 
#' @param improvements List from calculate_likelihood_improvements
#' @param method Method used (GLOBETROTTER or DATES)
#' @param fst FST between admixing populations
#' @param cv Coefficient of variation for scaling
#' @return Recommended k value
#' @export
apply_threshold_selection <- function(improvements, method, fst, cv) {
  cat("\n>>> FINAL ANALYSIS WITH LIKELIHOOD THRESHOLDS <<<\n")
  
  # Use internal thresholds_config data
  # Calculate FST scaling
  fst_scaling <- calculate_fst_scaling(fst, method)
  
  cat(sprintf("FST-based threshold scaling:\n"))
  cat(sprintf("  Input FST: %.3f\n", fst_scaling$input_fst))
  cat(sprintf("  Baseline FST: %.3f\n", fst_scaling$baseline_fst))
  cat(sprintf("  Scaling factor: %.3f\n", fst_scaling$scaling_factor))
  
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
  
  # Get final thresholds (currently not applying FST scaling, but structure is ready)
  threshold_k1_to_k2_for_k2 <- base_thresholds$k1_to_k2_for_k2
  threshold_k1_to_k2_for_k3 <- base_thresholds$k1_to_k2_for_k3
  threshold_k2_to_k3 <- base_thresholds$k2_to_k3
  threshold_k3_to_k4 <- base_thresholds$k3_to_k4
  
  cat(sprintf("Method: %s\n", method))
  cat(sprintf("Thresholds: k1 to k2(k=2): %.2f, k1 to k2(k=3): %.2f, k2 to k3: %.2f, k3 to k4: %.2f\n", 
              threshold_k1_to_k2_for_k2, threshold_k1_to_k2_for_k3, 
              threshold_k2_to_k3, threshold_k3_to_k4))
  
  recommended_k <- "1"  # Start with k=1 (as character to allow "4+")
  
  # Check k1 to k2 (basic threshold for k=2)
  if ("k1_to_k2" %in% names(scaled_improvements)) {
    if (scaled_improvements$k1_to_k2 > threshold_k1_to_k2_for_k2) {
      recommended_k <- "2"
      cat(sprintf("k1 to k2 improvement (%.4f) > threshold for k=2 (%.2f): Consider k=2\n", 
                  scaled_improvements$k1_to_k2, threshold_k1_to_k2_for_k2))
      
      # Check k2 to k3
      if ("k2_to_k3" %in% names(scaled_improvements)) {
        if (scaled_improvements$k2_to_k3 > threshold_k2_to_k3) {
          # For k=3, we need BOTH conditions: higher k1 to k2 AND good k2 to k3
          if (scaled_improvements$k1_to_k2 > threshold_k1_to_k2_for_k3) {
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
            cat(sprintf("k2 to k3 improvement (%.4f) > threshold (%.2f) BUT k1 to k2 (%.4f) <= k=3 threshold (%.2f): Stop at k=2\n", 
                        scaled_improvements$k2_to_k3, threshold_k2_to_k3,
                        scaled_improvements$k1_to_k2, threshold_k1_to_k2_for_k3))
          }
        } else {
          cat(sprintf("k2 to k3 improvement (%.4f) <= threshold (%.2f): Stop at k=2\n", 
                      scaled_improvements$k2_to_k3, threshold_k2_to_k3))
        }
      }
    } else {
      cat(sprintf("k1 to k2 improvement (%.4f) <= threshold (%.2f): Stop at k=1\n", 
                  scaled_improvements$k1_to_k2, threshold_k1_to_k2_for_k2))
    }
  }
  
  cat(sprintf("\n*** RECOMMENDED K: %s ***\n", recommended_k))
  
  list(
    recommended_k = recommended_k,
    thresholds_used = base_thresholds,
    fst_scaling = fst_scaling
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