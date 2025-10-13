# tests/testthat/test-threshold-selection.R (COMPLETE CORRECTED VERSION)
# Test the threshold selection logic

# Fixed create_mock_mcmc_results function that actually creates the strength you want:

create_mock_mcmc_results <- function(k2_clustering_strength = 0.8) {
  # Create a proper co-clustering matrix that will give the exact clustering strength specified
  
  # For a 5x5 matrix with 2 clusters (first 3 vs last 2), we need:
  # - Within-cluster probabilities = high (close to 1)
  # - Between-cluster probabilities = low (close to 0)
  
  # Calculate what the within/between probabilities should be to get target strength
  # Clustering strength = mean(abs(values - 0.5) * 2)
  # For perfect separation: within = 1, between = 0, strength = 1
  # For no separation: within = between = 0.5, strength = 0
  
  within_cluster_prob <- 0.5 + (k2_clustering_strength / 2)
  between_cluster_prob <- 0.5 - (k2_clustering_strength / 2)
  
  # Create the matrix with proper structure
  co_clustering_matrix <- matrix(c(
    1, within_cluster_prob, within_cluster_prob, between_cluster_prob, between_cluster_prob,
    within_cluster_prob, 1, within_cluster_prob, between_cluster_prob, between_cluster_prob,
    within_cluster_prob, within_cluster_prob, 1, between_cluster_prob, between_cluster_prob,
    between_cluster_prob, between_cluster_prob, between_cluster_prob, 1, within_cluster_prob,
    between_cluster_prob, between_cluster_prob, between_cluster_prob, within_cluster_prob, 1
  ), nrow = 5, byrow = TRUE)
  
  list(
    "1" = list(
      result = list(
        final_log_likelihood_best_sampling_ages = c(-100),
        final_log_likelihood_best_sampling_ages_perind = c(-10, -10, -10, -10, -10),
        cluster = rep(1, 5)
      )
    ),
    "2" = list(
      result = list(
        final_log_likelihood_best_sampling_ages = c(-80),
        final_log_likelihood_best_sampling_ages_perind = c(-8, -8, -8, -8, -8),
        cluster = c(1, 1, 1, 2, 2)
      ),
      co_clustering_matrix = co_clustering_matrix
    )
  )
}



create_mock_improvements <- function() {
  list(
    completed_ks = c("1", "2"),
    likelihoods = c("1" = -100, "2" = -80),
    likelihoods_per_ind = c("1" = -20, "2" = -16),
    improvements = list(
      k1_to_k2 = 4  # Per individual improvement
    )
  )
}

test_that("apply_threshold_selection recommends k=1 for weak clustering", {
  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.5)  # Weak
  mock_improvements <- create_mock_improvements()
  
  result <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER", 
    cv = 2,
    all_mcmc_results = mock_results
  )
  
  expect_equal(result$recommended_k, "1")
  expect_true(result$clustering_strength_k2 < 0.77)
})

test_that("apply_threshold_selection recommends k=2 for strong clustering", {
  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.9)  # Strong
  mock_improvements <- create_mock_improvements()
  
  result <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER",
    cv = 2, 
    all_mcmc_results = mock_results
  )
  
  expect_equal(result$recommended_k, "2")
  expect_true(result$clustering_strength_k2 > 0.77)
})

test_that("apply_threshold_selection handles low CV correctly", {
  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.9)
  mock_improvements <- create_mock_improvements()
  
  # With low CV, should ignore likelihood thresholds
  result <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER",
    cv = 0.5,  # Low CV
    all_mcmc_results = mock_results
  )
  
  # Should recommend k=2 based on clustering strength alone
  expect_equal(result$recommended_k, "2")
})

test_that("apply_threshold_selection uses correct method thresholds", {
  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.8)
  mock_improvements <- create_mock_improvements()
  
  # Test DATES method with CORRECT expected values
  result_dates <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "DATES",
    cv = 2,
    all_mcmc_results = mock_results
  )
  
  # Use the actual values from your thresholds_config
  expect_equal(result_dates$thresholds_used$k1_to_k2_for_k2, 0.97)  # CORRECTED
  expect_equal(result_dates$thresholds_used$k1_to_k2_for_k3, 4.05)  # CORRECTED
  expect_equal(result_dates$thresholds_used$k2_to_k3, 0.58)         # CORRECTED
  expect_equal(result_dates$thresholds_used$k3_to_k4, 0.42)         # CORRECTED
  expect_equal(result_dates$thresholds_used$clustering_strength, 0.99)
  
  # Test GLOBETROTTER method  
  result_gt <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER",
    cv = 2,
    all_mcmc_results = mock_results
  )
  
  expect_equal(result_gt$thresholds_used$k1_to_k2_for_k2, 0.65)
  expect_equal(result_gt$thresholds_used$k1_to_k2_for_k3, 1.21)
  expect_equal(result_gt$thresholds_used$k2_to_k3, 0.34)
  expect_equal(result_gt$thresholds_used$k3_to_k4, 0.18)
  expect_equal(result_gt$thresholds_used$clustering_strength, 0.77)
})

test_that("apply_threshold_selection handles missing k=2 results", {
  # No k=2 results available
  mock_results <- list(
    "1" = list(
      result = list(
        final_log_likelihood_best_sampling_ages = c(-100),
        cluster = rep(1, 5)
      )
    )
  )
  mock_improvements <- list(
    improvements = list(k1_to_k2 = 0.5)  # Low improvement
  )
  
  result <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER",
    cv = 2,
    all_mcmc_results = mock_results
  )
  
  # Should fall back to likelihood-based decision
  expect_equal(result$recommended_k, "1")  # Low improvement should give k=1
  expect_true(is.na(result$clustering_strength_k2))
})

test_that("DATES method uses higher clustering strength threshold", {
  # Test that DATES method uses 0.99 clustering threshold (stricter than GLOBETROTTER's 0.77)
  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.95)  # Between 0.77 and 0.99
  mock_improvements <- create_mock_improvements()
  
  result_gt <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER",
    cv = 2,
    all_mcmc_results = mock_results
  )
  
  result_dates <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "DATES", 
    cv = 2,
    all_mcmc_results = mock_results
  )
  
  # GLOBETROTTER should recommend k=2 (0.85 > 0.77)
  expect_equal(result_gt$recommended_k, "2")
  
  # DATES should recommend k=1 (0.85 < 0.99)  
  expect_equal(result_dates$recommended_k, "1")
})