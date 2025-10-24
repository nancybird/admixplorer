# tests/testthat/test-threshold-selection.R (ENHANCED VERSION)

# Helper functions (keep your existing ones)
create_mock_mcmc_results <- function(k2_clustering_strength = 0.8) {
  within_cluster_prob <- 0.5 + (k2_clustering_strength / 2)
  between_cluster_prob <- 0.5 - (k2_clustering_strength / 2)

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
    ),
    "3" = list(
      result = list(
        final_log_likelihood_best_sampling_ages = c(-60),
        final_log_likelihood_best_sampling_ages_perind = c(-6, -6, -6, -6, -6),
        cluster = c(1, 1, 2, 3, 3)
      )
    ),
    "4" = list(
      result = list(
        final_log_likelihood_best_sampling_ages = c(-50),
        final_log_likelihood_best_sampling_ages_perind = c(-5, -5, -5, -5, -5),
        cluster = c(1, 2, 3, 4, 4)
      )
    )
  )
}

create_mock_improvements_weak <- function() {
  list(
    completed_ks = c("1", "2"),
    likelihoods = c("1" = -100, "2" = -99),
    likelihoods_per_ind = c("1" = -20, "2" = -19.8),
    improvements = list(
      k1_to_k2 = 0.2
    )
  )
}

create_mock_improvements_strong <- function() {
  list(
    completed_ks = c("1", "2"),
    likelihoods = c("1" = -100, "2" = -80),
    likelihoods_per_ind = c("1" = -20, "2" = -16),
    improvements = list(
      k1_to_k2 = 4
    )
  )
}

# NEW: Mock improvements for k=3 testing
create_mock_improvements_to_k3 <- function() {
  list(
    completed_ks = c("1", "2", "3"),
    likelihoods = c("1" = -100, "2" = -80, "3" = -60),
    likelihoods_per_ind = c("1" = -20, "2" = -16, "3" = -12),
    improvements = list(
      k1_to_k2 = 4,
      k2_to_k3 = 4
    )
  )
}

# NEW: Mock improvements for k=4+ testing
create_mock_improvements_to_k4 <- function() {
  list(
    completed_ks = c("1", "2", "3", "4"),
    likelihoods = c("1" = -100, "2" = -80, "3" = -60, "4" = -50),
    likelihoods_per_ind = c("1" = -20, "2" = -16, "3" = -12, "4" = -10),
    improvements = list(
      k1_to_k2 = 4,
      k2_to_k3 = 4,
      k3_to_k4 = 2
    )
  )
}

# NEW: Mock for missing k=1 with higher k available
create_mock_improvements_no_k1 <- function() {
  list(
    completed_ks = c("2", "3"),
    likelihoods = c("2" = -80, "3" = -60),
    likelihoods_per_ind = c("2" = -16, "3" = -12),
    improvements = list(
      k1_to_k2 = 4,    # Still provide this for cascade logic
      k2_to_k3 = 4
    )
  )
}

# Original tests (keep all of these)
test_that("apply_threshold_selection recommends k=1 for weak clustering", {
  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.7)
  mock_improvements <- create_mock_improvements_weak()

  result <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER",
    cv = 2,
    all_mcmc_results = mock_results
  )

  expect_equal(result$recommended_k, "1")
  expect_true(result$clustering_strength_k2 < 0.78)
})

test_that("apply_threshold_selection recommends k=2 for strong clustering", {
  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.9)
  mock_improvements <- create_mock_improvements_strong()

  result <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER",
    cv = 2,
    all_mcmc_results = mock_results
  )

  expect_equal(result$recommended_k, "2")
  expect_true(result$clustering_strength_k2 > 0.78)
})

test_that("apply_threshold_selection handles low CV correctly", {
  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.9)
  mock_improvements <- create_mock_improvements_strong()

  result <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER",
    cv = 0.5,
    all_mcmc_results = mock_results
  )

  expect_equal(result$recommended_k, "2")
})



test_that("apply_threshold_selection handles missing k=2 results", {
  mock_results <- list(
    "1" = list(
      result = list(
        final_log_likelihood_best_sampling_ages = c(-100),
        cluster = rep(1, 5)
      )
    )
  )
  mock_improvements <- list(
    completed_ks = c("1"),
    improvements = list()
  )

  result <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER",
    cv = 2,
    all_mcmc_results = mock_results
  )

  expect_equal(result$recommended_k, "1")
})

test_that("DATES method uses higher clustering strength threshold", {
  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.85)
  mock_improvements <- create_mock_improvements_weak()

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

  expect_equal(result_gt$recommended_k, "2")
  expect_equal(result_dates$recommended_k, "1")
})

# ====== NEW TESTS FOR ENHANCED FUNCTIONALITY ======

test_that("apply_threshold_selection cascades to k=3 with strong improvements", {
  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.9)
  mock_improvements <- create_mock_improvements_to_k3()

  result <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER",
    cv = 2,
    all_mcmc_results = mock_results
  )

  expect_equal(result$recommended_k, "3")
  expect_true("decision_log" %in% names(result))
  expect_true(result$decision_log$pass_k2_to_k3)
})

test_that("apply_threshold_selection stops at k=2 when k2_to_k3 weak", {
  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.9)
  mock_improvements <- list(
    completed_ks = c("1", "2", "3"),
    improvements = list(
      k1_to_k2 = 4,
      k2_to_k3 = 0.1  # Below threshold
    )
  )

  result <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER",
    cv = 2,
    all_mcmc_results = mock_results
  )

  expect_equal(result$recommended_k, "2")
  expect_false(result$decision_log$pass_k2_to_k3)
})

test_that("apply_threshold_selection reaches k=4+ with all strong improvements", {
  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.9)
  mock_improvements <- create_mock_improvements_to_k4()

  result <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER",
    cv = 2,
    all_mcmc_results = mock_results
  )

  expect_equal(result$recommended_k, "4+")
  expect_true(result$decision_log$pass_k3_to_k4)
})

test_that("missing k=1 can still cascade to k=3", {
  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.9)
  mock_improvements <- create_mock_improvements_no_k1()

  result <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER",
    cv = 2,
    all_mcmc_results = mock_results
  )

  # Should recommend k=3 even without k=1
  expect_equal(result$recommended_k, "3")
})

test_that("missing k=1 stops at k=2 when k2_to_k3 weak", {
  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.9)
  mock_improvements <- list(
    completed_ks = c("2", "3"),
    improvements = list(
      k1_to_k2 = 4,
      k2_to_k3 = 0.1  # Weak
    )
  )

  result <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER",
    cv = 2,
    all_mcmc_results = mock_results
  )

  expect_equal(result$recommended_k, "2")
})

test_that("decision_log contains all threshold checks", {
  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.9)
  mock_improvements <- create_mock_improvements_to_k4()

  result <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER",
    cv = 2,
    all_mcmc_results = mock_results
  )

  expect_true("decision_log" %in% names(result))
  expect_true("k2_strength" %in% names(result$decision_log))
  expect_true("pass_k1_to_k2" %in% names(result$decision_log))
  expect_true("pass_k2_to_k3" %in% names(result$decision_log))
  expect_true("pass_k3_to_k4" %in% names(result$decision_log))
  expect_equal(result$decision_log$k2_strength, 0.9, tolerance = 0.01)
})

test_that("low CV blocks cascade beyond k=2 even with strong improvements", {
  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.9)
  mock_improvements <- create_mock_improvements_to_k3()

  result <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER",
    cv = 0.5,  # Low CV
    all_mcmc_results = mock_results
  )

  # Should stop at k=2 because CV < 1 disables likelihood cascade
  expect_equal(result$recommended_k, "2")
})

test_that("warnings issued for missing improvement metrics", {
  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.9)
  mock_improvements <- list(
    completed_ks = c("1", "2", "3"),
    improvements = list(
      k1_to_k2 = 4
      # k2_to_k3 missing
    )
  )

  expect_warning(
    apply_threshold_selection(
      improvements = mock_improvements,
      method = "GLOBETROTTER",
      cv = 2,
      all_mcmc_results = mock_results
    ),
    "No k2_to_k3 improvement available"
  )
})

test_that("empty completed_ks throws error", {
  mock_results <- create_mock_mcmc_results()
  mock_improvements <- list(
    completed_ks = character(0),
    improvements = list()
  )

  expect_error(
    apply_threshold_selection(
      improvements = mock_improvements,
      method = "GLOBETROTTER",
      cv = 2,
      all_mcmc_results = mock_results
    ),
    "No feasible solutions found for any k values"
  )
})

test_that("configurable CV threshold works if implemented", {
  # This test assumes you've added cv_threshold to thresholds_config
  # Skip if not yet implemented
  skip_if_not(exists("thresholds_config") &&
                "cv_threshold" %in% names(thresholds_config[["GLOBETROTTER"]]))

  mock_results <- create_mock_mcmc_results(k2_clustering_strength = 0.9)
  mock_improvements <- create_mock_improvements_to_k3()

  result <- apply_threshold_selection(
    improvements = mock_improvements,
    method = "GLOBETROTTER",
    cv = 0.8,
    all_mcmc_results = mock_results
  )

  # Behavior depends on whether cv_threshold is < or > 0.8
  expect_true(result$recommended_k %in% c("2", "3"))
})

