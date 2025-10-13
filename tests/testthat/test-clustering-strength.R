# tests/testthat/test-clustering-strength.R
# Test the clustering strength calculation

test_that("clustering_strength calculates correctly", {
  # Perfect clustering: all 1s and 0s
  perfect_matrix <- matrix(c(
    1, 1, 0, 0,
    1, 1, 0, 0, 
    0, 0, 1, 1,
    0, 0, 1, 1
  ), nrow = 4, byrow = TRUE)
  
  strength <- clustering_strength(perfect_matrix)
  expect_equal(strength, 1.0)  # Perfect certainty
})

test_that("clustering_strength handles uncertain clustering", {
  # All 0.5s = completely uncertain
  uncertain_matrix <- matrix(0.5, nrow = 4, ncol = 4)
  diag(uncertain_matrix) <- 1  # Diagonal should be 1
  
  strength <- clustering_strength(uncertain_matrix)
  expect_equal(strength, 0.0)  # No certainty
})

test_that("clustering_strength handles mixed clustering", {
  # Mix of certain and uncertain
  mixed_matrix <- matrix(c(
    1, 0.8, 0.2, 0.1,
    0.8, 1, 0.3, 0.2,
    0.2, 0.3, 1, 0.9, 
    0.1, 0.2, 0.9, 1
  ), nrow = 4, byrow = TRUE)
  
  strength <- clustering_strength(mixed_matrix)
  expect_true(strength > 0 && strength < 1)  # Somewhere in between
  expect_type(strength, "double")
})

test_that("clustering_strength returns valid range", {
  # Generate random co-clustering matrix
  set.seed(123)
  n <- 6
  random_matrix <- matrix(runif(n*n, 0, 1), nrow = n)
  diag(random_matrix) <- 1
  # Make symmetric
  random_matrix[lower.tri(random_matrix)] <- t(random_matrix)[lower.tri(random_matrix)]
  
  strength <- clustering_strength(random_matrix)
  expect_true(strength >= 0 && strength <= 1)
  expect_type(strength, "double")
  expect_length(strength, 1)
})