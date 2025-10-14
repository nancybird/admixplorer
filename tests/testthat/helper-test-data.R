# Start with these test files - run usethis::use_testthat(3) first

# 1. TEST DATA SETUP
# Create this file: tests/testthat/helper-test-data.R
# (Files starting with "helper-" are loaded before tests run)

create_test_data_1cluster <- function() {
  # Simple 1-cluster data that should return k=1
  data.frame(
    V1 = paste0("ind", 1:10),  # Individual IDs
    V2 = rep(0, 10),           # Min age (modern samples)
    V3 = rep(0, 10),           # Max age (modern samples)
    V4 = rnorm(10, 45, 2),     # Admixture dates around 45 generations
    V5 = rnorm(10, 10, 0.2),    # Standard errors around 3
    stringsAsFactors = FALSE
  )
}

create_test_data_2cluster <- function() {
  # 2-cluster data: half at ~30 gen, half at ~60 gen
  data.frame(
    V1 = paste0("ind", 1:10),
    V2 = rep(0, 10),
    V3 = rep(0, 10),
    V4 = c(rep(30, 5), rep(60, 5)),  # Two distinct clusters
    V5 = rep(2, 10),                 # Low error for clear separation
    stringsAsFactors = FALSE
  )
}

create_bad_data_wrong_columns <- function() {
  # Missing columns
  data.frame(
    V1 = paste0("ind", 1:5),
    V2 = rep(0, 5),
    V3 = rep(0, 5),
    V4 = rep(45, 5)
    # Missing V5 column
  )
}

create_bad_data_negatives <- function() {
  # Negative dates/errors
  data.frame(
    V1 = paste0("ind", 1:5),
    V2 = rep(0, 5),
    V3 = rep(0, 5),
    V4 = c(45, -10, 50, 40, 35),  # Negative date
    V5 = c(3, 2, -1, 4, 2),       # Negative error
    stringsAsFactors = FALSE
  )
}

create_time_traveler_data <- function() {
  # One time traveler (ind2): date=1, sampling age=1000, small error
  # Others normal: date ~ 50, sampling age=0
  data.frame(
    V1 = paste0("ind", 1:5),
    V2 = c(0, 1000, 0, 0, 0),      # Min ages (sampling age)
    V3 = c(0, 1000, 0, 0, 0),      # Max ages = same
    V4 = c(45, 1, 50, 55, 60),     # ind2 has date=1
    V5 = c(3, 1, 3, 3, 3),         # ind2 has small error
    stringsAsFactors = FALSE
  )
}

# In tests/testthat/helper-test-data.R

create_feasible_ancient_data <- function() {
  # Two ancient samples (ind2, ind4) with moderate ages
  # and admixture dates that still allow feasible solutions
  data.frame(
    V1 = paste0("ind", 1:6),
    V2 = c(0, 50, 0, 80, 0, 0),     # ind2:50, ind4:80 sampling ages
    V3 = c(0, 50, 0, 80, 0, 0),     # same max ages
    V4 = c(45, 20, 50, 60, 55, 40), # dates > sampling ages for ind2/4
    V5 = c(3, 2, 3, 4, 3, 3),       # low errors for convergence
    stringsAsFactors = FALSE
  )
}



create_low_cv_data <- function() {
  # Data with cv < 1 (high errors relative to dates)
  data.frame(
    V1 = paste0("ind", 1:10),
    V2 = rep(0, 10),
    V3 = rep(0, 10),
    V4 = rep(10, 10),         # Low dates
    V5 = rep(15, 10),         # High errors (cv = 10/15 = 0.67)
    stringsAsFactors = FALSE
  )
}
