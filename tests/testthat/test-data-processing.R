# tests/testthat/test-data-processing.R (FINAL CORRECTED VERSION)
# Test the data input and filtering functions

test_that("read_and_filter_data works with valid input", {
  # Create temporary test file
  test_data <- create_test_data_1cluster()
  temp_file <- tempfile(fileext = ".txt")
  write.table(test_data, temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)

  result <- read_and_filter_data(temp_file)

  expect_type(result, "list")
  expect_named(result, c("data", "original_inds", "removed_inds"))
  expect_s3_class(result$data, "data.frame")
  expect_equal(nrow(result$data), 10)  # All rows should pass filters
  expect_length(result$removed_inds, 0)  # No outliers removed

  unlink(temp_file)
})

test_that("read_and_filter_data removes outliers correctly", {
  # Test the exact data from your debug script
  outlier_data <- data.frame(
    V1 = paste0("ind", 1:6),
    V2 = rep(0, 6),
    V3 = rep(0, 6),
    V4 = c(45, 1, 50, 1, 40, 35),  # ind2: V4=1 (fails V4>2), ind4: V4=200 (fails V4<150)
    V5 = c(3, 2, 4, 3, 50, 2),       # ind5: V5=50, V4=40, 50<400 so passes
    stringsAsFactors = FALSE
  )

  temp_file <- tempfile(fileext = ".txt")
  write.table(outlier_data, temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)

  result <- read_and_filter_data(temp_file)

  # Based on your actual filters:
  # ind2 removed: V4=1 fails V4 > 2
  # ind4 removed: V4=200 fails V4 < 150
  # ind5 kept: V5=50, V4=40, 50 < 40*10=400 (passes)
  expect_equal(nrow(result$data), 4)  # 6 - 2 = 4 remaining
  expect_length(result$removed_inds, 2)  # Only 2 outliers removed
  expect_setequal(result$removed_inds, c("ind2", "ind4"))  # Exactly these two
  expect_true("ind5" %in% result$data$V1)  # ind5 should remain

  unlink(temp_file)
})

test_that("read_and_filter_data filters work individually", {
  # Test V4 > 2 filter
  low_date_data <- data.frame(
    V1 = c("ind1", "ind2"),
    V2 = c(0, 0),
    V3 = c(0, 0),
    V4 = c(10, 1),  # ind2 has V4=1 < 2
    V5 = c(2, 1),
    stringsAsFactors = FALSE
  )

  temp_file <- tempfile(fileext = ".txt")
  write.table(low_date_data, temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)

  result <- read_and_filter_data(temp_file)
  expect_equal(nrow(result$data), 1)  # Only ind1 remains
  expect_equal(result$removed_inds, "ind2")

  unlink(temp_file)

  # Test V4 < 150 filter
  high_date_data <- data.frame(
    V1 = c("ind1", "ind2"),
    V2 = c(0, 0),
    V3 = c(0, 0),
    V4 = c(100, 1),  # ind2 has V4=200 > 150
    V5 = c(5, 5),
    stringsAsFactors = FALSE
  )

  temp_file <- tempfile(fileext = ".txt")
  write.table(high_date_data, temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)

  result <- read_and_filter_data(temp_file)
  expect_equal(nrow(result$data), 1)  # Only ind1 remains
  expect_equal(result$removed_inds, "ind2")

  unlink(temp_file)

  # Test V5 < V4 * 10 filter
  high_error_data <- data.frame(
    V1 = c("ind1", "ind2"),
    V2 = c(0, 0),
    V3 = c(0, 0),
    V4 = c(20, 5),     # V4 * 10 = 200, 50 respectively
    V5 = c(10, 60),    # ind2: V5=60 > V4*10=50
    stringsAsFactors = FALSE
  )

  temp_file <- tempfile(fileext = ".txt")
  write.table(high_error_data, temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)

  result <- read_and_filter_data(temp_file)
  expect_equal(nrow(result$data), 1)  # Only ind1 remains
  expect_equal(result$removed_inds, "ind2")

  unlink(temp_file)
})

test_that("prepare_clustering_data works correctly", {
  test_data <- create_test_data_1cluster()

  result <- prepare_clustering_data(test_data, sample_age_est = TRUE)

  expect_type(result, "list")
  expect_true(all(c("pop.vec", "age.matrix", "dates_original", "std_errors") %in% names(result)))
  expect_length(result$pop.vec, 10)
  expect_equal(nrow(result$age.matrix), 10)
  expect_length(result$dates_original, 10)
  expect_length(result$std_errors, 10)
})

test_that("prepare_clustering_data handles sample_age_est = FALSE", {
  test_data <- create_test_data_1cluster()

  result <- prepare_clustering_data(test_data, sample_age_est = FALSE)

  # When sample_age_est = FALSE, age matrix should be set to midpoint
  expect_true(all(result$age.matrix[,1] == result$age.matrix[,2]))
})
