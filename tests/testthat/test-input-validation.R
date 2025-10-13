# Now update your tests to expect proper errors:

# tests/testthat/test-input-validation.R (PROPER VERSION)

test_that("admixplorer handles missing file", {
  expect_error(
    admixplorer("nonexistent_file.txt", "output"),
    "does not exist"
  )
})

test_that("admixplorer validates k values", {
  test_data <- create_test_data_1cluster()
  temp_file <- tempfile(fileext = ".txt")
  write.table(test_data, temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  expect_error(
    admixplorer(temp_file, "output", ks = "0,1,2"),
    "All k values must be positive integers"
  )
  
  unlink(temp_file)
})

test_that("admixplorer validates mcmc_chains", {
  test_data <- create_test_data_1cluster()
  temp_file <- tempfile(fileext = ".txt")
  write.table(test_data, temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  expect_error(
    admixplorer(temp_file, "output", mcmc_chains = 0),
    "Number of MCMC chains must be a positive integer"
  )
  
  unlink(temp_file)
})

test_that("admixplorer validates negative dates", {
  bad_data <- data.frame(
    V1 = paste0("ind", 1:3),
    V2 = rep(0, 3),
    V3 = rep(0, 3), 
    V4 = c(45, -10, 50),  # Negative date
    V5 = c(3, 2, 4),
    stringsAsFactors = FALSE
  )
  
  temp_file <- tempfile(fileext = ".txt")
  write.table(bad_data, temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  expect_error(
    admixplorer(temp_file, "output"),
    "Negative or zero admixture dates found"
  )
  
  unlink(temp_file)
})

test_that("admixplorer validates negative standard errors", {
  bad_data <- data.frame(
    V1 = paste0("ind", 1:3),
    V2 = rep(0, 3),
    V3 = rep(0, 3), 
    V4 = c(45, 40, 50),
    V5 = c(3, -2, 4),     # Negative error
    stringsAsFactors = FALSE
  )
  
  temp_file <- tempfile(fileext = ".txt")
  write.table(bad_data, temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  expect_error(
    admixplorer(temp_file, "output"),
    "Negative or zero standard errors found"
  )
  
  unlink(temp_file)
})

test_that("admixplorer validates age ranges", {
  bad_data <- data.frame(
    V1 = paste0("ind", 1:3),
    V2 = c(0, 100, 0),    # ind2: min age = 100
    V3 = c(0, 50, 0),     # ind2: max age = 50 (max < min!)
    V4 = c(45, 40, 50),
    V5 = c(3, 2, 4),
    stringsAsFactors = FALSE
  )
  
  temp_file <- tempfile(fileext = ".txt")
  write.table(bad_data, temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  expect_error(
    admixplorer(temp_file, "output"),
    "Invalid age ranges found"
  )
  
  unlink(temp_file)
})

test_that("admixplorer validates negative sampling ages", {
  bad_data <- data.frame(
    V1 = paste0("ind", 1:3),
    V2 = c(0, -100, 0),   # Negative min age
    V3 = c(0, 0, 0),
    V4 = c(45, 40, 50),
    V5 = c(3, 2, 4),
    stringsAsFactors = FALSE
  )
  
  temp_file <- tempfile(fileext = ".txt")
  write.table(bad_data, temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  expect_error(
    admixplorer(temp_file, "output"),
    "Negative sampling ages found"
  )
  
  unlink(temp_file)
})

test_that("admixplorer handles all data filtered out", {
  # Create data where all will be filtered during outlier removal
  all_bad_data <- data.frame(
    V1 = paste0("ind", 1:3),
    V2 = rep(0, 3),
    V3 = rep(0, 3),
    V4 = c(1, 200, 1),     # All fail filters
    V5 = c(1, 1, 1),
    stringsAsFactors = FALSE
  )
  
  temp_file <- tempfile(fileext = ".txt")
  write.table(all_bad_data, temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  expect_error(
    admixplorer(temp_file, "output"),
    "No data remaining after filtering"
  )
  
  unlink(temp_file)
})

test_that("admixplorer works with valid data", {
  test_data <- create_test_data_1cluster()
  temp_file <- tempfile(fileext = ".txt")
  temp_dir <- tempdir()
  output_prefix <- file.path(temp_dir, "test_output")
  
  write.table(test_data, temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Should work without errors
  result <- admixplorer(temp_file, output_prefix, ks = "1,2", mcmc_chains = 1)
  
  expect_s3_class(result, "data.frame")
  expect_true("recommended_k" %in% colnames(result))
  expect_true(file.exists(paste0(output_prefix, ".output.txt")))
  
  # Clean up
  unlink(temp_file)
  unlink(paste0(output_prefix, ".output.txt"))
})