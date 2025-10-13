# tests/testthat/test-admixplorer-main.R
# End-to-end tests for the main admixplorer() function

test_that("admixplorer works end-to-end with 1-cluster data", {
  test_data <- create_test_data_1cluster()
  temp_file <- tempfile(fileext = ".txt")
  temp_dir <- tempdir()
  output_prefix <- file.path(temp_dir, "test_1cluster")
  
  write.table(test_data, temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  result <- admixplorer(
    infile = temp_file,
    outfile = output_prefix,
    method = "GLOBETROTTER",
    ks = "1,2",  # Just test k=1,2 for speed
    mcmc_chains = 1  # Reduce for testing speed
  )
  
  # Check basic output structure
  expect_s3_class(result, "data.frame")
  expect_true("model_k" %in% colnames(result))
  expect_true("recommended_k" %in% colnames(result))
  expect_true("pop" %in% colnames(result))
  
  # Should recommend k=1 for 1-cluster data
  expect_equal(unique(result$recommended_k), "1")
  
  # Check output file exists
  expect_true(file.exists(paste0(output_prefix, ".output.txt")))
  
  # Clean up
  unlink(temp_file)
  unlink(paste0(output_prefix, ".output.txt"))
})

test_that("admixplorer works end-to-end with 2-cluster data", {
  test_data <- create_test_data_2cluster()
  temp_file <- tempfile(fileext = ".txt") 
  temp_dir <- tempdir()
  output_prefix <- file.path(temp_dir, "test_2cluster")
  
  write.table(test_data, temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  result <- admixplorer(
    infile = temp_file,
    outfile = output_prefix, 
    method = "GLOBETROTTER",
    ks = "1,2",
    mcmc_chains = 1
  )
  
  # Should recommend k>=2 for 2-cluster data (hopefully!)
  recommended <- unique(result$recommended_k)
  expect_true(recommended %in% c("2", "3", "4+"))
  
  unlink(temp_file)
  unlink(paste0(output_prefix, ".output.txt"))
})

test_that("admixplorer handles low CV data correctly", {
  test_data <- create_low_cv_data()  # cv < 1
  temp_file <- tempfile(fileext = ".txt")
  temp_dir <- tempdir() 
  output_prefix <- file.path(temp_dir, "test_lowcv")
  
  write.table(test_data, temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Should not error and should use only clustering strength
  expect_no_error(
    result <- admixplorer(
      infile = temp_file,
      outfile = output_prefix,
      ks = "1,2",
      mcmc_chains = 1
    )
  )
  
  unlink(temp_file)
  unlink(paste0(output_prefix, ".output.txt"))
})

test_that("admixplorer respects method parameter", {
  test_data <- create_test_data_1cluster()
  temp_file <- tempfile(fileext = ".txt")
  temp_dir <- tempdir()
  output_prefix <- file.path(temp_dir, "test_dates_method")
  
  write.table(test_data, temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  result <- admixplorer(
    infile = temp_file,
    outfile = output_prefix,
    method = "DATES",  # Different method
    ks = "1,2",
    mcmc_chains = 1
  )
  
  # Should complete without error
  expect_s3_class(result, "data.frame")
  
  # Check output file mentions DATES method
  output_content <- readLines(paste0(output_prefix, ".output.txt"))
  expect_true(any(grepl("DATES", output_content)))
  
  unlink(temp_file)
  unlink(paste0(output_prefix, ".output.txt"))
})