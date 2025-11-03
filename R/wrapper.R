#' Main wrapper function for admixplorer analysis
#'
#' @param infile Path to input file
#' @param outfile Output file prefix
#' @param method Method to use ("GLOBETROTTER" or "DATES")
#' @param ks Comma-separated string of k values to try (default "1,2,3,4")
#' @param mcmc_chains Number of MCMC chains per k (default 3)
#' @param sample_age_est Whether to estimate sample ages (default TRUE)
#' @param plot Whether to create plots (default TRUE)
#' @param apply_date_filter Whether to filter out dates > 200 gen (default TRUE)
#' @return Combined output dataframe
#' @export
admixplorer <- function(infile, outfile, method = "GLOBETROTTER",
                        ks = "1,2,3,4", mcmc_chains = 3,
                        sample_age_est = TRUE, plot = TRUE,
                        apply_date_filter = TRUE) {

  # INPUT VALIDATION

  # Check file exists
  if (!file.exists(infile)) {
    stop("Input file '", infile, "' does not exist. Please check the file path.")
  }

  output_dir <- dirname(outfile)
  if (output_dir != "." && !dir.exists(output_dir)) {
    stop("Output directory '", output_dir, "' does not exist. Please create it first or use a valid path.")
  }

  # Test if output directory is writable
  test_file <- file.path(output_dir, paste0("test_write_", Sys.getpid(), ".tmp"))
  tryCatch({
    writeLines("test", test_file)
    unlink(test_file)
  }, error = function(e) {
    stop("Cannot write to output directory '", output_dir, "'. Please check permissions.")
  })

  # Parse and validate k values
  ks_requested <- sort(unique(as.integer(unlist(strsplit(ks, ",")))))
  if (any(ks_requested <= 0)) {
    stop("All k values must be positive integers. Found invalid k values: ",
         paste(ks_requested[ks_requested <= 0], collapse = ", "))
  }

  # Validate MCMC chains
  chains <- as.numeric(mcmc_chains)
  if (chains <= 0) {
    stop("Number of MCMC chains must be a positive integer. Got: ", mcmc_chains)
  }

  # Check required packages
  required_packages <- c("ggplot2", "dplyr", "stats")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed. Please install it with: install.packages('", pkg, "')")
    }
  }

  cat("=== ADMIXPLORER ANALYSIS ===\n")
  cat("Input file:", infile, "\n")
  cat("Output prefix:", outfile, "\n")
  cat("Method:", method, "\n")
  cat("K values:", paste(ks_requested, collapse = ", "), "\n")
  cat("MCMC chains:", chains, "\n")
  cat("Sample age estimation:", sample_age_est, "\n")
  cat("Create plots:", plot, "\n\n")
  cat("Apply date filters:", apply_date_filter, "\n\n")
  # Step 1: Read RAW data first for validation
  print("READING IN DATA")
  raw_data <- utils::read.table(infile, as.is = TRUE)

  # Check for negative dates
  if (any(raw_data$V4 <= 0)) {
    negative_dates <- which(raw_data$V4 <= 0)
    stop("Negative or zero admixture dates found in individuals: ",
         paste(raw_data$V1[negative_dates], collapse = ", "),
         ". All admixture dates must be positive.")
  }

  # Check for negative standard errors
  if (any(raw_data$V5 <= 0)) {
    negative_errors <- which(raw_data$V5 <= 0)
    stop("Negative or zero standard errors found in individuals: ",
         paste(raw_data$V1[negative_errors], collapse = ", "),
         ". All standard errors must be positive.")
  }

  # Check for invalid age ranges (V3 < V2)
  if (any(raw_data$V3 < raw_data$V2)) {
    invalid_ages <- which(raw_data$V3 < raw_data$V2)
    stop("Invalid age ranges found (max age < min age) in individuals: ",
         paste(raw_data$V1[invalid_ages], collapse = ", "),
         ". Maximum sampling age must be >= minimum sampling age.")
  }

  # Check for negative sampling ages
  if (any(raw_data$V2 < 0) || any(raw_data$V3 < 0)) {
    negative_sampling <- which(raw_data$V2 < 0 | raw_data$V3 < 0)
    stop("Negative sampling ages found in individuals: ",
         paste(raw_data$V1[negative_sampling], collapse = ", "),
         ". All sampling ages must be non-negative.")
  }

  filtered_data <- read_and_filter_data(infile, apply_date_filter = apply_date_filter)

  # Check if any data remains after filtering
  if (nrow(filtered_data$data) == 0) {
    stop("No data remaining after filtering. All individuals were removed as outliers. ",
         "Please check your data quality or filtering criteria.")
  }

  # Step 2: Prepare data for clustering
  prepared_data <- prepare_clustering_data(filtered_data$data, sample_age_est)

  # Calculate CV for scaling
  cv <- mean(filtered_data$data$V4 / filtered_data$data$V5)
  cat("Coefficient of variation:", round(cv, 3), "\n")

  # Step 3: Run clustering analysis
  all_mcmc_results <- run_clustering_analysis(
    prepared_data = prepared_data,
    ks = ks_requested,
    chains = chains,
    sample_age_est = sample_age_est,
    plot = plot,
    outfile = outfile
  )

  completed_ks <- names(all_mcmc_results)[!sapply(all_mcmc_results, is.null)]

  if (length(completed_ks) == 0) {
    # ALL k values failed - stop completely
    warning(
      "No feasible solutions found for any k values (",
      paste(ks_requested, collapse = ", "),
      "). Stopping analysis. Consider different k values or checking data quality."
    )
    return(invisible(NULL))
  } else if (length(completed_ks) < length(ks_requested)) {
    # Some k values failed, but at least one succeeded - continue with a warning
    failed_ks <- setdiff(as.character(ks_requested), completed_ks)
    cat("Warning: No feasible solutions for k =", paste(failed_ks, collapse = ", "),
        ". Continuing with k =", paste(completed_ks, collapse = ", "), "\n")
  }

  # Step 4: Calculate likelihood improvements
  likelihood_analysis <- calculate_likelihood_improvements(
    all_mcmc_results = all_mcmc_results,
    n_individuals = as.numeric(nrow(filtered_data$data))
  )

  # Step 5: Apply threshold selection
  threshold_results <- apply_threshold_selection(
    improvements = likelihood_analysis,
    method = method,
    cv = cv,
    all_mcmc_results = all_mcmc_results
  )

  # Step 6: Generate final output
  final_output <- generate_final_output(
    all_mcmc_results = all_mcmc_results,
    recommended_k = threshold_results$recommended_k,
    original_data = filtered_data$data,
    pop.vec = prepared_data$pop.vec,
    dates_original = prepared_data$dates_original,
    sample_age_est = sample_age_est,
    outfile = outfile,
    method = method,
    thresholds_info = threshold_results
  )

  return(invisible(final_output))
}
