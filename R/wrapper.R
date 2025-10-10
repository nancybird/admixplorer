#' Main wrapper function for admixplorer analysis
#'
#' @param infile Path to input file
#' @param outfile Output file prefix
#' @param method Method to use ("GLOBETROTTER" or "DATES")
#' @param ks Comma-separated string of k values to try (default "1,2,3,4")
#' @param mcmc_chains Number of MCMC chains per k (default 3)
#' @param sample_age_est Whether to estimate sample ages (default TRUE)
#' @param plot Whether to create plots (default TRUE)
#' @return Combined output dataframe
#' @export
admixplorer <- function(infile, outfile, method = "GLOBETROTTER",
                               ks = "1,2,3,4", mcmc_chains = 3,
                               sample_age_est = TRUE, plot = TRUE) {

  # Check required packages
  required_packages <- c("ggplot2", "dplyr", "stats")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {

    }
  }

  # Parse k values
  ks_requested <- sort(unique(as.integer(unlist(strsplit(ks, ",")))))
  chains <- as.numeric(mcmc_chains)

  cat("=== ADMIXCRAFT ANALYSIS ===\n")
  cat("Input file:", infile, "\n")
  cat("Output prefix:", outfile, "\n")
  cat("Method:", method, "\n")
  cat("K values:", paste(ks_requested, collapse = ", "), "\n")
  cat("MCMC chains:", chains, "\n")
  cat("Sample age estimation:", sample_age_est, "\n")
  cat("Create plots:", plot, "\n\n")

  # Step 1: Read and filter data
  print("READING IN DATA")
  filtered_data <- read_and_filter_data(infile)

  # Step 2: Prepare data for clustering
  prepared_data <- prepare_clustering_data(filtered_data$data, sample_age_est)

  # Calculate CV for scaling
  cv <- mean(filtered_data$data$V4 / filtered_data$data$V5)

  # Step 3: Run clustering analysis
  all_mcmc_results <- run_clustering_analysis(
    prepared_data = prepared_data,
    ks = ks_requested,
    chains = chains,
    sample_age_est = sample_age_est,
    plot = plot,
    outfile = outfile
  )

  # Step 4: Calculate likelihood improvements
  likelihood_analysis <- calculate_likelihood_improvements(
    all_mcmc_results = all_mcmc_results,
    n_individuals = nrow(filtered_data$data)
  )

threshold_results <- apply_threshold_selection(
  improvements = likelihood_analysis,
  method = method,
  cv = cv,
  all_mcmc_results = all_mcmc_results  # ADD THIS LINE
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

  return(final_output)
}


#' Check and install required packages
#'
#' @param pkg Package name
#' @keywords internal
check_install_packages <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    warning(paste("Package", pkg, "is not installed. Please install it before running"))
    quit(status = 1)
  }
}
