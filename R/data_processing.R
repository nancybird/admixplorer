#' Read and filter admixture data
#'
#' @param filename Path to input file
#' @return List containing filtered dataframe and original individuals
#' @export
#' @importFrom utils read.table
#' @importFrom dplyr filter

#' Read and filter admixture data
#'
#' @param filename Path to input file
#' @param apply_date_filter Logical, whether to apply date-based outlier filters (default TRUE)
#' @return List containing filtered dataframe and original individuals
#' @export
#' @importFrom utils read.table
#' @importFrom dplyr filter

read_and_filter_data <- function(filename, apply_date_filter = TRUE) {
  # Read the file
  file.read <- utils::read.table(filename, as.is = TRUE)
  original_inds <- file.read$V1

  # Apply filters
  if (apply_date_filter) {
    cat("Applying date-based outlier filters...\n")
    file.read <- file.read[file.read$V4 > 2, ]
    file.read <- file.read[file.read$V5 < file.read$V4 * 10, ]
    file.read <- file.read[file.read$V4 < 200, ]
  } else {
    cat("Skipping date-based outlier filters (apply_date_filter = FALSE)\n")
  }

  # Find which individuals were removed
  remaining_inds <- file.read$V1
  removed_inds <- setdiff(original_inds, remaining_inds)

  # Print outliers removed with their names
  cat("Outliers removed:", length(removed_inds), "\n")
  if(length(removed_inds) > 0) {
    cat("Removed individuals:", paste(removed_inds, collapse = ", "), "\n")
  }

  list(
    data = file.read,
    original_inds = original_inds,
    removed_inds = removed_inds
  )
}

#' Prepare clustering data
#'
#' @param data Filtered dataframe from read_and_filter_data
#' @param sample_age_est Whether to estimate sample age or use mean
#' @return List of prepared data for clustering
#' @export
prepare_clustering_data <- function(data, sample_age_est = TRUE) {
  pop.vec <- as.character(data[,1])
  age.matrix <- matrix(as.matrix(as.double(unlist(data[,2:3]))), ncol = 2)
  clusters <- rep(1, nrow(data))

  if (sample_age_est == FALSE) {
    age.matrix[,1:2] <- (age.matrix[,1] + age.matrix[,2])/2
  }

  age.matrix.individual <- age.matrix
  age.matrix.individual[,1:2] <- (age.matrix.individual[,1] + age.matrix.individual[,2])/2
  dates_original <- as.numeric(data[,4])
  dates <- dates_original + age.matrix.individual[,1]
  std_errors <- as.numeric(data[,5])

  clusters <- rep(1, length(dates))
  midpoint_ages <- age.matrix.individual[,1]

  # Calculate initial weights and joint date
  weights <- 1 / std_errors^2
  joint_date <- sum(dates * weights) / sum(weights)

  list(
    pop.vec = pop.vec,
    age.matrix = age.matrix,
    age.matrix.individual = age.matrix.individual,
    dates_original = dates_original,
    dates = dates,
    std_errors = std_errors,
    clusters = clusters,
    midpoint_ages = midpoint_ages,
    joint_date = joint_date
  )
}
