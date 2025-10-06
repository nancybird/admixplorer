# Cluster Statistics Helpers

# Function to compute weighted cluster means
compute_cluster_means <- function(assignments, dates, std_errors, k) {
  means <- numeric(k)
  for (c in 1:k) {
    inds <- which(assignments == c)
    if (length(inds) > 0) {
      weights <- 1 / (std_errors[inds]^2)
      means[c] <- sum(dates[inds] * weights) / sum(weights)
    } else {
      means[c] <- NA  # Avoid division by zero
    }
  }
  return(means)
}
compute_cluster_sds <- function(assignments, std_errors, n, k) {
  sds <- numeric(k)
  for (c in 1:k) {
    inds <- which(assignments == c)
    if (length(inds) > 0) {
      weights <- 1 / (std_errors[inds]^2)
      sds[c] <- sqrt(1 / sum(weights))
    } else {
      # Default value for empty clusters
      sds[c] <- mean(std_errors)
    }
  }
  return(sds)
}


