# Validate the complete state that will be passed to MCMC.
# Preserves the existing sampling-age LOWER-bound constraint.
valid_mcmc_start <- function(
    assignments, means, sampling_ages, age_ranges, k
) {
  n <- length(sampling_ages)

  if (
    n < k ||
    length(assignments) != n ||
    length(means) != k ||
    any(!is.finite(assignments)) ||
    any(assignments != floor(assignments)) ||
    any(assignments < 1L | assignments > k) ||
    any(!is.finite(means)) ||
    any(!is.finite(sampling_ages))
  ) {
    return(FALSE)
  }

  # Every requested cluster must contain at least one individual.
  if (any(tabulate(assignments, nbins = k) == 0L)) {
    return(FALSE)
  }

  if (!is.null(age_ranges)) {
    if (
      !is.matrix(age_ranges) ||
      !all(dim(age_ranges) == c(n, 2L)) ||
      any(!is.finite(age_ranges))
    ) {
      return(FALSE)
    }

    # Starting ages must lie within their allowed ranges.
    if (any(
      sampling_ages < age_ranges[, 1] |
      sampling_ages > age_ranges[, 2]
    )) {
      return(FALSE)
    }

    # An assigned admixture date cannot be more recent than
    # the individual's sampling-age lower bound.
    if (any(means[assignments] < age_ranges[, 1])) {
      return(FALSE)
    }
  }

  TRUE
}
