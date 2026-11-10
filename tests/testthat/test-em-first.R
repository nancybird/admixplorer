# Save as tests/testthat/test-em-first.R (replace the old standalone script).
# Uses the loaded package and the existing helper-test-data.R fixture.
# EM and MCMC are stubbed to test initialisation without running long chains.

.em_first_harness <- function(em_assignments, chains = 1L,
                              feasible = TRUE, variable_ages = FALSE,
                              bad_fallback = FALSE) {
  prepared <- getFromNamespace("prepare_clustering_data", "admixplorer")(
    create_feasible_ancient_data()
  )
  if (variable_ages) {
    prepared$age.matrix[2L, 2L] <- 60
    prepared$age.matrix[4L, 2L] <- 100
    prepared$midpoint_ages <- rowMeans(prepared$age.matrix)
  }

  process <- getFromNamespace("process_single_k", "admixplorer")
  validate <- getFromNamespace("valid_mcmc_start", "admixplorer")
  means_fn <- getFromNamespace("compute_cluster_means", "admixplorer")
  sandbox <- new.env(parent = environment(process))
  environment(process) <- sandbox

  seen <- new.env(parent = emptyenv())
  seen$events <- character()
  seen$starts <- list()
  seen$feasibility_calls <- 0L
  record <- function(event) seen$events <- c(seen$events, event)

  valid_assignments <- c(1L, 1L, 1L, 2L, 1L, 1L)
  fallback_ages <- prepared$age.matrix[, 1L]
  if (variable_ages) {
    fallback_ages[2L] <- 53
    fallback_ages[4L] <- 87
  }
  fallback_assignments <- if (bad_fallback) {
    c(1L, 1L, 1L, 1L, 2L, 2L)
  } else {
    valid_assignments
  }
  solution <- list(assignments = fallback_assignments,
                   sampling_ages = fallback_ages)
  # Two entries also avoid R's sample(single_number, 1) special case.
  pool <- if (feasible) list(solution, solution) else list()

  sandbox$em_clustering_original <- function(dates, std_errors, k, ...) {
    record("EM")
    responsibilities <- matrix(0, nrow = length(dates), ncol = k)
    responsibilities[cbind(seq_along(em_assignments), em_assignments)] <- 1
    list(
      responsibilities = responsibilities,
      cluster_means = means_fn(em_assignments, dates, std_errors, k)
    )
  }

  sandbox$valid_mcmc_start <- function(...) {
    record("validate")
    validate(...)
  }

  sandbox$mc_feasibility_with_solutions <- function(
    dates, std_errors, k, age_ranges, n_draws = 200, ...) {
    record("feasibility")
    seen$feasibility_calls <- seen$feasibility_calls + 1L
    list(
      n_draws = n_draws,
      n_feasible = length(pool),
      prop_feasible = length(pool) / n_draws,
      n_draws_with_infeasible_assignments = n_draws,
      prop_draws_with_infeasible = 1,
      total_infeasible_assignments_tested = n_draws,
      infeasible_counts_per_draw = rep(1L, n_draws),
      solutions = pool
    )
  }

  sandbox$mcmc_clustering <- function(
    dates, std_errors, init_cluster_means, cluster_assignments,
    k, sampling_ages, sample_ages, age_ranges, ...) {
    record("MCMC")
    start <- list(assignments = cluster_assignments,
                  means = init_cluster_means,
                  sampling_ages = sampling_ages)
    seen$starts[[length(seen$starts) + 1L]] <- start
    if (!isTRUE(validate(cluster_assignments, init_cluster_means,
                         sampling_ages, age_ranges, k))) {
      stop("TEST: an invalid starting state reached MCMC.")
    }
    list(result = data.frame(final_log_likelihood_best_sampling_ages = -1))
  }

  list(
    run = function(k = 2L) {
      result <- NULL
      invisible(capture.output(
        result <- process(k_current = k, prepared_data = prepared,
                          chains = chains, sample_age_est = TRUE,
                          plot = FALSE, outfile = tempfile())
      ))
      result
    },
    seen = seen,
    prepared = prepared,
    solution = solution,
    means_fn = means_fn
  )
}

test_that("the start validator rejects time travellers and incomplete states", {
  p <- getFromNamespace("prepare_clustering_data", "admixplorer")(
    create_feasible_ancient_data()
  )
  validate <- getFromNamespace("valid_mcmc_start", "admixplorer")
  means_fn <- getFromNamespace("compute_cluster_means", "admixplorer")
  ages <- p$midpoint_ages
  good <- c(1L, 1L, 1L, 2L, 1L, 1L)
  bad <- c(1L, 1L, 1L, 1L, 2L, 2L)
  good_means <- means_fn(good, p$dates_original + ages, p$std_errors, 2L)
  bad_means <- means_fn(bad, p$dates_original + ages, p$std_errors, 2L)

  expect_true(validate(good, good_means, ages, p$age.matrix, 2L))
  expect_false(validate(bad, bad_means, ages, p$age.matrix, 2L))
  expect_false(validate(rep(1L, 6L), c(good_means[1L], NA_real_),
                        ages, p$age.matrix, 2L))
  expect_false(validate(good, c(NA_real_, good_means[2L]),
                        ages, p$age.matrix, 2L))
  outside <- ages
  outside[2L] <- outside[2L] + 1
  expect_false(validate(good, good_means, outside, p$age.matrix, 2L))
  expect_false(validate(good[-1L], good_means, ages, p$age.matrix, 2L))
})

test_that("a valid EM start bypasses feasibility testing", {
  h <- .em_first_harness(c(1L, 1L, 1L, 2L, 1L, 1L), chains = 2L)
  expect_no_error(h$run())
  expect_identical(h$seen$events[1L], "EM")
  expect_identical(sum(h$seen$events == "EM"), 1L)
  expect_identical(h$seen$feasibility_calls, 0L)
  expect_length(h$seen$starts, 2L)
  expect_true(match("validate", h$seen$events) < match("MCMC", h$seen$events))
  for (start in h$seen$starts) {
    expect_equal(start$assignments, c(1L, 1L, 1L, 2L, 1L, 1L))
  }
})

test_that("invalid EM starts trigger one cached feasibility search before MCMC", {
  h <- .em_first_harness(c(1L, 1L, 1L, 1L, 2L, 2L), chains = 3L,
                         variable_ages = TRUE)
  expect_no_error(h$run())
  events <- h$seen$events
  expect_identical(events[1L], "EM")
  expect_true(match("validate", events) < match("feasibility", events))
  expect_true(match("feasibility", events) < match("MCMC", events))
  expect_identical(h$seen$feasibility_calls, 1L)
  expect_length(h$seen$starts, 3L)
  for (start in h$seen$starts) {
    expect_equal(start$assignments, h$solution$assignments)
    expect_equal(start$sampling_ages, h$solution$sampling_ages)
    expect_equal(start$means, h$means_fn(
      h$solution$assignments,
      h$prepared$dates_original + h$solution$sampling_ages,
      h$prepared$std_errors, 2L
    ))
  }
})

test_that("an empty EM cluster triggers feasibility fallback", {
  h <- .em_first_harness(rep(1L, 6L))
  expect_no_error(h$run())
  expect_identical(h$seen$events[1L], "EM")
  expect_identical(h$seen$feasibility_calls, 1L)
  expect_length(h$seen$starts, 1L)
})

test_that("k is skipped without MCMC when no feasible solution is found", {
  h <- .em_first_harness(c(1L, 1L, 1L, 1L, 2L, 2L), feasible = FALSE)
  expect_null(h$run())
  expect_identical(h$seen$events[1L], "EM")
  expect_identical(h$seen$feasibility_calls, 1L)
  expect_length(h$seen$starts, 0L)
})

test_that("an invalid stored fallback never reaches MCMC", {
  h <- .em_first_harness(c(1L, 1L, 1L, 1L, 2L, 2L), bad_fallback = TRUE)
  # Either skipping k or explicitly stopping is a safe failure policy.
  result <- tryCatch(h$run(), error = identity)
  expect_true(is.null(result) || inherits(result, "error"))
  expect_identical(h$seen$feasibility_calls, 1L)
  expect_length(h$seen$starts, 0L)
})

test_that("EM-derived means use each chain's actual starting sampling ages", {
  h <- .em_first_harness(c(1L, 1L, 1L, 2L, 1L, 1L), chains = 2L,
                         variable_ages = TRUE)
  expect_no_error(h$run())
  expect_identical(h$seen$feasibility_calls, 0L)
  expect_length(h$seen$starts, 2L)
  for (start in h$seen$starts) {
    expect_equal(start$means, h$means_fn(
      start$assignments,
      h$prepared$dates_original + start$sampling_ages,
      h$prepared$std_errors, 2L
    ))
    expect_true(all(start$sampling_ages >= h$prepared$age.matrix[, 1L]))
    expect_true(all(start$sampling_ages <= h$prepared$age.matrix[, 2L]))
  }
})

test_that("n smaller than k is skipped before EM or feasibility testing", {
  h <- .em_first_harness(c(1L, 1L, 1L, 2L, 1L, 1L))
  expect_null(h$run(k = 7L))
  expect_length(h$seen$events, 0L)
  expect_length(h$seen$starts, 0L)
})

test_that("real MCMC rejects an invalid state at entry", {
  p <- getFromNamespace("prepare_clustering_data", "admixplorer")(
    create_feasible_ancient_data()
  )
  means_fn <- getFromNamespace("compute_cluster_means", "admixplorer")
  mcmc <- getFromNamespace("mcmc_clustering", "admixplorer")
  a <- c(1L, 1L, 1L, 1L, 2L, 2L)
  m <- means_fn(a, p$dates_original + p$midpoint_ages, p$std_errors, 2L)
  expect_error(
    mcmc(dates = p$dates_original, std_errors = p$std_errors,
         init_cluster_means = m, cluster_assignments = a,
         num_iter = 4L, burn_in = 2L, thin = 1L, k = 2L,
         sample_ages = TRUE, age_ranges = p$age.matrix,
         sampling_ages = p$midpoint_ages),
    "[Ii]nvalid.*[Ii]nitial|[Ii]nvalid.*start"
  )
})

test_that("the likelihood also excludes time travellers at k equals one", {
  p <- getFromNamespace("prepare_clustering_data", "admixplorer")(
    create_feasible_ancient_data()
  )
  likelihood <- getFromNamespace("nig_log_likelihood", "admixplorer")
  ll <- likelihood(rep(1L, length(p$dates_original)),
                   cluster_means = min(p$age.matrix[, 1L]) - 1,
                   std_errors = p$std_errors,
                   dates = p$dates_original + p$midpoint_ages,
                   n = length(p$dates_original),
                   age_ranges = p$age.matrix)
  expect_identical(ll$total, -Inf)
})

test_that("the real feasibility search returns usable occupied-cluster starts", {
  p <- getFromNamespace("prepare_clustering_data", "admixplorer")(
    create_feasible_ancient_data()
  )
  search <- getFromNamespace("mc_feasibility_with_solutions", "admixplorer")
  means_fn <- getFromNamespace("compute_cluster_means", "admixplorer")
  validate <- getFromNamespace("valid_mcmc_start", "admixplorer")
  result <- NULL
  invisible(capture.output(
    result <- search(dates = p$dates_original, std_errors = p$std_errors,
                     k = 2L, age_ranges = p$age.matrix,
                     n_draws = 3L, verbose = FALSE)
  ))
  expect_gt(result$n_feasible, 0L)
  expect_length(result$solutions, result$n_feasible)
  for (solution in result$solutions) {
    m <- means_fn(solution$assignments,
                  p$dates_original + solution$sampling_ages,
                  p$std_errors, 2L)
    expect_true(validate(solution$assignments, m, solution$sampling_ages,
                         p$age.matrix, 2L))
  }
})
