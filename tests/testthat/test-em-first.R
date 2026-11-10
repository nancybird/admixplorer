# Tests the requested EM -> greedy repair if needed -> gate -> fallback flow.

.greedy_gate_harness <- function(labels, modern = FALSE, chains = 1L,
                                 no_solution = FALSE, invalid_pool = FALSE) {
  ns <- asNamespace("admixplorer")
  prepare <- get("prepare_clustering_data", ns)
  mean_fn <- get("compute_cluster_means", ns)
  validate <- get("valid_mcmc_start", ns)
  p <- prepare(if (modern) create_test_data_1cluster() else create_feasible_ancient_data())
  k <- max(2L, max(labels))
  runner <- get("process_single_k", ns)
  sandbox <- new.env(parent = environment(runner))
  environment(runner) <- sandbox
  seen <- new.env(parent = emptyenv())
  seen$events <- character()
  seen$starts <- list()
  seen$candidates <- list()
  event <- function(x) seen$events <- c(seen$events, x)

  fallback <- list(
    assignments = c(1L, 1L, 1L, 2L, 1L, 1L),
    sampling_ages = p$midpoint_ages
  )
  if (invalid_pool) fallback$assignments <- c(1L, 1L, 1L, 1L, 2L, 2L)
  pool <- if (no_solution) list() else list(fallback)

  sandbox$em_clustering_original <- function(dates, std_errors, k, ...) {
    event("EM")
    r <- matrix(0, length(dates), k)
    r[cbind(seq_along(labels), labels)] <- 1
    list(responsibilities = r, cluster_means = mean_fn(labels, dates, std_errors, k))
  }
  # Force the same original greedy random choices to make branch tests stable.
  sandbox$sample <- function(x, size, replace = FALSE, ...) {
    event("greedy")
    rep(x[1L], size)
  }
  sandbox$compute_cluster_means <- function(assignments, dates, std_errors, k) {
    seen$candidates[[length(seen$candidates) + 1L]] <- assignments
    mean_fn(assignments, dates, std_errors, k)
  }
  sandbox$valid_mcmc_start <- function(...) {
    event("validate")
    validate(...)
  }
  sandbox$mc_feasibility_with_solutions <- function(..., n_draws = 200) {
    event("feasibility")
    list(n_feasible = length(pool), n_draws = n_draws, solutions = pool)
  }
  sandbox$mcmc_clustering <- function(dates, std_errors, init_cluster_means,
                                      cluster_assignments, k, sampling_ages,
                                      sample_ages, age_ranges, ...) {
    event("MCMC")
    if (!isTRUE(validate(cluster_assignments, init_cluster_means,
                         sampling_ages, age_ranges, k))) {
      stop("TEST: invalid state reached MCMC.")
    }
    seen$starts[[length(seen$starts) + 1L]] <- list(
      assignments = cluster_assignments, ages = sampling_ages,
      means = init_cluster_means
    )
    list(result = data.frame(final_log_likelihood_best_sampling_ages = -1))
  }
  list(
    seen = seen, p = p, mean_fn = mean_fn, fallback = fallback,
    run = function(k_requested = k) {
      result <- NULL
      invisible(capture.output(result <- runner(
        k_current = k_requested, prepared_data = p, chains = chains,
        sample_age_est = TRUE, plot = FALSE, outfile = tempfile()
      )))
      result
    }
  )
}

test_that("valid EM assignments go straight through the gate to MCMC", {
  a <- c(1L, 1L, 1L, 2L, 1L, 1L)
  h <- .greedy_gate_harness(a, chains = 2L)
  expect_no_error(h$run())
  expect_identical(h$seen$events[1L], "EM")
  expect_identical(sum(h$seen$events == "EM"), 1L)
  expect_false("greedy" %in% h$seen$events)
  expect_false("feasibility" %in% h$seen$events)
  expect_length(h$seen$starts, 2L)
  expect_true(match("validate", h$seen$events) < match("MCMC", h$seen$events))
  expect_equal(h$seen$starts[[1L]]$assignments, a)
})

test_that("collapsed EM assignments are repaired before validation", {
  h <- .greedy_gate_harness(rep(1L, 10L), modern = TRUE)
  expect_no_error(h$run())
  events <- h$seen$events
  expect_identical(events[1L], "EM")
  expect_true(match("greedy", events) < match("validate", events))
  expect_false("feasibility" %in% events)
  expect_length(h$seen$starts, 1L)
  expected <- rep(1L, 10L)
  ord <- order(h$p$dates_original)
  expected[ord[1:2]] <- 1:2
  expect_equal(h$seen$starts[[1L]]$assignments, expected)
  expect_true(all(tabulate(h$seen$starts[[1L]]$assignments, nbins = 2L) > 0L))
})

test_that("a partially empty k equals three EM start is also repaired", {
  h <- .greedy_gate_harness(c(rep(1L, 9L), 3L), modern = TRUE)
  expect_no_error(h$run())
  expect_true("greedy" %in% h$seen$events)
  expect_false("feasibility" %in% h$seen$events)
  expect_true(all(tabulate(h$seen$starts[[1L]]$assignments, nbins = 3L) > 0L))
})

test_that("invalid greedy repair triggers feasibility before MCMC", {
  h <- .greedy_gate_harness(rep(1L, 6L), chains = 2L)
  expect_no_error(h$run())
  events <- h$seen$events
  expect_true(match("EM", events) < match("greedy", events))
  expect_true(match("greedy", events) < match("validate", events))
  expect_true(match("validate", events) < match("feasibility", events))
  expect_true(match("feasibility", events) < match("MCMC", events))
  expect_identical(sum(events == "feasibility"), 1L)
  expect_length(h$seen$starts, 2L)
  for (start in h$seen$starts) {
    expect_equal(start$assignments, h$fallback$assignments)
    expect_equal(start$ages, h$fallback$sampling_ages)
    expect_equal(start$means, h$mean_fn(
      start$assignments, h$p$dates_original + start$ages, h$p$std_errors, 2L
    ))
  }
})

test_that("occupied but time-travelling EM assignments use feasibility directly", {
  h <- .greedy_gate_harness(c(1L, 1L, 1L, 1L, 2L, 2L))
  expect_no_error(h$run())
  expect_false("greedy" %in% h$seen$events)
  expect_identical(h$seen$events[1L], "EM")
  expect_true(match("validate", h$seen$events) < match("feasibility", h$seen$events))
  expect_true(match("feasibility", h$seen$events) < match("MCMC", h$seen$events))
  expect_equal(h$seen$starts[[1L]]$assignments, h$fallback$assignments)
})

test_that("no feasible start means no MCMC", {
  h <- .greedy_gate_harness(c(1L, 1L, 1L, 1L, 2L, 2L), no_solution = TRUE)
  expect_null(h$run())
  expect_true("feasibility" %in% h$seen$events)
  expect_false("MCMC" %in% h$seen$events)
})

test_that("invalid stored feasibility solutions are not passed to MCMC", {
  h <- .greedy_gate_harness(c(1L, 1L, 1L, 1L, 2L, 2L), invalid_pool = TRUE)
  expect_null(h$run())
  expect_true("feasibility" %in% h$seen$events)
  expect_false("MCMC" %in% h$seen$events)
})

test_that("n less than k is skipped before any initialisation", {
  h <- .greedy_gate_harness(rep(1L, 6L))
  expect_null(h$run(k_requested = 7L))
  expect_length(h$seen$events, 0L)
})
