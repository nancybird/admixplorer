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
  # Choose the first eligible individual, making repair tests deterministic.
  # This harness supplies only one fallback solution, so its selection
  # uses n = 1 and is not recorded as an empty-cluster repair.
  sandbox$sample.int <- function(n, size, replace = FALSE, ...) {
    if (n > 1L) {
      event("greedy")
    }
    rep.int(1L, size)
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

test_that("collapsed EM assignments move only one individual into the empty cluster", {
  original <- rep(1L, 10L)

  h <- .greedy_gate_harness(
    original,
    modern = TRUE
  )

  expect_no_error(h$run())

  events <- h$seen$events

  expect_identical(events[1L], "EM")
  expect_true(
    match("greedy", events) < match("validate", events)
  )
  expect_false("feasibility" %in% events)
  expect_length(h$seen$starts, 1L)

  repaired <- h$seen$starts[[1L]]$assignments
  counts <- tabulate(repaired, nbins = 2L)

  # Exactly one individual moves; all other EM assignments survive.
  expect_equal(sum(repaired != original), 1L)

  # The formerly empty cluster gets exactly one individual.
  expect_equal(counts, c(9L, 1L))

  # The deterministic mock selects the first eligible individual.
  expected <- original
  expected[1L] <- 2L
  expect_equal(repaired, expected)

  # The actual repaired state is validated before entering MCMC.
  expect_true(
    match("validate", events) < match("MCMC", events)
  )
})

test_that("empty-cluster repair preserves occupied singleton clusters", {
  # Cluster 2 is empty; cluster 3 already contains a singleton.
  original <- c(rep(1L, 9L), 3L)

  h <- .greedy_gate_harness(
    original,
    modern = TRUE
  )

  expect_no_error(h$run())
  expect_true("greedy" %in% h$seen$events)
  expect_false("feasibility" %in% h$seen$events)
  expect_length(h$seen$starts, 1L)

  repaired <- h$seen$starts[[1L]]$assignments

  # Only one individual should move into the empty cluster.
  expect_equal(sum(repaired != original), 1L)
  expect_equal(
    tabulate(repaired, nbins = 3L),
    c(8L, 1L, 1L)
  )

  # The existing singleton must not be used as a donor.
  expect_equal(repaired[10L], 3L)

  expected <- original
  expected[1L] <- 2L
  expect_equal(repaired, expected)
})

test_that("several empty clusters require only one move each", {
  # Fit k=5, but EM occupies only clusters 1 and 5.
  # Cluster 5 is a singleton and must remain untouched.
  original <- c(rep(1L, 9L), 5L)

  h <- .greedy_gate_harness(
    original,
    modern = TRUE
  )

  expect_no_error(h$run())
  expect_false("feasibility" %in% h$seen$events)
  expect_length(h$seen$starts, 1L)

  repaired <- h$seen$starts[[1L]]$assignments
  counts <- tabulate(repaired, nbins = 5L)

  # Three empty clusters: exactly three individuals move.
  expect_equal(sum(repaired != original), 3L)

  # Every requested cluster is occupied.
  expect_true(all(counts > 0L))
  expect_equal(counts, c(6L, 1L, 1L, 1L, 1L))

  # Existing singleton is retained.
  expect_equal(repaired[10L], 5L)

  # With the deterministic mock, individuals 1, 2 and 3 move.
  expected <- original
  expected[1:3] <- 2:4
  expect_equal(repaired, expected)
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
