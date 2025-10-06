# Feasibility and Assignment Testing
#' @importFrom stats runif quantile kmeans


# Helper: test one assignment for ancient time-traveler violations
test_valid_assignment <- function(
    full_assign, dates, std_errors, age_ranges, sampling_ages, verbose=FALSE
) {
  k <- max(full_assign)
  means <- numeric(k)
  bad_inds <- integer(0)
  
  for (j in seq_len(k)) {
    members <- which(full_assign == j)
    # Weighted mean over all members
    w    <- 1 / (std_errors[members]^2)
    dates_adj <- dates[members] + sampling_ages[members]
    m_j  <- sum(dates_adj * w) / sum(w)
    means[j] <- m_j
    
    # Check only ancient members
    anc <- members[age_ranges[members,2] > 0]
    viol <- anc[m_j < age_ranges[anc,1]]
    if (length(viol)) {
      bad_inds <- c(bad_inds, viol)
      if (verbose) {
        cat(sprintf("  Cluster %d mean = %.2f  Time travelers: %s\n",
                    j, m_j, paste(viol, collapse=",")))
      }
    } else if (verbose) {
      cat(sprintf("  Cluster %d mean = %.2f  Valid\n", j, m_j))
    }
  }
  
  ok <- length(bad_inds) == 0
  if (ok && verbose) {
    cat("  SUCCESS: Valid clustering found!\n")
    cat("  Cluster means:", paste(round(means,2), collapse=" "), "\n")
  }
  list(ok = ok, bad = sort(unique(bad_inds)), means = means)
}

# Comprehensive feasibility tester across multiple strategies
test_clustering_feasibility <- function(
    dates, std_errors, k, age_ranges, sampling_ages, max_attempts=100, verbose=FALSE
) {
  n <- length(dates)
  if (verbose) cat("=== Testing clustering feasibility for k =", k, "===\n")
  
  ancient <- which(age_ranges[,2] > 0)
  modern  <- setdiff(seq_len(n), ancient)
  if (verbose) cat("Ancient:", length(ancient), "| Modern:", length(modern), "\n")
  
  if (!length(ancient)) {
    if (verbose) cat("No ancient samples  trivially feasible\n")
    return(list(feasible=TRUE, strategy="none_needed"))
  }
  
  # Prepare adjusted dates
  dates_adj_all <- dates + sampling_ages
  
  # 1) Age-based bins on ancient sampling upper bounds
  if (verbose) cat("\nStrategy 1: age-based\n")
  ages_up <- age_ranges[ancient,2]
  breaks_age <- unique(quantile(ages_up, probs=seq(0,1,length.out=k+1)))
  if (length(breaks_age) <= k) {
    asn <- rep(1,n)
    for (i in seq_len(min(k-1,length(ancient)))) {
      asn[ancient[i]] <- i+1
    }
  } else {
    bins <- cut(ages_up, breaks=breaks_age, include.lowest=TRUE, labels=FALSE)
    asn <- rep(1,n)
    asn[ancient] <- bins
  }
  res1 <- test_valid_assignment(asn, dates, std_errors, age_ranges, sampling_ages, verbose)
  if (res1$ok) {
    if (verbose) {
      cat(" age-based works!\nFull assignment:\n")
      print(asn)
    }
    return(list(feasible=TRUE, strategy="age_based", assignments=asn))
  }
  
  # 2) Data-based bins on adjusted ancient dates
  if (verbose) cat("\nStrategy 2: data-based\n")
  dat_a <- dates_adj_all[ancient]
  breaks_dat <- unique(quantile(dat_a, probs=seq(0,1,length.out=k+1)))
  if (length(breaks_dat) <= k) {
    asn <- sample(1:k, n, replace=TRUE)
  } else {
    bins <- cut(dat_a, breaks=breaks_dat, include.lowest=TRUE, labels=FALSE)
    asn <- rep(1,n)
    asn[ancient] <- bins
  }
  res2 <- test_valid_assignment(asn, dates, std_errors, age_ranges, sampling_ages, verbose)
  if (res2$ok) {
    if (verbose) {
      cat(" data-based works!\nFull assignment:\n")
      print(asn)
    }
    return(list(feasible=TRUE, strategy="data_based", assignments=asn))
  }
  
  # 3) Constraint-aware: isolate the k-1 oldest ancient samples
  if (verbose) cat("\nStrategy 3: constraint-aware\n")
  sorted_old <- order(age_ranges[ancient,2], decreasing=TRUE)
  asn <- rep(1,n)
  for (i in seq_len(min(k-1,length(sorted_old)))) {
    asn[ancient[sorted_old[i]]] <- i+1
  }
  res3 <- test_valid_assignment(asn, dates, std_errors, age_ranges, sampling_ages, verbose)
  if (res3$ok) {
    if (verbose) {
      cat(" constraint-aware works!\nFull assignment:\n")
      print(asn)
    }
    return(list(feasible=TRUE, strategy="constraint_aware", assignments=asn))
  }
  
  # 4) Reverse constraint-aware
  if (verbose) cat("\nStrategy 4: reverse_constraint\n")
  sorted_young <- order(age_ranges[ancient,2], decreasing=FALSE)
  asn <- rep(1,n)
  for (i in seq_len(min(k-1,length(sorted_young)))) {
    asn[ancient[sorted_young[i]]] <- i+1
  }
  res4 <- test_valid_assignment(asn, dates, std_errors, age_ranges, sampling_ages, verbose)
  if (res4$ok) {
    if (verbose) {
      cat(" reverse_constraint works!\nFull assignment:\n")
      print(asn)
    }
    return(list(feasible=TRUE, strategy="reverse_constraint", assignments=asn))
  }
  
  # 5) k-means on all adjusted dates
  if (verbose) cat("\nStrategy 5: kmeans\n")
  if (n >= k) {
    km <- try(kmeans(dates_adj_all, centers=k, nstart=5)$cluster, silent=TRUE)
    if (!inherits(km,"try-error")) {
      res5 <- test_valid_assignment(km, dates, std_errors, age_ranges, sampling_ages, verbose)
      if (res5$ok) {
        if (verbose) {
          cat(" kmeans works!\nFull assignment:\n")
          print(km)
        }
        return(list(feasible=TRUE, strategy="kmeans", assignments=km))
      }
    }
  }
  
  # 6) Random attempts
  if (verbose) cat("\nStrategy 6: random\n")
  for (i in seq_len(max_attempts)) {
    rnd <- sample(1:k, n, replace=TRUE)
    rr  <- test_valid_assignment(rnd, dates, std_errors, age_ranges, sampling_ages, verbose=FALSE)
    if (rr$ok) {
      if (verbose) {
        cat(" random works on attempt", i, "\nFull assignment:\n")
        print(rnd)
      }
      return(list(feasible=TRUE, strategy=paste0("random_",i), assignments=rnd))
    }
  }
  
  # 7) Exhaustive (small n<=8,k<=3)
  if (n <= 8 && k <= 3) {
    if (verbose) cat("\nStrategy 7: exhaustive\n")
    combos <- expand.grid(rep(list(1:k), n))
    for (i in seq_len(min(nrow(combos),1000))) {
      asn <- as.numeric(combos[i,])
      rr  <- test_valid_assignment(asn, dates, std_errors, age_ranges, sampling_ages, verbose=FALSE)
      if (rr$ok) {
        if (verbose) {
          cat(" exhaustive works!\nFull assignment:\n")
          print(asn)
        }
        return(list(feasible=TRUE, strategy="exhaustive", assignments=asn))
      }
    }
  }
  
  # No valid clustering
  if (verbose) cat("\nFAILURE: no valid clustering for k =", k, "\n")
  return(list(
    feasible = FALSE,
    reason = "none_worked",
    problematic_inds = if (exists("rr")) rr$bad else integer(0)
  ))
}














mc_feasibility_with_solutions <- function(
    dates, std_errors, k,
    age_ranges, n_draws = 500,
    test_fn = test_clustering_feasibility,
    ...
) {
  n <- length(dates)
  feasible_flags  <- logical(n_draws)
  solutions_ages  <- vector("list", n_draws)
  solutions_assign <- vector("list", n_draws)
  
  for (i in seq_len(n_draws)) {
    # 1) Draw one sampling-age vector
    samp_ages <- runif(n, age_ranges[,1], age_ranges[,2])
    
    # 2) Test feasibility  also return assignments if feasible
    res <- test_fn(dates, std_errors, k, age_ranges, samp_ages, ...)
    
    feasible_flags[i] <- res$feasible
    if (res$feasible) {
      solutions_ages[[i]]  <- samp_ages
      solutions_assign[[i]] <- res$assignments
    }
  }
  
  # Keep only the successful draws
  idx_ok <- which(feasible_flags)
  list(
    n_draws         = n_draws,
    n_feasible      = length(idx_ok),
    prop_feasible   = length(idx_ok)/n_draws,
    solutions       = lapply(idx_ok, function(i) {
      list(
        sampling_ages = solutions_ages[[i]],
        assignments   = solutions_assign[[i]]
      )
    })
  )
}

