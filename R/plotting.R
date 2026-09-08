# Plotting Functions
#' @importFrom stats sd
#' @importFrom grDevices pdf dev.off rainbow colorRampPalette
#' @importFrom ggplot2 ggplot geom_point aes geom_hline annotate xlab ylab theme_minimal .data
#' @importFrom gplots heatmap.2
#' @importFrom utils write.csv


plot_mcmc_results <- function(mcmc_result, k, outfile_prefix, plot = TRUE) {
  if (!plot) return()
  if (k==1) return()
  # Define colors for up to 8 clusters
  colors <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#F0E442", "#56B4E9", "#999999")

  if (k > length(colors)) {
    # Generate more colors if needed
    colors <- rainbow(k)
  }

  # MCMC trace plot
  if (!is.null(mcmc_result$cluster_mean_samples)) {
    result_df <- as.data.frame(mcmc_result$cluster_mean_samples)
    output.outfile <- paste0(outfile_prefix, ".", k, "clust.pdf")
    pdf(output.outfile)

    # Create the base plot
    p <- ggplot(result_df)

    # Add points and horizontal lines for each cluster
    for (i in 1:k) {
      col_name <- paste0("V", i)
      if (col_name %in% names(result_df)) {
        p <- p +
          geom_point(aes(x = 1:nrow(result_df), y = .data[[col_name]]),
                     colour = colors[i]) +
          geom_hline(aes(yintercept = mean(.data[[col_name]], na.rm = TRUE)),
                     colour = colors[i])
      }
    }

    # Create cluster labels for annotation
    cluster_labels <- character(k)
    cluster_names <- c("blue", "orange", "green", "red", "pink", "yellow", "light blue", "gray")

    for (i in 1:k) {
      col_name <- paste0("V", i)
      if (col_name %in% names(result_df)) {
        color_name <- if (i <= length(cluster_names)) cluster_names[i] else paste("color", i)
        cluster_labels[i] <- sprintf("Cluster %d (%s): mean = %.2f, sd = %.2f",
                                     i, color_name,
                                     mean(result_df[[col_name]], na.rm = TRUE),
                                     sd(result_df[[col_name]], na.rm = TRUE))
      }
    }

    # Add annotation
    p <- p +
      annotate("text",
               x = Inf, y = Inf,
               label = paste(cluster_labels, collapse = "\n"),
               hjust = 1.1, vjust = 1.2, size = 6) +
      xlab("Iteration") +
      ylab("Cluster mean") +
      theme_minimal()

    print(p)
    dev.off()
  }

  if (!is.null(mcmc_result$log_likelihood_trace_lambda10)) {
    ll_df <- data.frame(
      iter = seq_along(mcmc_result$log_likelihood_trace_lambda10),
      ll = mcmc_result$log_likelihood_trace_lambda10
    )

    keep <- seq(1, nrow(ll_df), by = 100)
    ll_df_sub <- ll_df[keep, , drop = FALSE]

    output.outfile <- paste0(outfile_prefix, ".", k, "clust.loglik.lambda10.pdf")
    pdf(output.outfile, width = 9, height = 5)

    p <- ggplot(ll_df_sub, aes(x = iter, y = ll)) +
      geom_point(alpha = 0.6) +
      xlab("Iteration") +
      ylab("Log-likelihood (lambda = 10)") +
      theme_minimal()

    print(p)
    dev.off()
  }
  # Co-clustering heatmap
  if (!is.null(mcmc_result$co_clustering_matrix)) {
    coincidence <- as.matrix(mcmc_result$co_clustering_matrix)
    clusters <- mcmc_result$result$cluster
    ordering <- order(clusters)
    ordered_matrix <- coincidence[ordering, ordering]

    clusters_ord <- clusters[ordering]                       # vector of sorted cluster labels
    cluster_sizes <- as.numeric(table(clusters_ord))         # tabulate in the same order
    cluster_bounds <- cumsum(c(0, cluster_sizes))
    output.outfile <- paste0(outfile_prefix, ".", k, "clust.coincidence.pdf")
    pdf(output.outfile, width = 8, height = 6)
    heatmap.2(
      ordered_matrix,
      Rowv = FALSE,
      Colv = FALSE,
      dendrogram = "none",
      colsep = cluster_bounds,
      rowsep = cluster_bounds,
      trace = "none",
      col = colorRampPalette(c("white", "red"))(100),
      breaks = seq(0, 1, length.out = 101),
      margins = c(5, 10),
      cexRow = 0.7,
      cexCol = 0.7,
      key = TRUE,
      density.info = "none",
      main = paste("Co-clustering Matrix (k =", k, ")")
    )
    dev.off()

    # Save CSV
    output.outfile <- paste0(outfile_prefix, ".", k, "clust.coincidence.csv")
    write.csv(coincidence, output.outfile)
  }
}

#' Plot per-individual likelihood heatmap across k
#'
#' @param all_mcmc_results List of MCMC results (one entry per k, as used in generate_final_output)
#' @param pop.vec Population vector, same order as individuals in mcmc_result
#' @param outfile Output file prefix
#' @param outlier_sd_thresh SD threshold below the row mean to flag an individual as an outlier (default 2)
#' @return Data frame of flagged outlier individuals (invisibly)
#' @export
#' @importFrom gplots heatmap.2
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom utils write.csv
plot_individual_likelihood_heatmap <- function(all_mcmc_results, pop.vec, outfile,
                                               outlier_sd_thresh = 2) {

  completed_ks <- names(all_mcmc_results)[!sapply(all_mcmc_results, is.null)]
  completed_ks <- completed_ks[order(as.numeric(completed_ks))]

  if (length(completed_ks) < 2) {
    cat("Need at least two completed k values to build likelihood heatmap - skipping.\n")
    return(invisible(NULL))
  }

  # Build matrix: rows = k (pulses), columns = individuals
  ll_list <- lapply(completed_ks, function(k) {
    mcmc_result <- all_mcmc_results[[k]]
    mcmc_result$result$final_log_likelihood_best_sampling_ages_perind
  })

  n_ind <- length(ll_list[[1]])
  ll_matrix <- do.call(rbind, ll_list)
  rownames(ll_matrix) <- paste0("k=", completed_ks)
  colnames(ll_matrix) <- if (!is.null(pop.vec)) {
    make.unique(paste0(seq_len(n_ind), "_", pop.vec))
  } else {
    as.character(seq_len(n_ind))
  }

  # Flag outliers: individuals whose likelihood is far below the row mean
  # in every single k (i.e. consistently poor fit regardless of model complexity)
  row_z <- t(scale(t(ll_matrix)))  # z-score within each k row
  is_low_everywhere <- apply(row_z, 2, function(col) all(col < -outlier_sd_thresh))
  outlier_inds <- colnames(ll_matrix)[is_low_everywhere]

  output.outfile <- paste0(outfile, ".individual_likelihood_heatmap.pdf")
  pdf(output.outfile, width = 10, height = 6)

  heatmap.2(
    ll_matrix,
    Rowv = FALSE,
    Colv = FALSE,
    dendrogram = "none",
    trace = "none",
    col = colorRampPalette(c("darkred", "red", "white"))(100),
    margins = c(8, 8),
    cexRow = 0.9,
    cexCol = 0.6,
    key = TRUE,
    symbreaks = FALSE,
    symkey = FALSE,
    density.info = "none",
    key.xlab = "Log-likelihood",
    lhei = c(1.5, 5),
    main = "Per-individual log-likelihood across pulse number (k)",
    xlab = "Individual",
    ylab = "Number of pulses (k)"
  )
  dev.off()
  # Save CSV
  write.csv(ll_matrix, paste0(outfile, ".individual_likelihood_heatmap.csv"))

  if (length(outlier_inds) > 0) {
    cat("Potential outlier individuals (low likelihood across all tested k):",
        paste(outlier_inds, collapse = ", "), "\n")
  } else {
    cat("No individuals flagged as consistent outliers across tested k.\n")
  }

  invisible(data.frame(individual = outlier_inds, stringsAsFactors = FALSE))
}

