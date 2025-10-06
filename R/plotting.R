# Plotting Functions
#' @importFrom stats sd
#' @importFrom grDevices pdf dev.off rainbow colorRampPalette
#' @importFrom ggplot2 ggplot geom_point aes geom_hline annotate xlab ylab theme_minimal .data
#' @importFrom gplots heatmap.2
#' @importFrom utils write.csv

plot_mcmc_results <- function(mcmc_result, k, outfile_prefix, plot = TRUE) {
  if (!plot) return()
  
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
  
  # Co-clustering heatmap
  if (!is.null(mcmc_result$co_clustering_matrix)) {
    coincidence <- as.matrix(mcmc_result$co_clustering_matrix)
    clusters <- mcmc_result$result$cluster
    ordering <- order(clusters)
    ordered_matrix <- coincidence[ordering, ordering]
    cluster_counts <- table(clusters)
    cluster_sizes <- as.numeric(cluster_counts[order(unique(clusters))])
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

