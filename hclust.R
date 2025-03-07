

TSS_cluster_hclust <- function(DT, max_cluster_size = 100, use_weights = TRUE, plot_results = FALSE) {
  # DT must have columns: seqnames, strand, position, count
  required_cols <- c("seqnames", "strand", "position", "count")
  if (!all(required_cols %in% colnames(DT))) {
    stop("DT must contain the following columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Exclude zero counts
  DT <- DT[count > 0]
  
  # Process each strand separately
  clustered_DT <- DT[, {
    dt_strand <- .SD[order(position)]
    if (nrow(dt_strand) > 1) {
      d <- dist(dt_strand$position)
      hc <- hclust(d, method = "complete")
      clusters <- cutree(hc, h = max_cluster_size)
    } else {
      clusters <- 1
    }
    dt_strand[, cluster := clusters]
    
    # Compute cluster-level statistics
    dt_strand[, cluster_start := min(position), by = cluster]
    dt_strand[, cluster_end := max(position), by = cluster]
    dt_strand[, cluster_width := cluster_end - cluster_start + 1]
    dt_strand[, cluster_center := round(mean(position), 0), by = cluster]
    dt_strand[, cluster_peak := position[which.max(count)], by = cluster]
    dt_strand[, cluster_total_count := sum(count), by = cluster]
    
    # Enforce consistent types (as numeric/double)
    dt_strand[, c("position", "cluster", "cluster_start", "cluster_end", 
                  "cluster_width", "cluster_center", "cluster_peak", "cluster_total_count") :=
                lapply(.SD, as.numeric),
              .SDcols = c("position", "cluster", "cluster_start", "cluster_end", 
                          "cluster_width", "cluster_center", "cluster_peak", "cluster_total_count")]
    dt_strand
  }, by = .(seqnames, strand)]
  
  if (plot_results) {
    library(ggplot2)
    p <- ggplot(clustered_DT, aes(x = position, y = count, color = factor(cluster))) +
      geom_point() +
      geom_vline(aes(xintercept = cluster_center), linetype = "dashed", color = "black") +
      facet_wrap(~ strand, scales = "free_x") +
      theme_minimal() +
      labs(title = "Hierarchical Clustering of TSSs",
           x = "Genomic Position (bp)",
           y = "Read Count",
           color = "Cluster")
    print(p)
  }
  
  return(clustered_DT)
}



all_clusters <- prime5.counts[correct_tss == TRUE, {
  # Subset the necessary columns and rename 'pos' to 'position'
  dt_to_cluster <- .SD[, .(seqnames, strand, position = pos, count)]
  
  # Apply the hierarchical clustering function
  clustered <- TSS_cluster_hclust(dt_to_cluster, max_cluster_size = 100, use_weights = TRUE)
  
  # Get unique clusters per strand
  unique_clusters <- unique(clustered[, .(seqnames, strand, cluster, 
                                          cluster_start, cluster_end, 
                                          cluster_peak, cluster_width, 
                                          cluster_total_count)])
  
  # Filter clusters using your filtering function (assumed to be defined)
  filtered <- filter_clusters(unique_clusters, cluster_min_ratio = 0.0005, cluster_min_count = 0.0005)
  
  # Wrap the filtered result in a list so that each group's result is a single list-column entry
  list(clusters = list(filtered))
}, by = .(hpi, cell_line)]

