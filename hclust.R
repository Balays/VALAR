

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



filter_clusters <- function(hclust_res, cluster_min_ratio = NULL, cluster_min_count = NULL) {
  # Calculate cluster_ratio if not already present
  if (!"cluster_ratio" %in% colnames(hclust_res)) {
    hclust_res[, cluster_ratio := cluster_total_count / sum(cluster_total_count)]
  }
  
  # Start with the full set of clusters
  filtered_clusters <- copy(hclust_res)
  
  # Apply ratio filter if provided
  if (!is.null(cluster_min_ratio)) {
    filtered_clusters <- filtered_clusters[cluster_ratio >= cluster_min_ratio]
  }
  
  # Apply count filter if provided
  if (!is.null(cluster_min_count)) {
    filtered_clusters <- filtered_clusters[cluster_total_count >= cluster_min_count]
  }
  
  return(filtered_clusters)
}



all_clusters <- prime5.counts[correct_tss == TRUE, {
  # Subset the necessary columns and rename 'pos' to 'position'
  dt_to_cluster <- .SD[, .(seqnames, strand, position = pos, count)]
  
  # Apply the hierarchical clustering function
  clustered <- TSS_cluster_hclust(dt_to_cluster, max_cluster_size = 100, use_weights = TRUE)
  
}, by = .(hpi, cell_line)]

# Get unique clusters per strand
unique_clusters <- unique(all_clusters[, .(seqnames, strand, cluster, 
                                          cluster_start, cluster_end, 
                                          cluster_peak, cluster_width, 
                                          cluster_total_count)])



# Filter clusters using your filtering function (assumed to be defined)
filtered_clusters <- all_clusters_long[, {
  # Subset the necessary columns and rename 'pos' to 'position'
  dt_to_cluster <- .SD[, .(seqnames, strand, 
                           cluster, cluster_start, cluster_end, cluster_peak, cluster_width, cluster_total_count)]
  
  ## Filter
  filtered<- filter_clusters(dt_to_cluster, cluster_min_ratio = cluster_min_ratio, cluster_min_count = cluster_min_count)
  
 }, by = .(hpi, cell_line)]



# Get unique clusters per strand
unique_filt_clusters <- unique(filtered_clusters[, .(seqnames, strand, cluster, 
                                                cluster_start, cluster_end, 
                                                cluster_peak, cluster_width, 
                                                cluster_total_count, cluster_ratio)])




#### Meta Clusters

# Order the clusters
setorder(filtered_clusters, seqnames, strand, cluster_start)


# For each contig and strand, assign a meta_cluster ID by merging overlapping intervals.
filtered_clusters[, meta_cluster := {
  meta = integer(.N)
  current_meta = 1L
  current_end = cluster_end[1]
  meta[1] = current_meta
  if (.N > 1) {
    for (i in 2:.N) {
      # If this cluster overlaps with the previous meta-cluster (its start is <= current_end)
      if (cluster_start[i] <= current_end) {
        meta[i] = current_meta
        # Update the current meta-cluster end if needed
        current_end = max(current_end, cluster_end[i])
      } else {
        current_meta = current_meta + 1L
        meta[i] = current_meta
        current_end = cluster_end[i]
      }
    }
  }
  meta
}, by = .(seqnames, strand)]


fwrite(filtered_clusters, paste0(outdir, '/filtered_clusters.tsv'), sep='\t')

filtered_clusters <- fread(paste0(outdir, '/filtered_clusters.tsv'), na.strings = '')

meta_summary <- filtered_clusters[, .(
  meta_cluster_start = min(cluster_start),
  meta_cluster_end   = max(cluster_end),
  consensus_peak     = mean(cluster_peak),  # or weighted average if desired
  peak_shift         = max(cluster_peak) - min(cluster_peak),
  min_peak = min(cluster_peak),
  max_peak = max(cluster_peak),
  mean_peak = mean(cluster_peak),
  width_range        = max(cluster_width) - min(cluster_width),
  sample_count       = .N,
  hpi_samples        = paste(unique(hpi), collapse = ", "),
  cell_lines         = paste(unique(cell_line), collapse = ", "),
  total_read_count   = sum(cluster_total_count)
), by = .(seqnames, strand, meta_cluster)]

filtered_clusters <- merge(filtered_clusters, meta_summary, by=c('seqnames', 'strand', 'meta_cluster'))


meta_summary

ggplot(meta_summary) +
  geom_histogram(aes(width_range), binwidth =  1)

ggplot(filtered_clusters) +
  geom_histogram(aes(width_range), binwidth =  1)





# Compute the summary for each meta_cluster
meta_peaks <- all_clusters_long[, .(
  min_peak = min(cluster_peak),
  max_peak = max(cluster_peak),
  mean_peak = mean(cluster_peak)
), by = meta_cluster]

#all_clusters_long <- merge(all_clusters_long, meta_peaks, by=c('seqnames', 'strand', 'meta_cluster'))



### Plottin clusters

# --- Plot 1: Cluster Peak Positions across Meta Clusters ---
p1 <- ggplot(all_clusters_long, aes(x = factor(meta_cluster), y = cluster_peak, color = hpi, shape = cell_line)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 4, size = 3, color = "black") +
  labs(title = "Cluster Peak Positions across Meta Clusters",
       x = "Meta Cluster",
       y = "Peak Position (bp)") +
  theme_minimal()
print(p1)

# --- Plot 2: Cluster Widths across Meta Clusters ---
p2 <- ggplot(all_clusters_long, aes(x = factor(meta_cluster), y = cluster_width, color = hpi, shape = cell_line)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 4, size = 3, color = "black") +
  labs(title = "Cluster Widths across Meta Clusters",
       x = "Meta Cluster",
       y = "Cluster Width (bp)") +
  theme_minimal()
print(p2)


# --- Plot 3: Meta Cluster Peak Range with Error Bars ---

p3 <- ggplot() +
  # Error bars showing the min and max peak for each meta_cluster
  geom_errorbar(data = meta_peaks, aes(x = factor(meta_cluster), ymin = min_peak, ymax = max_peak),
                width = 0.2, color = "gray50") +
  # Jittered individual peak points, colored by hpi and shaped by cell_line
  geom_point(data = all_clusters_long, aes(x = factor(meta_cluster), y = cluster_peak, color = hpi, shape = cell_line),
             position = position_jitter(width = 0.2), size = 2, alpha = 0.7) +
  # Mean peak point for each meta_cluster
  geom_point(data = meta_peaks, aes(x = factor(meta_cluster), y = mean_peak),
             color = "black", size = 3) +
  labs(title = "Meta Cluster Peak Range",
       x = "Meta Cluster",
       y = "Cluster Peak (bp)") +
  theme_minimal()
print(p3)




p4 <- ggplot(all_clusters_long[hpi == '4h'], 
             aes(x = consensus_peak, y = log10(total_read_count), color = hpi, shape = cell_line)) +
    # Error bars showing the min and max peak for each meta_cluster
  geom_errorbar(aes(x = factor(meta_cluster), xmin = meta_cluster_start, xmax = meta_cluster_end),
                width = 0.2, color = "gray50") +
  
  #geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.7) +
  #stat_summary(fun = mean, geom = "point", shape = 4, size = 3, color = "black") +
  labs(title = "Cluster Widths across Meta Clusters",
       x = "Meta Cluster",
       y = "log10 meta cluster read count") +
  theme_minimal() +
  facet_nested(rows=vars(cell_line, hpi, strand), scales='free_y')

p4

