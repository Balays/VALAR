

library(mclust)
library(data.table)
library(ggplot2)

Mclust_stranded <- function(DT, sep_strands = TRUE, use_weights = TRUE, include_zeros = FALSE, plot_results = FALSE, G = NULL) {
  # DT: data.table with at least 'seqnames', 'strand', 'position', 'count'
  # sep_strands: if TRUE, cluster '+' and '-' strands separately
  # use_weights: if TRUE, use counts as weights in clustering
  # include_zeros: if TRUE, include positions with zero counts
  # plot_results: if TRUE, plot the clustering results
  # G: Number of mixture components (clusters). If NULL, let Mclust select the best number.

  # Check for required columns
  required_cols <- c('seqnames', 'strand', 'position', 'count')
  if (!all(required_cols %in% colnames(DT))) {
    stop(paste('DT must contain the following columns:', paste(required_cols, collapse = ', ')))
  }

  # Handle seqnames
  seqnames <- unique(DT[, seqnames])
  if (length(seqnames) > 1) {
    stop('There are more than 1 contigs in the data, which is not supported yet. Please run the function separately for each contig.')
  }

  if (!include_zeros) {
    DT <- DT[count > 0, ]
    message('Positions with a count > 0 will be used only.')
  } else {
    message('Including positions with zero counts.')
  }

  if (!use_weights) {
    DT[, count := 1]
    message('Weights will not be used; each position will have the same effect on the clustering.')
  } else {
    message('Using counts as weights in clustering.')
  }

  if (sep_strands) {
    message('Clustering will be performed separately for "+" and "-" strands.')

    ### Positive strand
    DT_pos <- DT[strand == '+', ]
    if (nrow(DT_pos) > 1) {
      mclust_result_pos <- Mclust(data = DT_pos$position, G = G, weights = DT_pos$count)
      DT_pos$cluster <- mclust_result_pos$classification
    } else if (nrow(DT_pos) == 1) {
      DT_pos$cluster <- 1
    } else {
      DT_pos <- NULL
    }

    ### Negative strand
    DT_neg <- DT[strand == '-', ]
    if (nrow(DT_neg) > 1) {
      mclust_result_neg <- Mclust(data = DT_neg$position, G = G, weights = DT_neg$count)
      DT_neg$cluster <- mclust_result_neg$classification
    } else if (nrow(DT_neg) == 1) {
      DT_neg$cluster <- 1
    } else {
      DT_neg <- NULL
    }

    # Combine
    DT <- rbind(DT_pos, DT_neg)

  } else {
    message('Clustering will be performed on both strands together.')
    if (nrow(DT) > 1) {
      mclust_result <- Mclust(data = DT$position, G = G, weights = DT$count)
      DT$cluster <- mclust_result$classification
    } else if (nrow(DT) == 1) {
      DT$cluster <- 1
    } else {
      stop('No data available for clustering.')
    }
  }

  # Assign cluster IDs
  DT[, cluster_ID := paste0('cluster_', cluster)]

  # Compute cluster statistics
  DT[, cluster_start := min(position), by = .(seqnames, strand, cluster_ID)]
  DT[, cluster_end := max(position), by = .(seqnames, strand, cluster_ID)]
  DT[, cluster_width := cluster_end - cluster_start + 1]
  DT[, cluster_center := round(mean(position), 0), by = .(seqnames, strand, cluster_ID)]
  DT[, cluster_peak := position[which.max(count)], by = .(seqnames, strand, cluster_ID)]

  # Order the data
  DT <- DT[order(position)]

  # Optionally plot the results
  if (plot_results) {
    ggp <- ggplot(DT, aes(x = position, y = count, color = cluster_ID)) +
      geom_point(size = 2) +
      geom_errorbarh(aes(xmin=cluster_start, xmax=cluster_end,  y = max(count)/10 *-1, color = cluster_ID)) +
      geom_point(data=DT[position == cluster_peak, ], aes(position, count), color='black', fill='lightgreen', size = 3, shape=21) +
      geom_vline(aes(xintercept=cluster_center)) +
      facet_wrap(~ strand, scales = 'free_x') +
      theme_minimal() +
      labs(title = 'Clustering Results using mclust',
           x = 'Position',
           y = 'Count',
           color = 'Cluster ID')
    print(ggp)
  }

  return(DT)
}

if(!exists('dontrun')) { dontrun <- T }

if(!dontrun) {

  DT <- trids.toclust[,':='(seqnames=seqnames, strand=strand, position=TR_end, count=sum_count, ID=seq_along(TR_end))][order(position)]

  # Assuming DT is your data.table with the necessary columns
  DT_mclust <- Mclust_stranded(
    DT, #[position >= 1908 & position <= 2012 & strand == '-'],
    sep_strands = FALSE,
    use_weights = TRUE,
    include_zeros = FALSE,
    plot_results = TRUE,
    G = NULL  # Let Mclust decide the optimal number of clusters
  )
}
