# consider between 2 and 10 clusters
k_values <- 2:12


wss_for_k <- function(data, k, hc) {
  clusters <- cutree(hc, k = k)
  # Compute cluster centroids
  centroids <- aggregate(data, by = list(cluster = clusters), FUN = mean)
  # Remove the cluster column
  centroid_coords <- centroids[ , -1, drop=FALSE]

  # For each point, get its cluster centroid and compute squared distance
  sum_sq <- 0
  for (i in seq_len(nrow(data))) {
    cluster_id <- clusters[i]
    centroid_vec <- as.numeric(centroid_coords[cluster_id, ])
    point_vec <- as.numeric(data[i, ])
    sum_sq <- sum_sq + sum((point_vec - centroid_vec)^2)
  }

  return(sum_sq)
}


wss_values <- sapply(k_values, function(k) wss_for_k(countData, k, hc))


plot(k_values, wss_values, type = "b", pch = 19,
     xlab = "Number of clusters (k)", ylab = "Within-Cluster Sum of Squares",
     main = "Elbow Plot for Hierarchical Clustering")
