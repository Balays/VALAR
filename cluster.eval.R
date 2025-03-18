## Clustering statistics

# Extract gene names and cut the tree
gene_clusters <- cutree(fit.gene$hclust, k = cluster_num)

# Create a data.table with genes and their cluster assignments
gene_cluster_dt <- data.table(
  gene    = rownames(countData),
  cluster = gene_clusters
)

# Extract statistics for the tree edges
edges_stats <- as.data.table(fit.gene$edges)

# Add the `au` values to the clusters (merge based on hierarchical structure if needed)
# This is indirect because `edges_stats` corresponds to branches, not directly to clusters.
# First, calculate the cluster ID at the tree cut level:
edge_cluster_map <- as.data.table(data.frame(cluster = cutree(fit.gene$hclust, k = cluster_num)))
edges_stats[, cluster := .I] # Assign cluster IDs to edges
merged_stats <- merge(gene_cluster_dt, edges_stats, by = "cluster", all.x = TRUE)

# View the resulting table
#print(merged_stats)

cluster_stats <- unique(merged_stats[,-2])


ggclust_sigstats <- ggplot(merged_stats, aes(x = gene, y = au, fill = factor(cluster))) +
  # Points for significance
  geom_point(shape=21, size=3) +
  scale_fill_manual(values = colorvec) +
  scale_y_continuous(
    limits = c(0.9, 1),
    name = "approximately unbiased (AU) p-value"
  ) +
  geom_hline(yintercept = 0.95, color='red', linetype='dashed') +
  labs(x = "Genes (ordered by silhouette width)",
       title = "AU value by Cluster") +
  facet_wrap(~ cluster, scales = "free_x", ncol=round((cluster_num/3),0)) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "italic", angle = 90, vjust = 0.5),
    axis.ticks.x = element_blank()
  )

## Silhouttes

library(cluster)

# Convert clusters to a factor for clarity
clusters_factor <- factor(clusters)

# 1. Compute Silhouette Widths
# silhouette() requires a vector of cluster assignments and a distance object.
sil <- silhouette(clusters, dist_matrix)
dimnames(sil)[[1]] <- names(clusters)

# Examine silhouette results
summary(sil)   # summary statistics of silhouette widths
#plot(sil, y = names(clusters))       # visualize silhouette widths for each gene

# The silhouette width (si) for each gene:
# si close to 1    -> gene is well-clustered
# si close to 0    -> gene lies between two clusters
# si negative      -> gene may be misclassified


# 2. Distance to Cluster Centroids
# Define a "centroid" of each cluster as the average expression profile across all replicates.
cluster_centroids <- sapply(levels(clusters_factor), function(cl) {
  # colMeans gives mean expression across samples (columns) for all genes in this cluster
  colMeans(countData[clusters_factor == cl, , drop=FALSE])
})

# Calculate Euclidean distance from each gene to its cluster centroid
gene_to_centroid_dist <- numeric(nrow(countData))
names(gene_to_centroid_dist) <- rownames(countData)

for (i in seq_len(nrow(countData))) {
  cl <- clusters_factor[i]
  gene_expr <- as.numeric(countData[i, ])
  centroid <- cluster_centroids[, cl]

  # Euclidean distance
  gene_to_centroid_dist[i] <- sqrt(sum((gene_expr - centroid)^2))
}


# gene_to_centroid_dist now holds a measure of how close each gene is to the centroid of its cluster.
# Lower distance means the gene is more "typical" of that cluster's expression pattern.




# Convert silhouette object to a data frame
sil_df <- as.data.frame(sil[, 1:3])
colnames(sil_df) <- c("cluster", "neighbor", "sil_width")

# Add gene names if available
sil_df$gene <- rownames(countData)

# Order by cluster and silhouette width
sil_df <- sil_df[order(sil_df$cluster, -sil_df$sil_width), ]

# Create a gene_order within each cluster for plotting
sil_df$gene_order <- ave(sil_df$sil_width, sil_df$cluster, FUN = function(x) seq_along(x))

# Add gene_to_centroid_dist to sil_df (assuming you've computed gene_to_centroid_dist as shown before)
# gene_to_centroid_dist is a numeric vector with names matching rownames(countData)
sil_df$gene_to_centroid_dist <- gene_to_centroid_dist[sil_df$gene]

# Scale distances for plotting on the same graph
max_dist <- max(sil_df$gene_to_centroid_dist, na.rm = TRUE)
sil_df$scaled_dist <- sil_df$gene_to_centroid_dist / max_dist

setDT(sil_df)

ggclust_stats <- ggplot(sil_df, aes(x = gene, y = sil_width, fill = factor(cluster))) +
  # Bars for silhouette widths
  geom_bar(stat="identity", width=0.8) +
  # Points for scaled distance
  geom_point(aes(y = scaled_dist), color = "red", size = 1.5) +
  facet_wrap(~ cluster, scales = "free_x", ncol=2) +
  # Primary axis for silhouette, secondary axis for distance
  scale_y_continuous(
    name = "Silhouette Width",
    sec.axis = sec_axis(~ . * max_dist, name = "Gene-to-Centroid Distance")
  ) +
  labs(x = "Genes (ordered by silhouette width)",
       title = "Silhouette Width by Cluster with Gene-to-Centroid Distance") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "italic", angle = 90, vjust = 0.5),
    axis.ticks.x = element_blank()
  )


