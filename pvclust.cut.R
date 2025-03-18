


plot(fit.gene)
hc <- fit.gene$hclust
#pvrect(fit.gene, alpha=0.85, pv="si")
rect.hclust(hc, k = cluster_num, border = 2:6) # draws colored rectangles around the 5 clusters

# Cut the tree where it makes sense biologically
clusters <- cutree(hc, k = cluster_num) # for example, 5 clusters


cluster_dt <- data.table(cluster=clusters, gene=rownames(countData))

## store the distance matrix as well for the silhoutte
dist_matrix <- pvclust:::dist.pvclust(t(countData), method=fit.gene$hclust$dist.method)
