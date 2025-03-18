# 3. Cluster Stability via pvclust to assess the statistical support for clusters:
## euclidean distance
#fit.gene <- pvclust(t(countData), method.hclust="complete", method.dist="euclidean", nboot=1000)
#plot(fit.gene)
#pvrect(fit.gene, alpha=0.85, pv="bp")

# Cut the tree where it makes sense biologically
#clusters <- cutree(fit.gene$hclust, k = 5) # for example, 5 clusters

## NOT GOOD:
#method.dist="manhattan"

# method.hclust= all yield very similar results if uncentered distance measure was used
# "average", "complete", "ward.D2"
##
fit.gene  <- pvclust(t(countData),
                     method.hclust = method.hclust,
                     method.dist   = method.dist,
                     nboot=1000)

#plot(fit.gene)
hc <- fit.gene$hclust
#pvrect(fit.gene, alpha=0.85, pv="si")
#rect.hclust(hc, k = cluster_num, border = 2:6) # draws colored rectangles around the 5 clusters

# Cut the tree where it makes sense biologically
#clusters <- cutree(hc, k = cluster_num) # for example, 5 clusters


#cluster_dt <- data.table(cluster=clusters, gene=rownames(countData))

## store the distance matrix as well for the silhoutte
#dist_matrix <- pvclust:::dist.pvclust(t(countData), method=fit.gene$hclust$dist.method)

# Remove genes with zero variance
## correlation distance -->> not good
#fit.gene2 <- pvclust(t( countData[apply(countData, 1, var) > 0, ]),
#                    method.hclust="complete", method.dist="correlation", nboot=1000)
#plot(fit.gene2)
#pvrect(fit.gene2)
