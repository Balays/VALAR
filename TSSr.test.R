


#### Testing TSS parameters


# Viral-Specific TSS Clustering
clusterTSS(myTSSr, method = "peakclu", 
           peakDistance = 50,      # Reduce to enforce closer peaks
           extensionDistance = 15, # Smaller extensions per cluster
           localThreshold = 0.05,  # Higher threshold for dominant signals
           clusterThreshold = 2   # Require stronger clusters
           #, useMultiCore = TRUE
           #, numCores = 14  # Enable parallelization
)


# Extract clusters
tagClusters <- rbindlist(myTSSr@tagClusters, idcol = 'group', use.names = T)
tagClusters[, cluster_width := end - start +1]
# 
tagClusters  <- merge(tagClusters, unique(prime5.counts[,.(group, hpi, cell_line, Time)]), by='group', all.x=T)

ggplot(tagClusters) +
  geom_boxplot(aes(cell_line, cluster_width)) +
  theme_minimal()


tagClusters.p1 <- tagClusters



#exportTSStoBedgraph(myTSSr, data = "processed", format = "bedGraph") 
#exportTSStoBedgraph(myTSSr, data = "processed", format = "BigWig")

#exportClustersTable(myTSSr, data = "tagClusters")
#exportClustersToBed(myTSSr, data = "tagClusters") 

#exportClustersTable(myTSSr, data = "consensusClusters")
