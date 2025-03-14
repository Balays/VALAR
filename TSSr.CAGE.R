library(TSSr)
library(BSgenome.PRV.MdBio.1.0)
library(data.table)

outdir <- CAGE_config$outdir


prime.counts.CAGE <- fread(paste0(outdir, '/prime.counts.CAGE.tsv'))


prime5.counts <- prime.counts.CAGE[endtype == 'prime5']

prime5.counts[, sample    := gsub('-', '_', sample)]
prime5.counts[, cell_line := gsub('-', '_', cell_line)]
prime5.counts[, group     := paste(cell_line, hpi, sep='_')]

## filter out C6 0.5 h samples
prime5.counts <- prime5.counts[!grepl('C6_0.5h', sample),]
## and dRNA
prime5.counts <- prime5.counts[!grepl('dRNA', sample),]


TssData <- dcast.data.table(prime5.counts, #[correct_tss == T], 
                            seqnames + pos + strand ~ sample, 
                            value.var = 'count', fill = 0)

colnames(TssData)[1] <- 'chr'

colSums(TssData[,-c(1:3)])

fwrite(TssData, paste0(outdir, '/TssData.CAGE.txt'), sep='\t')


sampleLabels  <- data.table(sampleLabels=colnames(TssData)[-c(1:3)])
sampleLabels[, sampleLabelsMerged := gsub('_[1-3]$', '', sampleLabels)]
sampleLabels[, mergeIndex         := .GRP, by=sampleLabelsMerged]



#set parameters
refSource <- "refgenome/Refgenome/LT934125.1-2.gff3"
organismName <- "PRV"
directory_path <- paste0(outdir, "/TSSr"); dir.create(directory_path)


#create TSSr object
myTSSr <- new("TSSr", genomeName = "BSgenome.PRV.MdBio.1.0",
              inputFiles = paste0(outdir, '/TssData.CAGE.txt'),
              inputFilesType = 'TSStable',
              sampleLabels = sampleLabels$sampleLabels,
              sampleLabelsMerged = unique(sampleLabels$sampleLabelsMerged),
              mergeIndex = sampleLabels$mergeIndex,
              refSource = refSource,
              organismName = organismName
)


getTSS(myTSSr)#, sequencingQualityThreshold = 10, mappingQualityThreshold = 20)

# export raw data
#exportTSStable(myTSSr, data = "raw", merged = "FALSE")

#create correlation plot
plotCorrelation(myTSSr, samples = "all")
# PCA
plotTssPCA(myTSSr, TSS.threshold=1)

# Merge replicates
TSSr::mergeSamples(myTSSr)

# export merged data
#exportTSStable(myTSSr, data = "raw", merged = "FALSE")

TSS.raw.dt  <- as.data.table(myTSSr@TSSprocessedMatrix)

fwrite(TSS.raw.dt, paste0(outdir, '/TSSr.CAGE.raw.txt'), sep='\t')


#myTSSr@librarySizes
# Filter and normalize counts
filterTSS(myTSSr, method = "poisson", normalization = T, pVal =0.01, tpmLow = 0.1)

#exportTSStable(myTSSr, data = "processed") 

TSS.dt <- myTSSr@TSSprocessedMatrix

fwrite(TSS.dt, paste0(outdir, '/TSSr.CAGE.processed.txt'), sep='\t')


# Viral-Specific TSS Clustering
clusterTSS(myTSSr, method = "peakclu", 
           peakDistance = 40,      # Reduce to enforce closer peaks
           extensionDistance = 10, # Smaller extensions per cluster
           localThreshold = 0.05,  # Higher threshold for dominant signals
           clusterThreshold = 2   # Require stronger clusters
           #, useMultiCore = TRUE
           #, numCores = 14  # Enable parallelization
)

# Extract clusters
tagClusters <- rbindlist(myTSSr@tagClusters, idcol = 'group', use.names = T)
tagClusters[, cluster_width := end - start +1]
# 
tagClusters  <- merge(tagClusters, unique(prime5.counts[,.(sample, group, hpi, cell_line, Time)]), by.x='group', by.y = 'sample', all.x=T)

fwrite(tagClusters, paste0(outdir, '/TSSr.CAGE.all.tagClusters.txt'), sep='\t')


#### -->> CHECK THE CLUSTERS ! 
## histogram
ggplot(tagClusters) +
  geom_histogram(aes(cluster_width, fill=cell_line), binwidth =  10) +
  #xlim(c(0,250))  +
  facet_nested(rows=vars(hpi)) +
  theme_minimal() +
  theme(legend.position = 'bottom')


# Consensus clusters (Core Promoters)
consensusCluster(myTSSr, dis = 25)
#
ConsClusters <- rbindlist(myTSSr@consensusClusters, idcol = 'group', use.names = T)
ConsClusters[, cluster_width := end - start +1]
# 
ConsClusters  <- merge(ConsClusters, unique(prime5.counts[,.(sample, group, hpi, cell_line, Time)]), by.x='group', by.y = 'sample', all.x=T)

fwrite(ConsClusters, paste0(outdir, '/TSSr.CAGE.all.ConsClusters.txt'), sep='\t')


# Cluster Shape of Consensus Clusters 
shapeCluster(myTSSr,clusters = "consensusClusters", method = "PSS",useMultiCore= FALSE, numCores = NULL)
#
ClusterShapes <- rbindlist(myTSSr@clusterShape, idcol = 'group', use.names = T)
ClusterShapes[, cluster_width := end - start +1]
# 
ClusterShapes  <- merge(ClusterShapes, unique(prime5.counts[,.(sample, group, hpi, cell_line, Time)]), by.x='group', by.y = 'sample', all.x=T)

fwrite(ClusterShapes, paste0(outdir, '/TSSr.CAGE.all.ClusterShapes.txt'), sep='\t')


#
saveRDS(myTSSr, paste0(CAGE_config$outdir, '/TSSr.CAGE.rds'))


# Read back
myTSSr <- readRDS(paste0(CAGE_config$outdir, '/TSSr.CAGE.rds'))

TSS.raw.dt    <- fread(paste0(outdir, '/TSSr.CAGE.raw.txt'))

TSS.dt        <- fread(paste0(outdir, '/TSSr.CAGE.processed.txt'))

tagClusters   <- fread(paste0(outdir, '/TSSr.CAGE.all.tagClusters.txt'))

ConsClusters  <- fread(paste0(outdir, '/TSSr.CAGE.all.ConsClusters.txt'))

ClusterShapes <- fread(paste0(outdir, '/TSSr.CAGE.all.ClusterShapes.txt'))



### Combine Clusters Across Samples

TSSr_clusters <- copy(ClusterShapes)

# TC Summaries
TC_summary <- TSSr_clusters[, .(
  TC.start = min(start), 
  TC.end   = max(end),   
  consensus_peak   = mean(dominant_tss),  # or weighted average if desired
  min_peak    = min(dominant_tss),
  max_peak    = max(dominant_tss),
  mean_shape  = min(shape.score),
  min_shape   = min(shape.score),
  max_shape   = max(shape.score),
  support     = uniqueN(group), # support is the number of sample (groups) in which the cluster was found
  score       = sum(tags), # score is the sum of the tags of the cluster in all the groups
  hpi_samples        = paste(unique(hpi), collapse = ", "),
  cell_lines         = paste(unique(cell_line), collapse = ", ")
), by = .(cluster)]

TSSr_clusters <- merge(TSSr_clusters, TC_summary, by=c('cluster'))

# Peak shifts
TSSr_clusters[, peak_shift := consensus_peak - dominant_tss, by = .(chr , strand, hpi, cell_line, cluster)]
TSSr_clusters[, TC.width   := max_peak - min_peak,  by = .(cluster)]


# unique clusters
TSSr_clusters_uni <- unique(TSSr_clusters[,.(seqnames = chr, strand, cluster, TC.start, TC.end, TC.width, support, score, 
                                             mean_shape, min_shape, max_shape,
                                             consensus_peak, min_peak, max_peak, 
                                             hpi_samples, cell_lines)])


#
fwrite(TSSr_clusters, paste0(outdir, '/TSSr.CAGE.all.Cluster.Results.txt'), sep='\t')
#
fwrite(TSSr_clusters_uni, paste0(outdir, '/TSSr.CAGE.uni.Cluster.Results.txt'), sep='\t')



### Peak shifts
ggplot(TSSr_clusters) +
  geom_boxplot(aes(cell_line, peak_shift)) +
  geom_point(aes(cell_line, peak_shift)) + 
  facet_nested(cols=vars(hpi)) +
  theme_minimal()




stop()
