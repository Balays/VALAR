library(TSSr)
library(BSgenome.PRV.MdBio.1.0)
library(data.table)

outdir <- config$outdir


prime5.counts <- fread(paste0(outdir, '/prime5.counts.tsv'), na.strings = '')

prime5.counts[, sample    := gsub('-', '_', sample)]
prime5.counts[, cell_line := gsub('-', '_', cell_line)]
prime5.counts[, group     := paste(cell_line, hpi, sep='_')]

## filter out C6 0.5 h samples
prime5.counts <- prime5.counts[!grepl('C6_0.5h', sample),]
## and dRNA
prime5.counts <- prime5.counts[!grepl('dRNA', sample),]


TssData <- dcast.data.table(prime5.counts[correct_tss == T], 
                            seqnames + pos + strand ~ sample, 
                            value.var = 'count', fill = 0)

colnames(TssData)[1] <- 'chr'

colSums(TssData[,-c(1:3)])

fwrite(TssData, paste0(outdir, '/TssData.dcDNA.txt'), sep='\t')


sampleLabels  <- data.table(sampleLabels=colnames(TssData)[-c(1:3)])
sampleLabels[, sampleLabelsMerged := gsub('_[1-3]$', '', sampleLabels)]
sampleLabels[, mergeIndex         := .GRP, by=sampleLabelsMerged]

  

#set parameters
refSource <- "refgenome/Refgenome/LT934125.1-2.gff3"
organismName <- "PRV"
directory_path <- paste0(outdir, "/TSSr"); dir.create(directory_path)


#create TSSr object
myTSSr <- new("TSSr", genomeName = "BSgenome.PRV.MdBio.1.0",
              inputFiles = paste0(outdir, '/TssData.dcDNA.txt'),
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
#ggp <- plotCorrelation(myTSSr, samples = "all")
# PCA
plotTssPCA(myTSSr, TSS.threshold=1)

# Merge replicates
TSSr::mergeSamples(myTSSr)

# export merged data
#exportTSStable(myTSSr, data = "raw", merged = "FALSE")

TSS.raw.dt  <- as.data.table(myTSSr@TSSprocessedMatrix)

fwrite(TSS.raw.dt, paste0(outdir, '/TSS.dcDNA.raw.txt'), sep='\t')


#myTSSr@librarySizes
# Filter and normalize counts
filterTSS(myTSSr, method = "poisson", normalization = T, pVal =0.01, tpmLow = 0.1)

#exportTSStable(myTSSr, data = "processed") 

TSS.dt <- myTSSr@TSSprocessedMatrix

fwrite(TSS.dt, paste0(outdir, '/TSS.dcDNA.processed.txt'), sep='\t')


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
tagClusters  <- merge(tagClusters, unique(prime5.counts[,.(group, hpi, cell_line, Time)]), by='group', all.x=T)

fwrite(tagClusters, paste0(outdir, '/TSSr.dcDNA.all.tagClusters.txt'), sep='\t')

#### -->> CHECK THE CLUSTERS ! 
## histogram
ggplot(tagClusters) +
  geom_boxplot(aes(hpi, cluster_width, color=cell_line), binwidth =  10) +
  #xlim(c(0,250))  +
  #facet_nested(rows=vars(hpi)) +
  theme_minimal() +
  theme(legend.position = 'bottom')


# Consensus clusters (Core Promoters)
consensusCluster(myTSSr, dis = 25)

ConsClusters <- rbindlist(myTSSr@consensusClusters, idcol = 'group', use.names = T)
ConsClusters[, cluster_width := end - start +1]
# 
ConsClusters  <- merge(ConsClusters, unique(prime5.counts[,.(group, hpi, cell_line, Time)]), by='group', all.x=T)

fwrite(ConsClusters, paste0(outdir, '/TSSr.dcDNA.all.ConsClusters.txt'), sep='\t')


# Cluster Shape of Consensus Clusters 
shapeCluster(myTSSr,clusters = "consensusClusters", method = "PSS",useMultiCore= FALSE, numCores = NULL)

ClusterShapes <- rbindlist(myTSSr@clusterShape, idcol = 'group', use.names = T)
ClusterShapes[, cluster_width := end - start +1]
# 
ClusterShapes  <- merge(ClusterShapes, unique(prime5.counts[,.(group, hpi, cell_line, Time)]), by='group', all.x=T)

fwrite(ClusterShapes, paste0(outdir, '/TSSr.dcDNA.all.ClusterShapes.txt'), sep='\t')


#
saveRDS(myTSSr, paste0(outdir, '/TSSr.dcDNA.rds'))

# Read back
myTSSr <- readRDS(paste0(outdir, '/TSSr.dcDNA.rds'))

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
fwrite(TSSr_clusters, paste0(outdir, '/TSSr.dcDNA.all.Cluster.Results.txt'), sep='\t')

#
fwrite(TSSr_clusters_uni, paste0(outdir, '/TSSr.dcDNA.uni.Cluster.Results.txt'), sep='\t')



### Peak shifts
ggplot(TSSr_clusters) +
  geom_boxplot(aes(cell_line, peak_shift)) +
  geom_point(aes(cell_line, peak_shift)) + 
  facet_nested(cols=vars(hpi)) +
  theme_minimal()






stop()

# Core promoter shifts
shiftPromoter(myTSSr,comparePairs=list(c("control","treat")), pval = 0.01)

# Analysis of enhancers
callEnhancer(myTSSr, flanking = 250)



##### 
#ConsClusters <- ConsClusters[!grepl('dRNA', group), ]

#### Meta Clusters
filtered_clusters <- ConsClusters[,.(cluster, group, 
                                     #hpi = unlist(stri_extract_last_regex(group, '[0-9\\.]*h')),
                                     #cell_line = stri_replace_first_regex(group, '_.*h', ''),
                                     cluster_start = start, cluster_end = end, seqnames = chr, strand,
                                     cluster_width = end - start + 1,
                                     cluster_peak = dominant_tss,
                                     cluster_total_count = tags )]

# 
filtered_clusters <- merge(filtered_clusters, unique(prime5.counts[,.(group, hpi, cell_line, Time)]), by='group', all.x=T)

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


fwrite(filtered_clusters, paste0(outdir, '/TSSr.filtered_clusters.tsv'), sep='\t')

filtered_clusters <- fread(paste0(outdir, '/TSSr.filtered_clusters.tsv'), na.strings = '')

meta_summary <- filtered_clusters[, .(
  meta_cluster_start = min(cluster_start),
  meta_cluster_end   = max(cluster_end),
  meta_cluster_width = max(cluster_end) - min(cluster_start) +1,
  consensus_peak     = mean(cluster_peak),  # or weighted average if desired
  min_peak = min(cluster_peak),
  max_peak = max(cluster_peak),
  mean_peak = mean(cluster_peak),
  #width_range        = max(cluster_width) - min(cluster_width),
  sample_count       = .N,
  hpi_samples        = paste(unique(hpi), collapse = ", "),
  cell_lines         = paste(unique(cell_line), collapse = ", "),
  total_read_count   = sum(cluster_total_count)
), by = .(seqnames, strand, meta_cluster)]

filtered_clusters <- merge(filtered_clusters, meta_summary, by=c('seqnames', 'strand', 'meta_cluster'))

filtered_clusters[, peak_shift := consensus_peak - cluster_peak, by = .(seqnames, strand, hpi, cell_line, cluster)]


### Cluster and Meta-cluster widths

ggplot(filtered_clusters) +
  geom_boxplot(aes(cell_line, cluster_width)) +
  facet_nested(cols=vars(hpi)) +
  theme_minimal()


ggplot(meta_summary) +
  geom_boxplot(aes(strand, meta_cluster_width))


## histogram
ggplot(tagClusters) +
  geom_histogram(aes(cluster_width, fill=cell_line), binwidth =  10) +
  xlim(c(0,1000))  +
  facet_nested(rows=vars(hpi)) +
  theme_minimal() +
  theme(legend.position = 'bottom')


cowplot::plot_grid(
  ggplot(filtered_clusters) +
    geom_histogram(aes(cluster_width, fill=cell_line), binwidth =  10) +
    xlim(c(0,1000))  +
    facet_nested(rows=vars(hpi)) +
    theme_minimal() +
    theme(legend.position = 'bottom')
  ,
  ggplot(meta_summary) +
    geom_histogram(aes(meta_cluster_width), binwidth =  10) +
    xlim(c(0,1000))  +
    facet_nested(rows=vars(seqnames)) +
    theme_minimal()
  ,
  ncol = 1, rel_heights = c(3,1), axis = 'tblr')

# the cluster and meta_cluster widths do not differ too much.
# so, the cluster refinement - parametric distribution fitting - will be carried out on the meta-clusters.


### Peak shifts
ggplot(filtered_clusters) +
  geom_boxplot(aes(cell_line, peak_shift)) +
  geom_point(aes(cell_line, peak_shift)) + 
  facet_nested(cols=vars(hpi)) +
  theme_minimal()


ggplot(meta_summary) +
  geom_bar(aes(peak_shift))























