library(TSSr)
library(BSgenome.PRV.MdBio.1.0)
library(data.table)

prime5.counts <- fread(paste0(outdir, '/prime5.counts.tsv'), na.strings = '')

prime5.counts[, sample := gsub('-', '_', sample)]

prime5.counts <- prime5.counts[!grepl('C6_0.5h', sample),]

TssData <- dcast.data.table(prime5.counts[correct_tss == T], 
                            seqnames + pos + strand ~ sample, 
                            value.var = 'count', fill = 0)

colnames(TssData)[1] <- 'chr'

colSums(TssData[,-c(1:3)])

fwrite(TssData, paste0(outdir, '/TssData.txt'), sep='\t')


sampleLabels  <- data.table(sampleLabels=colnames(TssData)[-c(1:3)])
sampleLabels[, sampleLabelsMerged := gsub('_[1-3]$', '', sampleLabels)]
sampleLabels[, mergeIndex         := .GRP, by=sampleLabelsMerged]

  

#set parameters
refSource <- "refgenome/Refgenome/LT934125.1-2.gff3"
organismName <- "PRV"
directory_path <- paste0(outdir, "/TSSr/")


#create TSSr object
myTSSr <- new("TSSr", genomeName = "BSgenome.PRV.MdBio.1.0",
              inputFiles = paste0(outdir, '/TssData.txt'),
              inputFilesType = 'TSStable',
              sampleLabels = sampleLabels$sampleLabels,
              sampleLabelsMerged = unique(sampleLabels$sampleLabelsMerged),
              mergeIndex = sampleLabels$mergeIndex,
              refSource = refSource,
              organismName = organismName
)


getTSS(myTSSr)#, sequencingQualityThreshold = 10, mappingQualityThreshold = 20)

# export raw data
exportTSStable(myTSSr, data = "raw", merged = "FALSE")

#create correlation plot
#ggp <- plotCorrelation(myTSSr, samples = "all")
# PCA
plotTssPCA(myTSSr, TSS.threshold=1)

# Merge replicates
TSSr::mergeSamples(myTSSr)

# export merged data
exportTSStable(myTSSr, data = "raw", merged = "FALSE")

TSS.raw.dt  <- as.data.table(myTSSr@TSSprocessedMatrix)


#myTSSr@librarySizes
# Filter and normalize counts
filterTSS(myTSSr, method = "poisson", normalization = T, pVal =0.01, tpmLow = 0.1)

exportTSStable(myTSSr, data = "processed") 

TSS.dt <- myTSSr@TSSprocessedMatrix

# Cluster into TSS Clusters
clusterTSS(myTSSr, method = "peakclu",peakDistance=100,extensionDistance=30
           ,localThreshold = 0.02, clusterThreshold = 1
           ,useMultiCore=FALSE, numCores=NULL)


# Extract clusters
tagClusters <- rbindlist(myTSSr@tagClusters, idcol = 'group', use.names = T)
fwrite(tagClusters, 'TSSr.dcDNA.all.tagClusters.txt', sep='\t')


# Consensus clusters
consensusCluster(myTSSr)

# Extract Consensus clusters
ConsClusters <- rbindlist(myTSSr@consensusClusters, idcol = 'group', use.names = T)
fwrite(ConsClusters, 'TSSr.dcDNA.all.ConsClusters.txt', sep='\t')


#exportTSStoBedgraph(myTSSr, data = "processed", format = "bedGraph") 
#exportTSStoBedgraph(myTSSr, data = "processed", format = "BigWig")

#exportClustersTable(myTSSr, data = "tagClusters")
#exportClustersToBed(myTSSr, data = "tagClusters") 

#exportClustersTable(myTSSr, data = "consensusClusters")


tagClusters <- fread('TSSr.dcDNA.all.tagClusters.txt')

ConsClusters <- fread('TSSr.dcDNA.all.ConsClusters.txt')



##### 
ConsClusters <- ConsClusters[!grepl('dRNA', group), ]

#### Meta Clusters
filtered_clusters <- ConsClusters[,.(hpi = unlist(stri_extract_first_regex('C6_12h', '[0-9]*h')),
                                     cell_line = gsub('_[0-9]*h', '', group),
                                     cluster_start = start, cluster_end = end, seqnames = chr, strand,
                                     cluster_width = end - start + 1,
                                     cluster_peak = dominant_tss,
                                     cluster_total_count = tags )]

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

ggplot(meta_summary) +
  geom_bar(aes(peak_shift))























