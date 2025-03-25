
library(TSSr)
library(BSgenome.PRV.MdBio.1.0)
library(data.table)

outdir <- config$outdir


prime3.counts <- fread(paste0(outdir, '/prime3.counts.tsv'), na.strings = '')

prime3.counts[, sample    := gsub('-', '_', sample)]
prime3.counts[, cell_line := gsub('-', '_', cell_line)]
prime3.counts[, group     := paste(cell_line, libtype, sep='_')]

## filter out C6 0.5 h samples
#prime3.counts <- prime3.counts[!grepl('C6_0.5h', sample),]
## but not dRNA, right?
prime3.counts <- prime3.counts[grepl('dRNA', sample),]


prime3.counts[, hpi  := libtype      ]
prime3.counts[, Time := libtype      ]

## Filtering of very-low coverage samples
samples_to_filt <- c('C6_2h_3', 'PC_12_8h_3', 'C6_0.5h_1', 'C6_0.5h_2')

prime3.counts <- prime3.counts[!prime3.counts$sample %in% samples_to_filt, ]

message(paste(samples_to_filt, ' was filtered out!', collapse = '\n') )


##
TesData <- dcast.data.table(prime3.counts[correct_tes == T], 
                            seqnames + pos + strand ~ sample, 
                            value.var = 'count', fill = 0)

colnames(TesData)[1] <- 'chr'

colSums(TesData[,-c(1:3)])



##
fwrite(TesData, paste0(outdir, '/TesData.dcDNA.txt'), sep='\t')


## MetaData 
sampleLabels  <- data.table(sampleLabels=colnames(TesData)[-c(1:3)])
sampleLabels[, sampleLabelsMerged := gsub('_[1-3]$', '', sampleLabels)]
sampleLabels[, mergeIndex         := .GRP, by=sampleLabelsMerged]

  

#set parameters
refSource <- "refgenome/LT934125.1.CDS.gff3"
organismName <- "PRV"
directory_path <- paste0(outdir, "/TSSr"); dir.create(directory_path)
#myTSSr@refSource <- refSource

#create TSSr object
myTSSr <- new("TSSr", genomeName = "BSgenome.PRV.MdBio.1.0",
              inputFiles = paste0(outdir, '/TesData.dcDNA.txt'),
              inputFilesType = 'TSStable',
              sampleLabels = sampleLabels$sampleLabels,
              sampleLabelsMerged = unique(sampleLabels$sampleLabelsMerged),
              mergeIndex = sampleLabels$mergeIndex,
              #refSource = refSource,
              refTable = data.table(genes.dt),
              organismName = organismName
)

myTSSr@refTable <- data.table(genes.dt)


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

TES.raw.dt  <- as.data.table(myTSSr@TSSprocessedMatrix)

fwrite(TES.raw.dt, paste0(outdir, '/TSSr.TES.dRNA.raw.txt'), sep='\t')


#myTSSr@librarySizes
# Filter and normalize counts
filterTSS(myTSSr, method = "poisson", normalization = T, pVal =0.01, tpmLow = 0.1)

#exportTSStable(myTSSr, data = "processed") 

TES.dt <- myTSSr@TSSprocessedMatrix

fwrite(TES.dt, paste0(outdir, '/TSSr.TES.dRNA.processed.txt'), sep='\t')


# Viral-Specific TES Clustering
clusterTSS(myTSSr, method = "peakclu", # Same logic applies for clustering peaks.
           peakDistance = 40, # Keep moderate distance since TES peaks are not always as sharp.
           extensionDistance = 15, # TES regions can be broader, allow a bit more extension.
           localThreshold = 0.01, # 	0.01 or lower	TES signals are often weaker, need to keep more low-signal regions.
           clusterThreshold = 1, # TES clusters might not have as strong cumulative signal as TSSs. 
           #useMultiCore = TRUE, numCores = 4
           )

# Extract clusters
tagClusters <- rbindlist(myTSSr@tagClusters, idcol = 'group', use.names = T)
tagClusters[, cluster_width := end - start +1]
# 
tagClusters  <- merge(tagClusters, unique(prime3.counts[,.(group, hpi, cell_line, Time)]), by='group', all.x=T)

fwrite(tagClusters, paste0(outdir, '/TSSr.TES.dRNA.all.tagClusters.txt'), sep='\t')

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
ConsClusters  <- merge(ConsClusters, unique(prime3.counts[,.(group, hpi, cell_line, Time)]), by='group', all.x=T)

fwrite(ConsClusters, paste0(outdir, '/TSSr.TES.dRNA.all.ConsClusters.txt'), sep='\t')


# Cluster Shape of Consensus Clusters 
shapeCluster(myTSSr,clusters = "consensusClusters", method = "PSS", useMultiCore= FALSE, numCores = NULL)

ClusterShapes <- rbindlist(myTSSr@clusterShape, idcol = 'group', use.names = T)
ClusterShapes[, cluster_width := end - start +1]
# 
ClusterShapes  <- merge(ClusterShapes, unique(prime3.counts[,.(group, hpi, cell_line, Time)]), by='group', all.x=T)

fwrite(ClusterShapes, paste0(outdir, '/TSSr.TES.dRNA.all.ClusterShapes.txt'), sep='\t')


####
# Annotation (Assigning TCs to genes)
annotateCluster(  myTSSr
                  , clusters = "consensusClusters"
                  , filterCluster = TRUE
                  , filterClusterThreshold = 0.02
                  , annotationType = "genes"
                  , upstream=50
                  , upstreamOverlap = 50
                  , downstream = 1000)


assignedClusters <- rbindlist(myTSSr@assignedClusters, idcol = 'group', use.names = T)
assignedClusters[, cluster_width := end - start +1]
# 
assignedClusters  <- merge(assignedClusters, unique(prime3.counts[,.(group, hpi, cell_line, Time)]), by='group', all.x=T)

# get the dominant TC per gene (highest tag)
assignedClusters[,  dominant_TC := max(tags), by=.(gene, group)]
assignedClusters[,  dominant_TC := fifelse(dominant_TC == tags,  cluster, 0)]

## still there are more TC for some genes, so we need to check the distance of the TC to gene and select the closer one
CDS.dt <- data.table(viral.CDS)[,.(CDS_start = min(start), CDS_end = max(end)), by=.(seqnames, strand, gene)]

## consider spliced CDSs!
CDS.dt[, CDS_prime5 := fifelse(strand == '+', CDS_start, CDS_end)]
CDS.dt[, CDS_prime3 := fifelse(strand == '-', CDS_start, CDS_end)]

dup(CDS.dt$gene)
assignedClusters <- merge(assignedClusters, CDS.dt, by.x=c('chr', 'strand', 'gene'), by.y=c('seqnames', 'strand', 'gene'), all.x = T)

assignedClusters[,  TSS_CDS_dist := CDS_prime3 - dominant_tss]

# this on a positve strand should be negative values, on the negative strand it should be a positive value!

ggplot(assignedClusters, aes(TSS_CDS_dist)) + geom_histogram() + facet_grid(cols=vars(strand))
## check! filter?


assignedClusters[,  closest_TC  := min(abs(TSS_CDS_dist)), by=.(gene, group)]
assignedClusters[,  closest_TC  := fifelse(closest_TC == abs(TSS_CDS_dist),  cluster, 0)]

##
fwrite(assignedClusters, paste0(outdir, '/TSSr.TES.dRNA.all.assignedClusters.txt'), sep='\t')


## select closest and most dominant cluster for each gene
assignedClusters  <- assignedClusters[ closest_TC  == cluster]

assignedClusters  <- assignedClusters[ dominant_TC > 0]

setnames(assignedClusters, colnames(assignedClusters), gsub('tss', 'tes', colnames(assignedClusters)) )

##
assignedClusters.sp <- dcast(assignedClusters[dominant_TC > 0], gene ~ group, value.var = 'dominant_tes')



## Gene-wise TES shifts from mean dominant_TES
assignedClusters[, gene_mean_dominant_TES  := round(mean(dominant_tes), 0), by=.(gene)]

assignedClusters[, gene_dominant_TES_shift := gene_mean_dominant_TES - dominant_tes, by = .(cell_line, hpi, group, Time, gene)]
assignedClusters[, gene_mean_shift_TES     := mean(abs(gene_dominant_TES_shift)), by=.(gene)]
assignedClusters[, gene_sd_shift_TES       := sd(abs(gene_dominant_TES_shift)),   by=.(gene)]
assignedClusters[, gene_varcoeff_shif_TES  := gene_sd_shift_TES / gene_mean_shift_TES,   by=.(gene)]

gene_TES_shifts_Sum <- unique(assignedClusters[, .(chr, strand, gene, gene_mean_shift_TES, gene_sd_shift_TES, gene_varcoeff_shif_TES)])
gene_TES_shifts_Sum <- gene_TES_shifts_Sum[order(gene_mean_shift_TES, decreasing = T)]
topShifts   <- gene_TES_shifts_Sum[1:20, gene]
topShifts   <- topShifts[!topShifts %in% c('US1_2', 'IE180_2')]


## Plot Gene-Wise Shifts
ggplot(assignedClusters[gene %in% topShifts
                      #& Time %in% c(1,2,4,6,8,12)  
                      ], 
       aes(gene_dominant_TES_shift, cell_line)) + 
  geom_line(aes(color = cell_line, group = cell_line)) + 
  geom_point(aes(color = cell_line)) + 
  #coord_flip() +
  facet_nested_wrap(~gene, ncol=4) +
  theme_minimal()


##
fwrite(assignedClusters, paste0(outdir, '/TSSr.TES.dRNA.best.assignedClusters.txt'), sep='\t')


# Analysis of enhancers
#callEnhancer(myTSSr, flanking = 200, dis2gene = 500)
flanking = 200
dis2gene = 500
object <- myTSSr
source('TSSr.Enhancers.R')
myTSSr <- object
#
# TSSr can identify this bidirectional cluster pairs and calculate a 
# sample-set wide directionality score D for each locus (Andersson et al., 2014). 
# D = (F-R)/(F+R), where F is the aggregated normalized tag counts in forward strand (p) and 
# R is the aggregated normalized tag counts in reverse strand (m). 
# Putative enhancers were then filtered with |D| < 0.8.

enhancers.dt <- rbindlist(myTSSr@enhancers, idcol = 'group', use.names = T)
#enhancers.dt[, cluster_width := end - start +1]

# 
enhancers.dt  <- merge(enhancers.dt, unique(prime3.counts[,.(group, hpi, cell_line, Time)]), by='group', all.x=T)


fwrite(enhancers.dt, paste0(outdir, '/TSSr.TES.dRNA.enhancers.txt'), sep='\t')


##### END ANALYSIS

#
saveRDS(myTSSr, paste0(outdir, '/TSSr.TES.dRNA.rds'))

# Read back
myTSSr <- readRDS(paste0(outdir, '/TSSr.TES.dRNA.rds'))

TSS.raw.dt    <- fread(paste0(outdir, '/TSSr.TES.dRNA.raw.txt'))

TSS.dt        <- fread(paste0(outdir, '/TSSr.TES.dRNA.processed.txt'))

tagClusters   <- fread(paste0(outdir, '/TSSr.TES.dRNA.all.tagClusters.txt'))

ConsClusters  <- fread(paste0(outdir, '/TSSr.TES.dRNA.all.ConsClusters.txt'))

ClusterShapes <- fread(paste0(outdir, '/TSSr.TES.dRNA.all.ClusterShapes.txt'))

##
assignedClusters <- fread(paste0(outdir, '/TSSr.TES.dRNA.best.assignedClusters.txt'))


### Combine Clusters Across Samples

TSSr_clusters <- copy(ClusterShapes)

# TC Summaries
TC_summary <- TSSr_clusters[, .(
  TC.start = min(start), 
  TC.end   = max(end),   
  consensus_peak   = round(mean(dominant_tss),0),  # or weighted average if desired
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

# Merge with Gene Annotation -> Dont need to
# TSSr_clusters <- merge(TSSr_clusters, assignedClusters, by=cols_by, all = T)
## OK!

# Peak shifts
TSSr_clusters[, peak_shift := consensus_peak - dominant_tss, by = .(chr , strand, hpi, cell_line, cluster)]
TSSr_clusters[, TC.width   := max_peak - min_peak,  by = .(cluster)]


# unique clusters
TSSr_clusters_uni <- unique(TSSr_clusters[,.(seqnames = chr, strand, cluster, TC.start, TC.end, TC.width, support, score, 
                                             mean_shape, min_shape, max_shape,
                                             consensus_peak, min_peak, max_peak, 
                                             hpi_samples, cell_lines)])

#
fwrite(TSSr_clusters, paste0(outdir, '/TSSr.TES.dRNA.all.Cluster.Results.txt'), sep='\t')

#
fwrite(TSSr_clusters_uni, paste0(outdir, '/TSSr.TES.dRNA.uni.Cluster.Results.txt'), sep='\t')



### Peak shifts
ggplot(TSSr_clusters, aes(cell_line, peak_shift, colour = cell_line)) +
  geom_boxplot() +
  geom_point() + 
  theme_minimal() +
  facet_nested(cols=vars(Time), rows =vars(strand))
 


stop()



## Plot Enhancers
gg <- ggplot(enhancers.dt, aes(cell_line, D)) + 
  geom_point(aes(fill = cell_line), shape = 21) + 
  coord_flip() +
  theme_minimal() +
  facet_nested(rows=vars(Time), scales='free', independent = F)


ggsave(paste0(outdir, '/TSSr.enhancers.D.jpg'), width = 20, height = 24)


enhancers.stranded <- melt(enhancers.dt, id.vars = c(1:7, 10:13) )
enhancers.stranded[variable == 'tags.m', value := value * -1]

gg <- ggplot(enhancers.stranded, aes(enhancer, value)) + 
  geom_col(aes(fill = variable), position = 'identity') + 
  coord_flip() +
  theme_minimal() +
  facet_nested(rows=vars(cell_line), cols=vars(Time), scales='free', independent = F)


ggsave(paste0(outdir, '/TSSr.enhancers.jpg'), width = 20, height = 24)



# Core promoter shifts
shiftPromoter(myTSSr,comparePairs=list(c("control","treat")), pval = 0.01)


# TSS plots - doesent work
plotTSS(myTSSr
        , samples=c("PC_12_6h", "C6_6h")
        , tssData = "processed"
        , clusters = "assigned"
        , clusterThreshold = 0.02 
        , genelist=c("US1")
        , up.dis =500, down.dis = 100, yFixed=TRUE)


##### 

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























