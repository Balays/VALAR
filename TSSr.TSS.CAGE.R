library(TSSr)
library(BSgenome.PRV.MdBio.1.0)
library(data.table)

outdir <- config$outdir


prime.counts.CAGE <- fread(paste0(CAGE_config$outdir, '/prime.counts.CAGE.tsv'))


prime5.counts <- prime.counts.CAGE[endtype == 'prime5']

prime5.counts[, cell_line := gsub('-', '_', cell_line)]
prime5.counts[, group     := paste(cell_line, hpi, sep='_')]
prime5.counts[, sample    := group]


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
refSource <- "refgenome/LT934125.1.CDS.gff3"
organismName <- "PRV"
directory_path <- paste0(outdir, "/TSSr"); dir.create(directory_path)
#myTSSr@refSource <- refSource

#create TSSr object
myTSSr <- new("TSSr", genomeName = "BSgenome.PRV.MdBio.1.0",
              inputFiles = paste0(outdir, '/TssData.CAGE.txt'),
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

TSS.raw.dt  <- as.data.table(myTSSr@TSSprocessedMatrix)

fwrite(TSS.raw.dt, paste0(outdir, '/TSSr.TSS.CAGE.raw.txt'), sep='\t')

#myTSSr@librarySizes


# Filter and normalize counts
filterTSS(myTSSr, method = CAGE_config$method, normalization = CAGE_config$normalization, pVal =  CAGE_config$pVal, tpmLow =  CAGE_config$tpmLow)

#exportTSStable(myTSSr, data = "processed") 

TSS.dt <- myTSSr@TSSprocessedMatrix

fwrite(TSS.dt, paste0(outdir, '/TSSr.TSS.CAGE.processed.txt'), sep='\t')


# Viral-Specific TSS Clustering
clusterTSS(myTSSr, method = "peakclu"
           , peakDistance = CAGE_config$peakDistance             # 40,    # Reduce to enforce closer peaks
           , extensionDistance = CAGE_config$extensionDistance   # 10,    # Smaller extensions per cluster
           , localThreshold = CAGE_config$localThreshold         # 0.05,  # Higher threshold for dominant signals
           , clusterThreshold = CAGE_config$clusterThreshold     # 2      # Require stronger clusters
           , useMultiCore = { if(.Platform$OS.type=="windows") F else T }
           , numCores = { if(.Platform$OS.type=="windows") 1 else nproc }  # Enable parallelization
)

# Extract clusters
tagClusters <- rbindlist(myTSSr@tagClusters, idcol = 'group', use.names = T)
tagClusters[, cluster_width := end - start +1]
# 
tagClusters  <- merge(tagClusters, unique(prime5.counts[,.(group, hpi, cell_line, Time)]), by='group', all.x=T)

fwrite(tagClusters, paste0(outdir, '/TSSr.TSS.CAGE.all.tagClusters.txt'), sep='\t')

#### -->> CHECK THE CLUSTERS ! 
## histogram
ggplot(tagClusters) +
  geom_boxplot(aes(hpi, cluster_width, color=cell_line)) +
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

fwrite(ConsClusters, paste0(outdir, '/TSSr.TSS.CAGE.all.ConsClusters.txt'), sep='\t')


# Cluster Shape of Consensus Clusters 
shapeCluster(myTSSr,clusters = "consensusClusters", method = "PSS",useMultiCore= FALSE, numCores = NULL)

ClusterShapes <- rbindlist(myTSSr@clusterShape, idcol = 'group', use.names = T)
ClusterShapes[, cluster_width := end - start +1]
# 
ClusterShapes  <- merge(ClusterShapes, unique(prime5.counts[,.(group, hpi, cell_line, Time)]), by='group', all.x=T)

fwrite(ClusterShapes, paste0(outdir, '/TSSr.TSS.CAGE.all.ClusterShapes.txt'), sep='\t')


####
# Annotation (Assigning TCs to genes)
annotateCluster(  myTSSr
                  , clusters = "consensusClusters"
                  , filterCluster = CAGE_config$filterCluster
                  , filterClusterThreshold = CAGE_config$filterClusterThreshold
                  , annotationType = "genes"
                  , upstream=CAGE_config$upstream
                  , upstreamOverlap = CAGE_config$upstreamOverlap
                  , downstream = CAGE_config$downstream)


assignedClusters <- rbindlist(myTSSr@assignedClusters, idcol = 'group', use.names = T)
assignedClusters[, cluster_width := end - start +1]
# 
assignedClusters  <- merge(assignedClusters, unique(prime5.counts[,.(group, hpi, cell_line, Time)]), by='group', all.x=T)

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

assignedClusters[,  TSS_CDS_dist := CDS_prime5 - dominant_tss]

ggplot(assignedClusters, aes(TSS_CDS_dist)) + geom_histogram() + facet_grid(cols=vars(strand))
## OK!

assignedClusters[,  closest_TC  := min(abs(TSS_CDS_dist)), by=.(gene, group)]
assignedClusters[,  closest_TC  := fifelse(closest_TC == abs(TSS_CDS_dist),  cluster, 0)]

##
fwrite(assignedClusters, paste0(outdir, '/TSSr.TSS.CAGE.all.assignedClusters.txt'), sep='\t')


## select closest and most dominant cluster for each gene
assignedClusters  <- assignedClusters[ closest_TC  == cluster]

assignedClusters  <- assignedClusters[ dominant_TC > 0]

##
assignedClusters.sp <- dcast(assignedClusters[dominant_TC > 0], gene ~ group, value.var = 'dominant_tss')



## Gene-wise TSS shifts from mean dominant_TSS
assignedClusters[, gene_mean_dominant_TSS  := round(mean(dominant_tss), 0), by=.(gene)]

assignedClusters[, gene_dominant_TSS_shift := gene_mean_dominant_TSS - dominant_tss, by = .(cell_line, hpi, group, Time, gene)]
assignedClusters[, gene_mean_shift_TSS     := mean(abs(gene_dominant_TSS_shift)), by=.(gene)]
assignedClusters[, gene_sd_shift_TSS       := sd(abs(gene_dominant_TSS_shift)),   by=.(gene)]
assignedClusters[, gene_varcoeff_shif_TSS  := gene_sd_shift_TSS / gene_mean_shift_TSS,   by=.(gene)]

gene_TSS_shifts_Sum <- unique(assignedClusters[, .(chr, strand, gene, gene_mean_shift_TSS, gene_sd_shift_TSS, gene_varcoeff_shif_TSS)])
gene_TSS_shifts_Sum <- gene_TSS_shifts_Sum[order(gene_mean_shift_TSS, decreasing = T)]
topShifts   <- gene_TSS_shifts_Sum[1:20, gene]
topShifts   <- topShifts[!topShifts %in% c('US1_2', 'IE180_2')]


## Plot Gene-Wise Shifts
ggplot(assignedClusters[gene %in% topShifts
                      & Time %in% c(1,2,4,6,8,12)  ], 
       aes(Time, gene_dominant_TSS_shift)) + 
  geom_line(aes(color = cell_line, group = cell_line)) + 
  geom_point(aes(color = cell_line)) + 
  #coord_flip() +
  facet_nested_wrap(~gene, ncol=2) +
  theme_minimal()


##
fwrite(assignedClusters, paste0(outdir, '/TSSr.TSS.CAGE.best.assignedClusters.txt'), sep='\t')


# Analysis of enhancers
#callEnhancer(myTSSr, flanking = 200, dis2gene = 500)
flanking = CAGE_config$flanking
dis2gene = CAGE_config$dis2gene
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
enhancers.dt  <- merge(enhancers.dt, unique(prime5.counts[,.(group, hpi, cell_line, Time)]), by='group', all.x=T)


fwrite(enhancers.dt, paste0(outdir, '/TSSr.TSS.CAGE.enhancers.txt'), sep='\t')


##### END ANALYSIS

#
saveRDS(myTSSr, paste0(outdir, '/TSSr.TSS.CAGE.rds'))

# Read back
myTSSr <- readRDS(paste0(outdir, '/TSSr.TSS.CAGE.rds'))

TSS.raw.dt    <- fread(paste0(outdir, '/TSSr.TSS.CAGE.raw.txt'))

TSS.dt        <- fread(paste0(outdir, '/TSSr.TSS.CAGE.processed.txt'))

tagClusters   <- fread(paste0(outdir, '/TSSr.TSS.CAGE.all.tagClusters.txt'))

ConsClusters  <- fread(paste0(outdir, '/TSSr.TSS.CAGE.all.ConsClusters.txt'))

ClusterShapes <- fread(paste0(outdir, '/TSSr.TSS.CAGE.all.ClusterShapes.txt'))

##
assignedClusters <- fread(paste0(outdir, '/TSSr.TSS.CAGE.best.assignedClusters.txt'))


### Combine Clusters Across Samples

TSSr_clusters <- copy(ClusterShapes)

# TC Summaries
TC_summary <- TSSr_clusters[, .(
  TC.start = min(start), 
  TC.end   = max(end),   
  consensus_peak  = round(mean(dominant_tss),  0),  # or weighted average if desired
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
fwrite(TSSr_clusters, paste0(outdir, '/TSSr.TSS.CAGE.all.Cluster.Results.txt'), sep='\t')

#
fwrite(TSSr_clusters_uni, paste0(outdir, '/TSSr.TSS.CAGE.uni.Cluster.Results.txt'), sep='\t')



### Peak shifts
gg <- ggplot(TSSr_clusters, aes(cell_line, peak_shift, colour = cell_line)) +
  geom_boxplot() +
  geom_point() + 
  theme_minimal() +
  facet_nested(cols=vars(Time), rows =vars(strand))
 

ggsave(paste0(outdir, '/TSSr.TSS.CAGE.PeakShifts.jpg'), width = 20, height = 24)




## Plot Enhancers
gg <- ggplot(enhancers.dt, aes(cell_line, D)) + 
  geom_point(aes(fill = cell_line), shape = 21) + 
  coord_flip() +
  theme_minimal() +
  facet_nested(rows=vars(Time), scales='free', independent = F)


ggsave(paste0(outdir, '/TSSr.TSS.CAGE.enhancers.D.jpg'), width = 20, height = 24)


enhancers.stranded <- melt(enhancers.dt, id.vars = c(1:7, 10:13) )
enhancers.stranded[variable == 'tags.m', value := value * -1]

gg <- ggplot(enhancers.stranded, aes(enhancer, value)) + 
  geom_col(aes(fill = variable), position = 'identity') + 
  coord_flip() +
  theme_minimal() +
  facet_nested(rows=vars(cell_line), cols=vars(Time), scales='free', independent = F)


ggsave(paste0(outdir, '/TSSr.TSS.CAGE.enhancers.jpg'), width = 20, height = 24)




















