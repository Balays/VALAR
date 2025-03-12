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

fwrite(tagClusters, 'all.tagClusters.txt', sep='\t')

#exportTSStoBedgraph(myTSSr, data = "processed", format = "bedGraph") 
#exportTSStoBedgraph(myTSSr, data = "processed", format = "BigWig")

#exportClustersTable(myTSSr, data = "tagClusters")
#exportClustersToBed(myTSSr, data = "tagClusters") 

tagClusters <- fread('all.tagClusters.txt')
