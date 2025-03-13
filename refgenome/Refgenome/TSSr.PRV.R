library(BSgenomeForge)
library(TSSr)
library(devtools)
library(BiocManager)
library(GenomeInfoDb)
library(BSgenome)


# Create BSgenome from seed file
#setwd("refgenome/Refgenome/")
#seed_file <- "D:/Armin/MPOX/refgenome/MPXV_MDBIO4.txt"
forgeBSgenomeDataPkg("refgenome/Refgenome/BSgenome.PRV-seed.txt", replace = T)

pkg_dir <- "refgenome/Refgenome/BSgenome.PRV.MdBio.1.0"

#### Change the unvalid 'K' character in the fasta file to 'N'
### create .2bit file in ubuntu and move the file into the correct folder
#file.copy("D:/Armin/PRV_3_cell_CAGE/Refgenome/LT934125.1.2bit",
#          "D:/Armin/PRV_3_cell_CAGE/Refgenome/BSgenome.PRV.MdBio.1.0/inst/extdata/single_sequences.2bit",
#          overwrite = TRUE)

devtools::check(pkg_dir)
devtools::build(pkg_dir)
devtools::install(pkg_dir)






stop()


#load BSgenome package
library(BSgenome.PRV.MdBio.1.0)
setwd("D:/Armin/PRV_3_cell_CAGE/TSSr/")

#set parameters
refSource <- "D:/Armin/PRV_3_cell_CAGE/Refgenome/LT934125.1-2.gff3"
organismName <- "PRV"
directory_path <- "D:/Armin/PRV_3_cell_CAGE/TSSr/"

inputFiles <- list.files(path = "D:/Armin/PRV_3_cell_CAGE/Gmail/", pattern = "\\.bam$", full.names = TRUE)
inputFilesType <- "bam"


sampleLabels <- c('1', '2', '3', '4', '5', '6', '7', '8')

sampleLabelsMerged <- "all" 
mergeIndex <- c(1, 1, 1, 1, 1, 1, 1, 1)              #c(5, 5, 5, 1, 1, 1, 6, 6, 6, 2, 2, 2, 3, 3, 3, 4, 4, 4)



#create TSSr object
myTSSr <- new("TSSr", genomeName = "BSgenome.PRV.MdBio.1.0",
              inputFiles = inputFiles,
              inputFilesType = inputFilesType,
              sampleLabels = sampleLabels,
              sampleLabelsMerged = sampleLabelsMerged,
              mergeIndex = mergeIndex,
              refSource = refSource,
              organismName = organismName
)


getTSS(myTSSr, sequencingQualityThreshold = 10, mappingQualityThreshold = 20)

# export raw data
exportTSStable(myTSSr, data = "raw", merged = "FALSE")

#create correlation plot
plotCorrelation(myTSSr, samples = "all")

plotTssPCA(myTSSr, TSS.threshold=10)


mergeSamples(myTSSr)

myTSSr@librarySizes

filterTSS(myTSSr, method = "poisson")

exportTSStable(myTSSr, data = "processed") 
exportTSStoBedgraph(myTSSr, data = "processed", format = "bedGraph") 
exportTSStoBedgraph(myTSSr, data = "processed", format = "BigWig")


clusterTSS(myTSSr, method = "peakclu",peakDistance=100,extensionDistance=30
           ,localThreshold = 0.02, clusterThreshold = 1
           ,useMultiCore=FALSE, numCores=NULL)


exportClustersTable(myTSSr, data = "tagClusters")
exportClustersToBed(myTSSr, data = "tagClusters") 