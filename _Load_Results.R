## Neccasseary packages
library(data.table)
library(Rsamtools)
library(stringi)

##
#### Before Starting ####

### Path for bamfiles
bamdir  <- "D:/data/PRV_3cell/rebasecall/mapped_v6_virus" 
if (.Platform$OS.type!="windows") {
  bamdir <- paste0('/mnt/', gsub(':', '/', tolower(gsub('/.*', '', bamdir))),
                   stri_replace_first_regex(bamdir, '.*:/', ''))
}


pattern  <- '_stranded_only.bam$' ## this will be subtracted from the filename and supplied as a sample column 
bamfiles <- list.files(bamdir, pattern, recursive = T, full.names = T)

### Project name and output directory
outdir  <- 'PRV-MDBIO-4cell'; try({ dir.create(outdir) })


### ####
##

#### Load configuration file
config <- readRDS(paste0(outdir, '_config.rds'))


##
#### WF Part 0. ####

##
### Packages, functions, settings, metadata and annotation
source('_WF.part0.R')


### Reference Transcripts
## Import and Format 
source('import.ref.TRs.R')
TR.Ref.data <- fread(paste0(outdir, '/TR.Ref.data.tsv'))



#### ####
##


#bam.all <- fread(paste0(outdir, '/bam.all.tsv'), na.strings = '')

mapped.cov       <- fread(paste0(outdir, '/mapped.cov.tsv'), na.strings = '')
merged_cov       <- fread(paste0(outdir, '/merged_cov.tsv'), na.strings = '')

win.cov.sum      <- fread(paste0(outdir, '/win.cov.sum.tsv'), na.strings = '')
win.cov.hpi.sum  <- fread(paste0(outdir, '/win.cov.hpi.sum.tsv'), na.strings = '')

norm_cov_summary <- fread(paste0(outdir, '/norm.cov.summary.tsv'), na.strings = '')

readcounts       <- fread(paste0(outdir, '/readcounts.tsv'), na.strings = '', header = T)


bam.filt <- fread(paste0(outdir, '/bam.filt.tsv'), na.strings = '')



### TransFrags

aln.uni   <- fread(paste0(outdir, '/aln.uni.tsv'),  na.strings = '')
TR.uni    <- fread(paste0(outdir, '/TR.uni.tsv'),  na.strings = '')
EX.uni    <- fread(paste0(outdir, '/EX.uni.tsv'),  na.strings = '')
TR.data   <- fread(paste0(outdir, '/TR.data.tsv'), na.strings = '')
TR.counts <- fread(paste0(outdir, '/TR.counts.tsv'), na.strings = '')


TR.gff <- data.table(as.data.frame(rtracklayer::import.gff2(config$TR.reads.gfffile)))
TR.EX  <- fread(paste0(outdir, '/TR.EX.tsv'))

## Adapter counts per TransFrag
TR.adapt.count <- aln.uni[,.(count=.N), by=.(seqnames, strand, TR_ID, TR_start, TR_end, correct_tss, correct_tes, sample)]

TR.counts.sp <- dcast(TR.adapt.count, TR_ID + correct_tss + correct_tes ~ sample, value.var = 'count', fill=0)
##


### GFF-compare
## all results
all.merged.result_gff.compare  <- fread(paste0(outdir, "/all.merged.result_gff.compare.tsv"), na.strings = '')

## best hits results
best.merged.result_gff.compare <- fread(paste0(outdir, "/best.merged.result_gff.compare.tsv"), na.strings = '')

## 
TR.gff.compare.merged.TR.counts.gt <- fread(paste0(outdir, "/TR.gff.compare.merged.TR.counts.gt.tsv"), na.strings = '')


## adapter counts per transfrags in each sample
adapt.TR.sp <- fread(file.path(outdir, "adapt.TR.sp.tsv"), sep = '\t')


### merge TransFrag count with Reference Transcript annotation for each Transfrag
all.merged.result_gff.compare      <- merge(all.merged.result_gff.compare,      adapt.TR.sp, by.x='transcript_id',  by.y='TR_ID') #, all = T)
best.merged.result_gff.compare     <- merge(best.merged.result_gff.compare,     adapt.TR.sp, by.x='transcript_id',  by.y='TR_ID') #, all = T)
TR.gff.compare.merged.TR.counts.gt <- merge(TR.gff.compare.merged.TR.counts.gt, adapt.TR.sp, by.x='transcript_id',  by.y='TR_ID') #, all = T)

