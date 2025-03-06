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

### Config file creation
config <- list(
  
  ## input bamfiles
  bamdir = bamdir,
  pattern = pattern,
  bamfiles = bamfiles,
  
  ## filter bamfiles for virus ?
  bam_filt_outdir = paste0(bamdir, "_filtered"),
  filter_bams = FALSE,
  
  ## results output - project folder
  outdir = outdir,
  
  ## Neccessary packages and functions
  misc_dir = 'C:/GitHub/Rlyeh/R',
  minitax_dir = 'C:/GitHub/minitax/R',
  
  ## Metadata
  metadata_from_bamfiles = T,
  metadata_file = 'PRV-MDBIO-4cell_metadata.tsv',
  metacols = c('sample', 'rep', 'hpi', 'Time', 'libtype', 'group', 'cell_line', 'sample_name'),
  
  ## Genome and annotation
  genome = 'LT934125.1',
  fasta_ref = 'refgenome/LT934125.1.fasta',
  virus = 'PRV-MDBIO',
  gff_file = 'refgenome/LT934125.1.gff3',
  create.ann.from.gff = T,
  feature.df.file = "refgenome/LT934125.1_feature.df.tsv",
  ## reference transcripts
  viral.ref.gff = "refgenome/LT934125.1_Torma_et_al_corrected.gff3",
  
  ## Other settings
  nproc = 48,
  
  write.all = TRUE,
  save.images = TRUE,
  
  rename_host_contigs = FALSE,
  fix.viral.contigs = TRUE,
  multiCellLines <- F,
  
  make.plots = TRUE,
  
  ## CAGE ?
  include.cage = T,
  
  ## Coverage analysis settings
  window_size = 50,
  window_step = 50,
  
  ## Alignment analysis settings
  is.lortia = TRUE,
  rm.gaps.in.aln = TRUE,
  
  ## Coverage and alignment Analysis settings
  flag  = scanBamFlag(isSupplementaryAlignment=FALSE),
  param = ScanBamParam(what=scanBamWhat(), flag=scanBamFlag(isSupplementaryAlignment=FALSE)),
  
  ## GFF-compare settings
  TR.reads.gfffile = paste0(outdir, '/TR.reads.gff2'),
  
  thresh.eq.prime5 = 10,
  thresh.eq.prime3 = 10,
  thresh.eq.junc   = 2
  
  
)

saveRDS(config, file = paste0(outdir, '_config.rds'))

#### ####
##

###### START WORKFLOW ! #######

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

## Export
export.gff2(as.data.frame(viral.ref), "refgenome/LT934125.1_Torma_et_al_corrected_w_CDS.genes.gff2")
export.gff3(as.data.frame(viral.ref), "refgenome/LT934125.1_Torma_et_al_corrected_w_CDS.genes.gff3")

fwrite(TR.Ref.data, paste0(outdir, '/TR.Ref.data.tsv'), sep = '\t')
TR.Ref.data <- fread(paste0(outdir, '/TR.Ref.data.tsv'))



#### ####
##


##
#### WF Part 1. ####

##
### Import alignments
source('_WF.part1.R')


##
### Import coverages
source('_WF.part2.R')

##
### CAGE Analysis
if (config$include.cage) {
  outdir  <- 'CAGE'; try({ dir.create(outdir) })
  
  config <- source('CAGE_config.R')
  
  is.lortia <- FALSE
  flag  <- scanBamFlag(isSupplementaryAlignment=FALSE)
  param <- ScanBamParam(what=scanBamWhat(), flag=scanBamFlag(isSupplementaryAlignment=FALSE))
  rm.gaps.in.aln <- TRUE
  
  
  ##
  ### Import CAGE alignments
  source('_WF.part1.R')
  
  ##
  ### Import CAGE coverages
  source('_WF.part2.R')
  
  
  #### Load MAIN configuration file again
  config <- readRDS(paste0(outdir, '_config.rds'))
  
  
}


##
### Load back
try({
  
  bam.all <- fread(paste0(outdir, '/bam.all.tsv'), na.strings = '')
  
  mapped.cov       <- fread(paste0(outdir, '/mapped.cov.tsv'), na.strings = '')
  merged_cov       <- fread(paste0(outdir, '/merged_cov.tsv'), na.strings = '')
  
  win.cov.sum      <- fread(paste0(outdir, '/win.cov.sum.tsv'), na.strings = '')
  win.cov.hpi.sum  <- fread(paste0(outdir, '/win.cov.hpi.sum.tsv'), na.strings = '')
  
  norm_cov_summary <- fread(paste0(outdir, '/norm.cov.summary.tsv'), na.strings = '')
  
  readcounts       <- fread(paste0(outdir, '/readcounts.tsv'), na.strings = '', header = T)
  
})

### remove
rm(list = c('bam.all.list', 'bam.cov.list', 'mapped.cov', 'merged_cov', 'norm.cov'))

#### ####
##

### Combine w CAGE
if(include.cage) {
  readcounts       <- rbind(readcounts, readcounts.CAGE)
  #bam.filt         <- rbind(bam.filt,   bam.filt.cage)
}


### SAVE / LOAD IMAGE
gc()
if(save.images) {
  save.image(paste0(outdir, '.P1.RData'))
  load(paste0(outdir, '.P1.RData'))
}
###

## add correct tags
if(config$is.lortia) {
  bam.all[,correct_tes := fifelse(grepl('correct', tag.l3) | grepl('correct', tag.r3), T, F)]
  bam.all[,correct_tss := fifelse(grepl('correct', tag.l5) | grepl('correct', tag.r5), T, F)]
}
fwrite(bam.all, paste0(outdir, '/bam.all.tsv'), sep='\t')


bam.all <- fread(paste0(outdir, '/bam.all.tsv'))
## filter out reads with supplementary alignments, or
## these were filtered out previously?
# bam.filt[qname]
bam.filt <- bam.all [!is.na(seqnames)]
bam.filt <- bam.filt[seqnames == genome]


#rm('bam.all')

fwrite(bam.filt, paste0(outdir, '/bam.filt.tsv'), sep='\t')

bam.filt <- fread(paste0(outdir, '/bam.filt.tsv'), na.strings = '')



##
#### WF Part 2. Adapter Checking ####

## check LoRTIA tags
source('LoRTIA_all.bam_barplot_v2.R')


## check the softclips
source('softclip_check.R')
head(bam.sc.correct)
bam.sc.correct[,.N,softclip_length][order(softclip_length)]

fwrite(bam.sc.correct, file.path(outdir, 'bam.sc.correct.tsv'), sep = '\t')

## adapter distance from the start of mapping
bam.sc.correct <- fread(file.path(outdir, 'bam.sc.correct.tsv'), sep = '\t')
source('adapter_check.R')
head(bam.sc.correct.adapt)
bam.sc.correct.adapt[,.N,by=.(adapter_found)]
bam.sc.correct.adapt[,.N,by=.(adapter_found, adapter_distance)][order(adapter_distance)]

fwrite(bam.sc.correct.adapt, file.path(outdir, 'bam.sc.correct.adapt.tsv'), sep = '\t')

## polyA size
bam.sc.correct.adapt <- fread(file.path(outdir, 'bam.sc.correct.adapt.tsv'), sep = '\t')
source('5polyA_check.R')
head(bam.sc.correct.adapt.polyA)
bam.sc.correct.adapt.polyA[,.N,by=.(polyA_distance)][order(polyA_distance)]
bam.sc.correct.adapt.polyA[,.N,by=.(polyA_length)]  [order(polyA_length)]

bam.sc.correct.adapt.polyA[,polyA_to_softclip_percent := round(polyA_to_softclip_ratio*100, 0)]

fwrite(bam.sc.correct.adapt.polyA, file.path(outdir, 'bam.sc.correct.adapt.polyA.tsv'), sep = '\t')


## make plots
bam.sc.correct.adapt       <- fread(file.path(outdir, 'bam.sc.correct.adapt.tsv'), sep = '\t')
bam.sc.correct.adapt.polyA <- fread(file.path(outdir, 'bam.sc.correct.adapt.polyA.tsv'), sep = '\t')

source('adapter_plots.R')


## filter 5-ends
source('correct_5p.r')

## merg results into bam.all
bam.all.m <- merge(bam.all, bam.sc.correct.adapt.polyA[,.(seqnames, strand, qname, correct_5p)], 
                 by=c('seqnames', 'strand', 'qname'), all=T)

if(nrow(bam.all) == nrow(bam.all.m)) {
  bam.all <- bam.all.m; rm(bam.all.m)  
} else { stop() }

fwrite(bam.all, file.path(outdir, 'bam.all.tsv'), sep = '\t')
saveRDS(bam.all, file.path(outdir, 'bam.all.rds'))


##
## filter bam.all to correct 5-prime OR 3-prime reads

## can there be reads whose 3-prime polyA is actually a 5-prime polyA? 
## or because of the 15 stretch of A in the LoRTIA command these will be filtered out?


bam.filt <- bam.all[correct_5p == T | 
                   (strand == '-' & grepl('correct', tag.l3)) |
                   (strand == '+' & grepl('correct', tag.r3))
                   ,]


## end statistics on a read level
reads.filt <- unique(bam.filt[,.(qname, strand, aln_ID, correct_5p, correct_tes)])
reads.filt[,.N, by=.(correct_5p, correct_tes)]

## NA 5-prime --> FALSE
bam.filt[,correct_5p := fifelse(is.na(correct_5p), F, correct_5p)]

## end statistics on a read level
reads.filt <- unique(bam.filt[,.(qname, strand, aln_ID, correct_5p, correct_tes)])
reads.filt[,.N, by=.(correct_5p, correct_tes)]


fwrite(bam.filt, file.path(outdir, 'bam.filt.tsv'), sep = '\t')

#### ####
##

bam.filt <- fread(file.path(outdir, 'bam.filt.tsv'), na.strings = '')


##
#### WF Part 3. Read clustering into TransFrags ####

## make Transcripts
source('makeTRs.R')
TR.uni[, TR_start := min(start), by =.(TR_ID)][, TR_end := max(end), by= .(TR_ID)]
TR.uni[,prime5.TR := fifelse(strand == '+', TR_start, TR_end)][,prime3.TR := fifelse(strand == '-', TR_start, TR_end)]
EX.uni <- TR.uni

TR.uni <- unique(TR.uni[,.(seqnames, TR_start, TR_end, strand, TR_ID, prime5.TR, prime3.TR)])

TR.data <- EX.uni
#
TR.data <- dcast(TR.data, seqnames + strand + TR_ID + TR_start + TR_end + exon_combination ~ exon_rank, value.var = 'exon_ID')
TR.data <- merge(TR.data, TR.counts, by='TR_ID', all=T)


## TR counts with adapter info per sample
aln.uni      <- unique(bam.TR[,.(seqnames, strand, aln_ID, TR_ID, TR_start, TR_end, correct_tss, correct_tes, sample)])

## write out
fwrite(aln.uni,   paste0(outdir, '/aln.uni.tsv'),  sep='\t')
fwrite(EX.uni,    paste0(outdir, '/EX.uni.tsv'),  sep='\t')
fwrite(TR.uni,    paste0(outdir, '/TR.uni.tsv'),  sep='\t')
fwrite(TR.data,   paste0(outdir, '/TR.data.tsv'), sep='\t')
fwrite(TR.counts, paste0(outdir, '/TR.counts.tsv'), sep='\t')

## load back
aln.uni   <- fread(paste0(outdir, '/aln.uni.tsv'),  na.strings = '')
TR.uni    <- fread(paste0(outdir, '/TR.uni.tsv'),  na.strings = '')
EX.uni    <- fread(paste0(outdir, '/EX.uni.tsv'),  na.strings = '')
TR.data   <- fread(paste0(outdir, '/TR.data.tsv'), na.strings = '')
TR.counts <- fread(paste0(outdir, '/TR.counts.tsv'), na.strings = '')


### Export Transfrags
source('export.gffs.R')

TR.gff <- data.table(as.data.frame(rtracklayer::import.gff2(config$TR.reads.gfffile)))
TR.EX  <- fread(paste0(outdir, '/TR.EX.tsv'))


## CAGE
#CAGE.TR.data <- TR.data
#


### Remove
rm(list = c('dt', 'reads'))

#### ####
##



### SAVE / LOAD IMAGE
gc()
if(save.images) {
  save.image(paste0(outdir, '.P2.RData'))
  load(paste0(outdir, '.P2.RData'))
}
###




##
#### WF Part 3. Validate 3-primes based on reference transcripts TES
valid.TES.win <- 10
validate_TES  <- T

if (validate_TES) {

  TR.ref[,transcript_prime3 := fifelse(strand == '+', transcript_end, transcript_start)]
  TR.ref[,transcript_prime5 := fifelse(strand == '-', transcript_end, transcript_start)]

  valid.prime3 <- unique(TR.ref[,.(seqnames, transcript_start, transcript_end, transcript_prime3, transcript_prime5, strand)])
  valid.prime3 <- valid.prime3[,.(seqnames, strand, start = transcript_prime3 - valid.TES.win, end = transcript_prime3 + valid.TES.win)]
  valid.prime3 <- unique(valid.prime3)
  valid.prime3[,valid_tes := T]

  aln.uni[,TR_prime3 := fifelse(strand == '+', TR_end, TR_start )]
  aln.uni[,TR_prime5 := fifelse(strand == '-', TR_end, TR_start )]
  aln.uni[,start := TR_prime3]
  aln.uni[,end   := TR_prime3]

  prime3.TR.ov <- foverlaps2(aln.uni, valid.prime3, by.x=c('seqnames', 'strand', 'start', 'end'), by.y=c('seqnames', 'strand', 'start', 'end'), minoverlap = 1)
  prime3.TR.ov <- prime3.TR.ov[,.(seqnames, strand, TR_start, TR_end, TR_prime3, TR_prime5, correct_tss, correct_tes, aln_ID, TR_ID, start, end, sample)]
  prime3.TR.ov <- merge(prime3.TR.ov, valid.prime3, by=c('seqnames', 'strand', 'start', 'end'), all.x=T)

  prime3.TR.ov[, start  := NULL]
  prime3.TR.ov[, end    := NULL]
  aln.uni[, start       := NULL]
  aln.uni[, end         := NULL]


  prime3.TR.ov <- merge(prime3.TR.ov, aln.uni, by=colnames(aln.uni), all=T)
  prime3.TR.ov[,valid_tes := fifelse(is.na(valid_tes), F, T)]

  prime3.valid.corr.freq <- prime3.TR.ov[,.N, by=.(sample, correct_tes, valid_tes)]

  ggvf <- ggplot(prime3.valid.corr.freq) +
    geom_col(aes(x=sample, y=N, fill=correct_tes), color='black') +
    coord_flip() +
    facet_wrap(~valid_tes, nrow=1) +
    theme_bw() +
    ggtitle('3-prime end result, according to ref mRNA TES')

  ggsave(file.path(outdir, 'Ref_mRNA_TES_validation.jpg'), ggvf, height = 12, width = 9)

  prime3.TR.ov[, correct_tes := fifelse(valid_tes == T | correct_tes == T, T, F)]
  prime3.TR.ov[, valid_tes   := NULL]

  keyby <- colnames(prime3.TR.ov)
  prime3.TR.ov <- unique(prime3.TR.ov[,]) ##.(), by=.()]

  if(nrow(aln.uni) == nrow(prime3.TR.ov)) {
    aln.uni      <- prime3.TR.ov
  } else { stop() }

}

TR.adapt.count <- aln.uni[,.(count=.N), by=.(seqnames, strand, TR_ID, TR_start, TR_end, correct_tss, correct_tes, sample)]
fwrite(TR.adapt.count, paste0(outdir, '/TR.adapt.count.tsv'), sep = '\t')

TR.adapt.count <- fread(paste0(outdir, '/TR.adapt.count.tsv'), na.strings = '')

TR.counts.sp <- dcast(TR.adapt.count, TR_ID + correct_tss + correct_tes ~ sample, value.var = 'count', fill=0)
##

#orf5 <- bam.TR[TR_end >= 46795 & TR_end <= 46995 & grepl('12h', sample) & strand == '-',]


### whats thiS?
#TR.data2 <- prime3.TR.ov[,.(count=.N), by=.(TR_ID, sample, correct_tss, correct_tes)]

#### ####
##


### SAVE / LOAD IMAGE
if(save.images) {
  save.image(paste0(outdir, '.P2.RData'))
  load(paste0(outdir, '.P2.RData'))
}



##
#### WF part 4. Ref Transcript count analysis ####

### Run GFF-compare on each ref TR separately,
## import results and calculate distances
source('run_GFF.COMPARE.R')

### Analyse GFF-compare results
## read in GFF-compare results
all.merged.result_gff.compare  <- fread(paste0(outdir, "/all.merged.result_gff.compare.tsv"), na.strings = '')

# Update the results as transcript IDs in the reference annotation was changed
#source('update_TR.ref.IDs.R')

## find the closest ref-TR for each query
## categorise non-equal matches
thresh.eq.prime5 <- config$thresh.eq.prime5
thresh.eq.prime3 <- config$thresh.eq.prime3
thresh.eq.junc   <- config$thresh.eq.junc

source('analyse_GFF.COMPARE.R')
## import results
best.merged.result_gff.compare <- fread(paste0(outdir, "/best.merged.result_gff.compare.tsv"), na.strings = '')
#best.merged.result_gff.compare[,transcript_id := as.integer(transcript_id) ]

### Summarise results
source('summarise_GFF.COMPARE.R')

TR.gff.compare.merged.TR.counts.gt <- fread(paste0(outdir, "/TR.gff.compare.merged.TR.counts.gt.tsv"), na.strings = '')

### include LoRTIA adapter info to Transcripts
adapt.TR <- aln.uni[,.N,by=.(TR_ID, correct_tss, correct_tes)]
adapt.TR[,correct_tss := paste0('correct_tss::', as.character(correct_tss))]
adapt.TR[,correct_tes := paste0('correct_tes::', as.character(correct_tes))]
adapt.TR[,adapter := paste0(correct_tes, ';', correct_tss)]
adapt.TR.sp <- dcast(adapt.TR, TR_ID ~ adapter, value.var = 'N', fill=0)


fwrite(adapt.TR.sp, file.path(outdir, "adapt.TR.sp.tsv"), sep = '\t')


### merge TransFrag count with Reference Transcript annotation for each Transfrag
all.merged.result_gff.compare      <- merge(all.merged.result_gff.compare,      adapt.TR.sp, by.x='transcript_id',  by.y='TR_ID') #, all = T)
best.merged.result_gff.compare     <- merge(best.merged.result_gff.compare,     adapt.TR.sp, by.x='transcript_id',  by.y='TR_ID') #, all = T)
TR.gff.compare.merged.TR.counts.gt <- merge(TR.gff.compare.merged.TR.counts.gt, adapt.TR.sp, by.x='transcript_id',  by.y='TR_ID') #, all = T)

# Note: those transfrags, that were not assigned to any Ref-TR by gff-compare, are not listed here !
# If we want to include those also, we need to change the code!

#### ####
##

### remove
rm(list=c('all.merged.result_gff.compare', 'TR.gff.compare.merged.TR', 'TR.gff.compare.merged.TR.counts'))
gc()

### SAVE / LOAD IMAGE
if(save.images) {
  save.image(paste0(outdir, '.P3.RData'))
  load(paste0(outdir, '.P3.RData'))
}




##
#### WF part 4. Count genes and Transcripts ####

#source('_WF.part0.R')

res.dir <- outdir; try({ dir.create(res.dir) })

## filter some genes ?
genes_to_filter <- 'CTO-L'

## treat spliced transcripts differently?
check.spliced.TRs <- F

### 1.) Count complete CDSs
TR.EX.df <- copy(TR.EX)
colnames(TR.EX.df)[3] <- 'TR_ID'
TR.EX.df[,EX_ID := paste0('EX_', .GRP), by =.(seqnames, start, end, strand)]

TR.CDS.OV.complete <- feature.OV.from.polyC.TR(feature.df, as.data.frame(TR.EX.df), EX.only = F, feature.colname = 'gene', antisense = F, type = 'within')
TR.CDS.OV.complete.sp <- TR.CDS.OV.complete[[4]]
TR.CDS.OV.complete    <- TR.CDS.OV.complete[[3]]
TR.CDS.OV.complete$TR_ID <- as.integer(TR.CDS.OV.complete$TR_ID)

## merge with transcript data (without counts)
TR.data.CDS.OV.complete    <- merge(
  unique(TR.data[,.SD, .SDcols = !c('sample', 'count')]), by.x= c('seqnames', 'strand', 'TR_ID'),
  TR.CDS.OV.complete,  by.y = c('seqnames', 'strand.TR', 'TR_ID'),
  all=T)


## merge with transcript data (with counts)
TR.data.CDS.OV.complete.counts  <- merge(
  unique(TR.data[,]), by.x= c('seqnames', 'strand', 'TR_ID'),
  TR.CDS.OV.complete,  by.y = c('seqnames', 'strand.TR', 'TR_ID'),
  all=T)

## put in adapter data ?


## summarise counts per gene and sample (first feature per transcript)

gene.sample_count <- TR.data.CDS.OV.complete.counts[feature_n == 'gene_1', .(gene_count = sum(count)), by=.(seqnames, strand, feature, sample)]
gene.sample_count <- merge(gene.sample_count, metafilt, by='sample')

## kinetic classes ? @Ármin

gene.sample_count[,Kinetic_class := 'unknown'][,gene_cluster := feature][,gene:=feature]

gene.sample_count.sp <- as.data.table(spread(gene.sample_count[,.(seqnames, strand, gene_cluster, gene, Kinetic_class, sample, gene_count)], 'sample', 'gene_count', fill=0))

fwrite(gene.sample_count.sp, paste0(res.dir, '/gene.sample.count.sp.tsv'), sep = '\t')



### 2.) Count genes based on TSS and/or TES

## Import canonic TSS and TES data
#gene.clusters.all <- fread('PRV.genes.TSS.TES.txt')
gene.clusters.all[,gene.region.start := as.integer(ifelse(strand == '+', TSS.canonic, TES.canonic))]
gene.clusters.all[,gene.region.end   := as.integer(ifelse(strand == '-', TSS.canonic, TES.canonic))]
### Dereplicate spliced genes !!!
gene.clusters.all <- unique(gene.clusters.all[,.(seqnames, TSS.canonic, TES.canonic, cluster_TES, strand, gene, gene_cluster, Kinetic_class, TSS.win.start, TSS.win.end, TES.win.start, TES.win.end)])


source('count.genes.v2.R')
## OK


######## CAGE
source('CAGEs.R')



##### ITT JÁROK, EDDIG OK


### 3.) Count transcripts and combine with gene (and cluster) counts
source('count.transcripts.R')




gc()

### SAVE / LOAD IMAGE
if(save.images) {
  save.image(paste0(outdir, '.P4.RData'))
  load(paste0(outdir, '.P4.RData'))
}

##### FINISH #######




## Count coverages of GENE REGIONS (FROM TSS TO TES)
source('count.gene.cov.R')
cov.gene.ov.counts <- fread(paste0(outdir, '/cov.gene.ov.counts.tsv'), sep = '\t')


## Count coverages of CDS
source('count.CDS.cov.R')
cov.CDS.ov.counts <- fread(paste0(outdir, '/cov.CDS.ov.counts.tsv'), sep = '\t')


## Count overlaps of reference transcripts
source('count.overlaps.R')


## Finding ANY 10 Overlaps between exons and GENE REGIONS (FROM TSS TO TES)
source('count.generegions.any10.R')
TR.any10.generegion.counts <- fread(paste0(outdir, '/TR.any10.generegion.counts.tsv'), sep = '\t')


## Finding ANY 10-nt Overlaps between exons and ORFs
source('count.ORFs.any10.R')
TR.any10.ORF.counts <- fread(paste0(outdir, '/TR.any10.ORF.counts.tsv'), sep = '\t')


## Finding WITHIN Overlaps between exons and ORFs
source('count.ORFs.within.R')
TR.within.ORF.counts <- fread(paste0(outdir, '/TR.within.ORF.counts.tsv'), sep = '\t')

#### Analysis ready ! ####
##

gc()

### SAVE / LOAD IMAGE
if(save.images) {
  save.image(paste0(outdir, '.P4.RData'))
  load(paste0(outdir, '.P4.RData'))
}

##
#### Plotting ####

### GFF-compare results
source('plot_GFF.COMPARE.R')

### Reference transcript's overlaps
source('plot.overlaps.R')

### Venn diagram of shared transcripts in the cell lines
source('plot.TR.Venn.R')

#### ####
##


##
#### Coverage and read end plots ####

### CAGE
CAGE.TR.data <- fread('CAGE/TR.data.tsv')
CAGE.TR.data <- merge(CAGE.TR.data, meta.cage, by='sample')
setnames(CAGE.TR.data, new=c('start', 'end'), old=c('TR_start', 'TR_end'),skip_absent=TRUE)

CAGE.TR.data[,prime3 := ifelse(strand == '+', end,   start)]
CAGE.TR.data[,prime5 := ifelse(strand == '+', start, end)]

cols_to_group <- c(metacols, 'seqnames', 'strand', 'prime5')
prime5.counts <- CAGE.TR.data[count>0, .(count=sum(count)), by=cols_to_group][order(seqnames, strand, prime5)]
setnames(prime5.counts, 'prime5', 'pos')
prime5.counts[,endtype := 'prime5']

cols_to_group <- c(metacols, 'seqnames', 'strand', 'prime3')
prime3.counts <- CAGE.TR.data[count>0, .(count=sum(count)), by=cols_to_group][order(seqnames, strand, prime3)]
setnames(prime3.counts, 'prime3', 'pos')
prime3.counts[,endtype := 'prime3']

prime.counts  <- rbind(prime5.counts, prime3.counts)
prime.counts.CAGE <- prime.counts
#prime.counts.CAGE$group <- paste0('CAGE_', prime.counts.CAGE$cell_line, '_', prime.counts.CAGE$hpi)

fwrite(prime.counts.CAGE, paste0(outdir, '/prime.counts.CAGE.tsv'))

### Read end counts from read-transcripts
TR.gff.compare.merged.TR.counts.gt <- fread(paste0(outdir, "/TR.gff.compare.merged.TR.counts.gt.tsv"))

TR.gff.compare.merged.TR.counts.gt[,prime3 := fifelse(strand == '+', end,   start)]
TR.gff.compare.merged.TR.counts.gt[,prime5 := fifelse(strand == '+', start, end)]

cols_to_group <- c(metacols, 'seqnames', 'strand', 'prime5')
prime5.counts <- TR.gff.compare.merged.TR.counts.gt[count>0,.(count=sum(count)),by=cols_to_group][order(seqnames, strand, prime5)]
setnames(prime5.counts, 'prime5', 'pos')
prime5.counts[,endtype := 'prime5']

cols_to_group <- c(metacols, 'seqnames', 'strand', 'prime3')
prime3.counts <- TR.gff.compare.merged.TR.counts.gt[count>0,.(count=sum(count)),by=cols_to_group][order(seqnames, strand, prime3)]
setnames(prime3.counts, 'prime3', 'pos')
prime3.counts[,endtype := 'prime3']

#### OR ! From LoRTIA outpu, to include adapter info !!!

### Read end counts from bamfiles, separated by adapters
bam.counts <- bam.TR
bam.counts[,prime3 := fifelse(strand == '+', end,   start)]
bam.counts[,prime5 := fifelse(strand == '+', start, end)]
bam.counts[,correct_tss := fifelse(grepl('correct', tag.r5) | grepl('correct', tag.l5), T,   F)]
bam.counts[,correct_tes := fifelse(grepl('correct', tag.r3) | grepl('correct', tag.l3), T,   F)]

cols_to_group <- c('sample', 'seqnames', 'TR_ID', 'aln_ID', 'strand', 'TR_start', 'TR_end', 'correct_tss', 'correct_tes')

aln.uni <- unique(bam.counts[,..cols_to_group])



#### START FROM HERE

aln.uni[,TR_prime3 := fifelse(strand == '+', TR_end,   TR_start)]
aln.uni[,TR_prime5 := fifelse(strand == '+', TR_start, TR_end)]

prime5.bam.counts <- aln.uni[,.(count=.N), by=.(seqnames, strand, TR_prime5, correct_tss, sample)]; setnames(prime5.bam.counts, 'TR_prime5', 'pos')
prime3.bam.counts <- aln.uni[,.(count=.N), by=.(seqnames, strand, TR_prime3, correct_tes, sample)]; setnames(prime3.bam.counts, 'TR_prime3', 'pos')

prime5.bam.counts <- merge(prime5.bam.counts, metafilt[,metacols], by='sample')
prime5.bam.counts[, endtype := 'prime5']

prime3.bam.counts <- merge(prime3.bam.counts, metafilt[,metacols], by='sample')
prime3.bam.counts[, endtype := 'prime3']

##
prime3.counts <- prime3.bam.counts
prime5.counts <- prime5.bam.counts

fwrite(prime3.counts, paste0(outdir, '/prime3.counts.tsv'))
fwrite(prime5.counts, paste0(outdir, '/prime5.counts.tsv'))


### Combine
prime.counts  <- rbind(prime5.counts, prime3.counts, fill=T)
fwrite(prime.counts, paste0(outdir, '/prime.counts.tsv'))


### Combine w CAGE ??
prime5.counts <-  rbind(prime.counts[endtype == 'prime5'], prime.counts.CAGE[endtype == 'prime5'])
prime3.counts <-  rbind(prime.counts[endtype == 'prime3'] ) #, prime.counts.CAGE[endtype == 'prime3'])
prime.counts  <-  rbind(prime.counts, prime.counts.CAGE)


#### Make 5- and 3-prime end plots
prime5.counts <- fread(paste0(outdir, '/prime5.counts.tsv'), na.strings = '')
prime3.counts <- fread(paste0(outdir, '/prime3.counts.tsv'), na.strings = '')

## Filter out false adaptered reads?
#prime5.counts <- prime5.counts[correct_tss == T,]
#prime3.counts <- prime3.counts[correct_tes == T,]

## OR !! use ref_mRNAs to accept 3-primes?
TR.ref[,transcript_prime3 := fifelse(strand == '+', transcript_end, transcript_start)]
TR.ref[,transcript_prime5 := fifelse(strand == '-', transcript_end, transcript_start)]

valid.prime3 <- unique(TR.ref[,.(seqnames, transcript_start, transcript_end, transcript_prime3, transcript_prime5, strand)])
valid.prime3 <- valid.prime3[,.(seqnames, strand, start = transcript_prime3 - 10, end = transcript_prime3 + 10)]

prime3.counts[,start := pos]
prime3.counts[,end   := pos]
prime3.TR.ov <- foverlaps2(prime3.counts, valid.prime3, by.x=c('seqnames', 'strand', 'start', 'end'), by.y=c('seqnames', 'strand', 'start', 'end'), minoverlap = 1)
prime3.TR.ov <- unique(prime3.TR.ov[,.(seqnames, strand, pos, correct_tes, endtype, sample, hpi, Time, cell_line, group, count)])

valid.prime3 <- unique(prime3.TR.ov[,.(seqnames, strand, pos, valid_tes=T)])

prime3.counts <- merge(prime3.counts, valid.prime3, by=c('seqnames', 'strand', 'pos'), all.x=T)
prime3.counts[,valid_tes := fifelse(is.na(valid_tes), F, T)]

prime3.valid.corr.freq <- prime3.counts[,.(sum_count = sum(count)), by=.(correct_tes, valid_tes)]

prime3.counts[,correct_tes := fifelse(valid_tes == T | correct_tes == T, T, F)]
prime3.counts[,valid_tes := NULL]


## Overwrite?
fwrite(prime3.counts, paste0(outdir, '/prime3.counts.tsv'), sep = '\t')
fwrite(prime5.counts, paste0(outdir, '/prime5.counts.tsv'), sep = '\t')

## Ovewrite TR count table?
valid.prime3.TR <- unique(merge(valid.prime3, TR.uni, by.x=c('seqnames', 'strand', 'pos'), by.y=c('seqnames', 'strand', 'prime3.TR'), all.x=T))
valid.prime3.TR <- valid.prime3.TR[,.(seqnames, strand, TR_ID, valid_tes)]
TR.adapt.count  <- merge(TR.adapt.count, valid.prime3.TR, by=c('seqnames', 'strand', 'TR_ID'), all.x=T)
TR.adapt.count[,valid_tes   := fifelse(is.na(valid_tes), F, T)]
TR.adapt.count[,correct_tes := fifelse(valid_tes == T | correct_tes == T, T, F)]
TR.adapt.count[,valid_tes   := NULL]

TR.adapt.count <- TR.adapt.count[,.(count=sum(count)), by=.(seqnames, strand, TR_ID, correct_tss, correct_tes, sample)]

TR.adapt.count <- merge(TR.uni, TR.adapt.count, by=c('seqnames', 'strand','TR_ID'))
TR.adapt.count[,prime5.TR := NULL]
TR.adapt.count[,prime3.TR := NULL]

fwrite(TR.adapt.count, paste0(outdir, '/TR.adapt.count.tsv'), sep = '\t')

## include adapter counts
TR.counts.sp <- dcast(TR.adapt.count, TR_ID + TR_start + TR_end + correct_tss + correct_tes ~ sample, value.var = 'count', fill = 0)


## use CAGE to accept 5-prime sites?
#     Answer: No.


##### NEW COUNT TABLE

## prime5 and prime3 counts from TR-table, considering adapters

prime.counts  <- merge(TR.adapt.count, metafilt, by='sample')
prime.counts[,start  := TR_start]
prime.counts[,end    := TR_end  ]
prime.counts[,prime5 := fifelse(strand == '+', start, end)]
prime.counts[,prime3 := fifelse(strand == '+', end,   start)  ]

prime3.counts <- prime.counts[,.(count=sum(count)), by=.(seqnames,	strand,	correct_tes, prime3, sample, hpi, Time,	cell_line,	group)]
prime3.counts[, endtype := 'prime3'][, pos := prime3]
prime3.counts[,start  := pos]
prime3.counts[,end    := pos]


prime5.counts <- prime.counts[,.(count=sum(count)), by=.(seqnames,	strand,	correct_tss, prime5, sample, hpi, Time,	cell_line,	group)]
prime5.counts[, endtype := 'prime5'][, pos := prime5]
prime5.counts[,start  := pos]
prime5.counts[,end    := pos]

## Mean coverage from stranded only bamfiles directly
cov.counts <- merged_cov[,.(count=mean(count)), by=.(seqnames,	strand,	pos, sample, hpi, Time,	cell_line,	group)]
cov.counts[,start  := pos]
cov.counts[,end    := pos]


## write
fwrite(prime3.counts, paste0(outdir, '/prime3.counts.tsv'), sep = '\t')
fwrite(prime5.counts, paste0(outdir, '/prime5.counts.tsv'), sep = '\t')
fwrite(cov.counts,    paste0(outdir, '/cov.counts.tsv'),    sep = '\t')



#### PLOTTING
source('plot_primes.R')





















##### START PLOTTING FROM HERE #######
## read back in
prime5.counts <- fread(paste0(outdir, '/prime5.counts.tsv'), na.strings = '')
prime3.counts <- fread(paste0(outdir, '/prime3.counts.tsv'), na.strings = '')
cov.counts    <- fread(paste0(outdir, '/cov.counts.tsv'),    na.strings = '')


## filter now
correct.only <- T
if(correct.only) {
  prime5.counts <- prime5.counts[correct_tss == T,]
  prime3.counts <- prime3.counts[correct_tes == T,]
}


### normalize for viral read counts?
## based on correct adaptered-only
norm <- F
if(norm) {
  colsby <- colnames(prime5.counts)[!colnames(prime5.counts) %in% c('pos', 'start', 'end', 'prime5', 'count')]
  prime5.counts[, sum_count := sum(count), by=colsby]
  prime5.counts[, count     := round((count/ sum_count) * 100, 3)]
  prime5.counts[, .(sum_count = sum(count)), by=colsby]
  prime5.counts[, sum_count := NULL]

  colsby <- colnames(prime3.counts)[!colnames(prime3.counts) %in% c('pos', 'start', 'end', 'prime3', 'count')]
  prime3.counts[, sum_count := sum(count), by=colsby]
  prime3.counts[, count     := round((count/ sum_count) * 100, 3)]
  prime3.counts[, .(sum_count = sum(count)), by=colsby]
  prime3.counts[, sum_count := NULL]

}

all(round(prime5.counts[,.(sum_count = (sum(count))), by=.(sample)][,sum_count], 0) == 1)

##calculate the mean? otherwise, the plot will SUM the counts based on the "group" column
calc.mean <- T
##
if(calc.mean) {
  prime3.counts <- prime3.counts[,.(count = round(mean(count), 2)), by=.(seqnames, strand, pos, correct_tes, hpi, Time, cell_line, group, endtype)]
  prime5.counts <- prime5.counts[,.(count = round(mean(count), 2)), by=.(seqnames, strand, pos, correct_tss, hpi, Time, cell_line, group, endtype)]

  #prime5.counts <- unique(prime5.counts[,.(seqnames, strand, pos, correct_tss, count, hpi, Time, cell_line, group, endtype)])
  #prime3.counts <- unique(prime3.counts[,.(seqnames, strand, pos, correct_tes, count, hpi, Time, cell_line, group, endtype)])

}


##### plotting settings
##
plot.all.together <- T

## grouping
combine.groups <- 'group' ## 'sample' ## OR:

## adapter settings
adapter.setting <- 'v4'


#### Normalized
file_name_suffix <- 'Mean_correct_norm'
##
ylim       <- NULL #c(0, 50) # c(-1000, 1000)

## prime5
source('plot_prime5_area_settings.R')
source('plot_prime5.R')

Fig1B <- Fig1

## prime3
source('plot_prime3_area_settings.R')
source('plot_prime3.R')

Fig2B <- Fig2



#### ylim500
file_name_suffix <- 'Mean_correct_ylim500'
##
ylim       <- c(0, 500) # c(-1000, 1000)

## prime5
source('plot_prime5_area_settings.R')
source('plot_prime5.R')

Fig1A <- Fig1

## prime3
source('plot_prime3_area_settings.R')
source('plot_prime3.R')

Fig2A <- Fig2



####### Supp Figs

#### ylim50
file_name_suffix <- 'Mean_correct_ylim50'
##
ylim       <- c(0, 50) # c(-1000, 1000)

## prime5
source('plot_prime5_area_settings.R')
source('plot_prime5.R')

SFig1A <- Fig1

## prime3
source('plot_prime3_area_settings.R')
source('plot_prime3.R')

SFig2A <- Fig2


#### ylim5000
file_name_suffix <- 'Mean_correct_ylim5000'
##
ylim       <- c(0, 5000) # c(-1000, 1000)

## prime5
source('plot_prime5_area_settings.R')
source('plot_prime5.R')

SFig1B <- Fig1

## prime3
source('plot_prime3_area_settings.R')
source('plot_prime3.R')

SFig2B <- Fig2

#
#### Coverage plots

##
file_name_suffix <- 'Mean_correct_ylim5000'
ylim       <- c(-5000, 5000) # c(-1000, 1000)

source('plot_coverage_area_settings.R')
source('plot_coverage.R')

SFig3A <- Fig3

##
file_name_suffix <- 'Mean_correct_ylim500'
ylim       <- c(-500, 500) # c(-1000, 1000)

source('plot_coverage_area_settings.R')
source('plot_coverage.R')

SFig3B <- Fig3

##
file_name_suffix <- 'Mean_correct_ylim50'
ylim       <- c(-50, 50) # c(-1000, 1000)

source('plot_coverage_area_settings.R')
source('plot_coverage.R')

SFig3C <- Fig3

##
file_name_suffix <- 'Mean_correct_norm'
ylim       <- NULL # c(-1000, 1000)

source('plot_coverage_area_settings.R')
source('plot_coverage.R')

SFig3D <- Fig3

#

## done
######


#### Small Figures for the article
fig.width  <- 30
fig.height <- 20*2

Fig1   <- cowplot::plot_grid(Fig1A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             Fig1B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Figure 1.small.jpg', Fig1, height = fig.height, width = fig.width, limitsize = F)

Fig2  <- cowplot::plot_grid(Fig2A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            Fig2B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Figure 2.small.jpg', Fig2, height = fig.height, width = fig.width, limitsize = F)


SFig1 <- cowplot::plot_grid(SFig1A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            SFig1B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Supp Fig 1.small.jpg', SFig1, height = fig.height, width = fig.width, limitsize = F)


SFig2 <- cowplot::plot_grid(SFig2A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            SFig2B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Supp Fig 2.small.jpg', SFig2, height = fig.height, width = fig.width, limitsize = F)



#### Normal Figures for the article
fig.width  <- 55
fig.height <- 20*2

Fig1   <- cowplot::plot_grid(Fig1A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             Fig1B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Figure 1.jpg',  Fig1, height = fig.height, width = fig.width, limitsize = F)

Fig2  <- cowplot::plot_grid(Fig2A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            Fig2B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Figure 2.jpg',  Fig2, height = fig.height, width = fig.width, limitsize = F)


SFig1 <- cowplot::plot_grid(SFig1A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            SFig1B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Supp Fig 1.jpg', SFig1, height = fig.height, width = fig.width, limitsize = F)


SFig2 <- cowplot::plot_grid(SFig2A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            SFig2B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Supp Fig 2.jpg', SFig2, height = fig.height, width = fig.width, limitsize = F)


##### Sup Fig 3
fig.width  <- 30
fig.height <- 20*3

SFig3  <- cowplot::plot_grid(SFig3A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             SFig3B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             SFig3C + draw_label("c", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Supp Fig 3.small.jpg', SFig3, height = fig.height, width = fig.width, limitsize = F)

fig.width  <- 55
fig.height <- 20*3

SFig3  <- cowplot::plot_grid(SFig3A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             SFig3B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             SFig3C + draw_label("c", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Supp Fig 3.jpg', SFig3, height = fig.height, width = fig.width, limitsize = F)























## CTO vicinity
source('plot_prime5_bar_settings.R')
source('plot_prime5.R')

source('plot_prime3_bar_settings.R')
source('plot_prime3.R')


###

## coverage
cov.counts <- data.table(data.frame(win.cov.sum[,.(pos = round(mean(c(window_end, window_start)),0), count=sum_coverage),
                                                by=.(seqnames, sample, strand, group, hpi, Time, window_start )]))

source('plot_coverage_area_bin50_settings.R')
source('plot_prime3_area_bin50.R')


#### CAGE anyalysis
source('CAGE.R')


#### DBSCAN clustering ####

TR.prime3 <- prime3.counts
TR.prime3 <- TR.prime3[,.(seqnames, strand, sample, count, TR.prime3 = pos, hpi, cell_line, group)]

##### Cluster the counts of the prime3, in th hpi12 samples with count as weight
## OMIT CAGE FOR NOW
TR.prime3.sum <- unique(TR.prime3[cell_line == 'PK-15' & hpi == '12h' & !grepl('CAGE', group),
                                  .(sum_count = sum(count)),
                                  by=.(seqnames, strand, TR.prime3, hpi)])

fwrite(TR.prime3.sum, paste0(outdir, '/DBScan.prime3.cluster.data.tsv'), sep='\t')

TR.prime3.sum <- fread(paste0(outdir, '/DBScan.prime3.cluster.data.tsv'), na.strings = '')


TR.prime5 <- prime5.counts
TR.prime5 <- TR.prime5[,.(seqnames, strand, sample, count, TR.prime5 = pos, hpi, cell_line, group)]

##### Cluster the counts of the prime5, in th hpi12 samples with count as weight
## OMIT CAGE FOR NOW
TR.prime5.sum <- unique(TR.prime5[cell_line == 'PK-15' & hpi == '12h' & !grepl('CAGE', group),
                                  .(sum_count = sum(count)),
                                  by=.(seqnames, strand, TR.prime5, hpi)])

fwrite(TR.prime5.sum, paste0(outdir, '/DBScan.prime5.cluster.data.tsv'), sep='\t')

TR.prime5.sum <- fread(paste0(outdir, '/DBScan.prime5.cluster.data.tsv'), na.strings = '')



#### Initial Clustering
## Parameters
eps    <- 20
minPts <- 5

DT <- unique(TR.prime3.sum[,.(seqnames, strand, position=TR.prime3, count=sum_count)])

source('DBScan_primary.R')


#### Secondary Clustering
## Parameters
eps_sec    <- 10
minPts_sec <- 5

DT.clust     <- fread(paste0(outdir, "/DBSCAN_primary_clusters_", "eps.", eps, "_minPts.", minPts, '.tsv'), na.strings = '')
DT.clust.uni <- fread(paste0(outdir, "/DBSCAN_UNIQ_primary_clusters_", "eps.", eps, "_minPts.", minPts, '.tsv'), na.strings = '')
DT.clust.uni <- DT.clust.uni[order(cluster_center)]


source('DBScan_secondary.R')

DT.sec.clust <- fread(paste0(outdir, "/DBSCAN_secondary_clusters_", "eps.", eps_sec, "_minPts.", minPts_sec, '.tsv'), na.strings = '')



#### ####
##


##
#### Write outputs ####

save.image('PRV.LoRTIA.all.RData')

#### ####
##






stop()



### Make the plots
source('knit.docs.R')

## 1.) Coverage + TES, TSS windows
source('cluster.plot1_coverages.and.windows.R')

## 2.) Coverage + TR-reads ## -->> How to plot reads in each cell line?
#source('cluster.plot2_coverages.and.TR.reads.R')

## 3.) Coverage + TR-annot
source('cluster.plot3_coverages.and.TR.ref.R')

## 4.) TR-annot + TR-reads ## -->> How to plot reads in each cell line?
#source('cluster.plot3_coverages.and.TR.ref.R')

## 5.) TES, TSS windows + TR-annot ##
source('cluster.plot5_windows.and.TR.ref.R')


## OK

source('TSS_abund.LoRTIA.Rmd')
## KNIT !

#### ####
##













### host reads
viral_reads <- melt(gene.sample_count.sp, variable.name = 'sample', value.name = 'count')
viral_reads <- viral_reads[,.(viral_read_count = sum(count)), by=sample]

total_reads <- fread('../fastq_read_counts.tsv')

all_read_counts <- merge(total_reads, viral_reads, by='sample')
all_read_counts[,host_read_count := read_count - viral_read_count]

fwrite(all_read_counts, 'all_read_counts.tsv', sep = '\t')

