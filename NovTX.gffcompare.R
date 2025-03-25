


#NOV.TXs.uni.sp <- fread(paste0(outdir, "/NOV.TXs.uni.sp.tsv"), na.strings = '')

#TR.gff.compare.merged.TR.counts.gt.NovTX <- fread(paste0(outdir, "/TR.gff.compare.merged.TR.counts.gt.NovTX.tsv.gz"), na.strings = '')



#### 
##
#### WF part 7. NovTX Analysis ####

### Run GFF-compare on each ref TR separately,
# Compare every TransFrag to Every Reference TX Iteratively
# then calculate distances
###
ref_mRNAs <- unique(viral.mrna[,transcript_id])
#
out_dir <-paste0(outdir, '/gffcompare_iterative_NovTX')
#
all.gff.compare.outfile <- "NovTXs.merged.result_gff.compare.tsv.gz"
#
query_gff <- paste0(outdir, '/NOV.TXs.gff2')
#
source('run_GFF.COMPARE.R')


## read in GFF-compare results
all.merged.result_gff.compare.NovTX  <- fread(paste0(outdir, "/", all.gff.compare.outfile), na.strings = '')



## Which classification result to use?
config$compare.transfrags.to <- 'closest.Ref'

if (config$compare.transfrags.to == 'canonic.TX') {
  
  
  ### 4.A Compare TransFrags to Canonical Transcripts (one for each gene)
  all.merged.result_gff.compare.canonic.NovTX  <- merge(all.merged.result_gff.compare.NovTX,
                                                  viral.ref[type == 'transcript' &
                                                              is.canonic == 'canonic', 
                                                            .(seqnames, strand, transcript_id)],
                                                  by.x = c('seqnames', 'strand', 'cmp_ref'), 
                                                  by.y = c('seqnames', 'strand', 'transcript_id') )
  
  all.merged.result_gff.compare.NovTX <- all.merged.result_gff.compare.canonic.NovTX
  
  
} else if (config$compare.transfrags.to == 'closest.Ref') {
  
  
  ### 4.B Classify TransFrags Compared to Closest Reference Transcripts (can be more than one for each gene)
  all.merged.result_gff.compare.NovTX <- all.merged.result_gff.compare.NovTX
  
}


### Analyse GFF-compare results:
## find the closest (canonic or all) ref-TR for each query,
## categorise non-equal matches
thresh.eq.prime5 <- config$thresh.eq.prime5
thresh.eq.prime3 <- config$thresh.eq.prime3
thresh.eq.junc   <- config$thresh.eq.junc
all.merged.result_gff.compare <- all.merged.result_gff.compare.NovTX
source('analyse_GFF.COMPARE.R')

#### Write out
fwrite(best.merged.result_gff.compare,     paste0(outdir, "/best.merged.result_gff.compare.NovTX.tsv.gz"),     sep = '\t')


## Import results back
best.merged.result_gff.compare <- fread(paste0(outdir, "/best.merged.result_gff.compare.NovTX.tsv.gz"), na.strings = '')
#TR.counts <- fread(paste0(outdir, "/TR.counts.tsv.gz"), na.strings = '')
#best.merged.result_gff.compare[,transcript_id := as.integer(transcript_id) ]


### Summarise results based on TransFrag class and metadata factors
TR.counts[,.(sample, transcript_id=as.character(TR_ID), count)]
TR.counts <- TransFrags.TSSr.filt.counts.gt[, .(count = sum(count)), by=.(NOV_TX_ID, sample)]
setnames(TR.counts, 'NOV_TX_ID', 'TR_ID')
source('summarise_GFF.COMPARE.R')
## Write out
fwrite(TR.gff.compare.merged.TR.counts.gt, paste0(outdir, "/TR.gff.compare.merged.TR.counts.gt.NovTX.tsv.gz"), sep = '\t')
fwrite(best.gff.compare.ref.TR.class.freq, paste0(outdir, "/best.gff.compare.ref.TR.class.freq.NovTX.tsv.gz"), sep = '\t')
fwrite(TR.gff.compare.merged.TR, paste0(outdir, "/TR.gff.compare.merged.TR.NovTX.tsv.gz"), sep = '\t')


## Import results back
TR.gff.compare.merged.TR.counts.gt <- fread(paste0(outdir, "/TR.gff.compare.merged.TR.counts.gt.NovTX.tsv.gz"), na.strings = '')


### Plot GFF-compare results
outdir <- paste0(outdir, '/NovTX') ;dir.create(outdir)
source('plot_GFF.COMPARE.R')



outdir  <- 'PRV-MDBIO-4cell'; try({ dir.create(outdir) })
