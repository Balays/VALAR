### Import gff_compare 'best' results
#best.merged.result_gff.compare <- fread(paste0(outdir, "/best.merged.result_gff.compare.tsv"))


seqs.to.keep <- genome

## GFF-COMPARE RESULTS TRANSCRIPT-LEVEL DISTANCES
TR.gff.compare.merged.TR <- best.merged.result_gff.compare
try({setnames(TR.gff.compare.merged.TR, old = c('distance_prime5.tr', 'distance_prime3.tr'), new = c('distance_prime5', 'distance_prime3'))})
TR.gff.compare.merged.TR$transcript_id <- as.character(TR.gff.compare.merged.TR$transcript_id)
##

### why these doesent match?
#sum(best.gff.compare.ref.TR.class.freq$N)
#length(unique(TR.reads$transcript_id))
#length(unique(TR.gene.ov.all$transcript_id))


### Include read counts 

TR.gff.compare.merged.TR.counts <- merge(TR.gff.compare.merged.TR, 
                                         TR.counts[,.(sample, transcript_id=as.character(TR_ID), count)], 
                                         by='transcript_id')

# since we used the melted count table and not the casted one, this will result in the same DT

TR.gff.compare.merged.TR.counts.gt <- melt(TR.gff.compare.merged.TR.counts, 
                                           id.vars = setdiff(names(TR.gff.compare.merged.TR.counts), metafilt$sample),
                                           variable.name = "sample", 
                                           value.name = "count")
#TR.gff.compare.merged.TR.counts %>% gather(sample, count, -c(1:20))
#setDT(TR.gff.compare.merged.TR.counts.gt)
#TR.gff.compare.merged.TR.counts.gt[,sample:=gsub('_pychopped', '', sample)]
TR.gff.compare.merged.TR.counts.gt <- merge(TR.gff.compare.merged.TR.counts.gt, metafilt, by='sample')


### Include average genome coverage
bym <- c('seqnames', metacols)
#TR.gff.compare.merged.TR.counts.gt[,rep := as.integer(rep), ]
TR.gff.compare.merged.TR.counts.gt <- merge(TR.gff.compare.merged.TR.counts.gt,
                                            norm_cov_summary,
                                            by=bym)
## normalize transcript counts to average genome coverage
TR.gff.compare.merged.TR.counts.gt[,norm_count := count / average_coverage]




#### Summarising results
## class code freq of best ref transcript per read
best.class_code.freq <- best.merged.result_gff.compare[
  , .N, by=.(seqnames, class_code)
]

## best class code per ref transcript
best.gff.compare.ref.TR.class.freq <- best.merged.result_gff.compare[
  , .N, by=.(seqnames, cmp_ref, class_code)
]

best.gff.compare.ref.TR.class.freq.sp <- dcast(best.gff.compare.ref.TR.class.freq,
                                               seqnames + cmp_ref ~ class_code, value.var = 'N'
)


## best class code per transcript
best.gff.compare.TR.class.freq <- best.merged.result_gff.compare[
  , .N, by=.(seqnames, transcript_id, class_code)
]

best.gff.compare.TR.class.freq.sp <- dcast(best.gff.compare.TR.class.freq,
                                      seqnames + transcript_id ~ class_code, value.var = 'N'
)



#### Write out
fwrite(TR.gff.compare.merged.TR.counts.gt, paste0(outdir, "/TR.gff.compare.merged.TR.counts.gt.tsv"), sep = '\t')
fwrite(best.gff.compare.ref.TR.class.freq, paste0(outdir, "/best.gff.compare.ref.TR.class.freq.tsv"), sep = '\t')
fwrite(TR.gff.compare.merged.TR, paste0(outdir, "/TR.gff.compare.merged.TR.tsv"), sep = '\t')
